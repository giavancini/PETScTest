#include "SolidDomain.h"
#include <algorithm>

static std::vector<std::string> split(std::string& str)
{
	std::istringstream is(str);
	std::vector<std::string> values;
	std::string token;
	while (getline(is, token, ' '))
		values.push_back(token);
	return values;
}

// Public methods
SolidDomain::SolidDomain(Geometry* geometry, const int& index)
	: index_(index),
	  dimension_(2),
	  numberOfDOFs_(0),
	  numberOfBlockedDOFs_(0),
	  numberOfIsolatedDOFs_(0),
	  numberOfNodes_(0),
	  numberOfBlockedNodes_(0),
	  numberOfIsolatedNodes_(0),
	  parameters_(new AnalysisParameters()),
	  geometry_(geometry),
	  remesh_(nullptr),
	  perm_(nullptr)
{
	int fail = system("mkdir -p ./results");
	fail = system("rm ./results/*.vtu 2> /dev/null");
}

SolidDomain::~SolidDomain() {}

void SolidDomain::setNumberOfSteps(const int numberOfSteps)
{
	parameters_->setNumberOfSteps(numberOfSteps);
}

void SolidDomain::setMaxNonlinearIterations(const int maxNonlinearIterations)
{
	parameters_->setMaxNonlinearIterations(maxNonlinearIterations);
}

void SolidDomain::setNonlinearTolerance(const double nonlinearTolerance)
{
	parameters_->setNonlinearTolerance(nonlinearTolerance);
}

void SolidDomain::setDeltat(const double deltat)
{
	parameters_->setDeltat(deltat);
}

void SolidDomain::setGravity(const double gravity_x,
                			 const double gravity_y,
                			 const double gravity_z)
{
	parameters_->setGravity(gravity_x, gravity_y, gravity_z);
}

void SolidDomain::setSpectralRadius(const double rhoInf)
{
	parameters_->setSpectralRadius(rhoInf);
}

void SolidDomain::setGeneralizedAlphas(const double alphaM, const double alphaF)
{
	parameters_->setGeneralizedAlphas(alphaM, alphaF);
}

void SolidDomain::setAlpha(const double alpha)
{
	parameters_->setAlpha(alpha);
}

void SolidDomain::setMeshLength(const double h)
{
	parameters_->setMeshLength(h);
}

void SolidDomain::setInitialAccel(const bool initialAccel)
{
	parameters_->setInitialAccel(initialAccel);
}

void SolidDomain::setExportFrequency(const int& freq)
{
	parameters_->setExportFrequency(freq);
}

void SolidDomain::setLumpedMass(const bool& useLumpedMass)
{
    parameters_->setLumpedMass(useLumpedMass);
}

void SolidDomain::setReferenceConfiguration(const ReferenceConfiguration reference)
{
	for (Element*& el : elements_)
		el->setReferenceConfiguration(reference);
}

void SolidDomain::applyMaterial(const std::vector<Line*> lines, Material*& material)
{
	materials_.push_back(material);
	for (Line* line : lines)
		line->setMaterial(material);
}

void SolidDomain::applyMaterial(const std::vector<Surface*> surfaces, Material*& material)
{
	materials_.push_back(material);
	for (Surface* surface : surfaces)
		surface->setMaterial(material);
}

void SolidDomain::applyMaterial(const std::vector<Volume*> volumes, Material*& material)
{
	materials_.push_back(material);
	for (Volume* volume : volumes)
		volume->setMaterial(material);
}

void SolidDomain::generateMesh(const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm, std::string geofile, const std::string& gmshPath,
							   const bool& plotMesh, const bool& showInfo)
{	
	int rank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	PetscPrintf(PETSC_COMM_WORLD, "...Starting the Pre-processing Procedures...\n");
	auto start_timer = std::chrono::high_resolution_clock::now();

	std::pair<std::string, bool> pair; pair.second = false;
	if (rank == 0)
	{
		pair = createMesh(geometry_, elementType, algorithm, geofile, gmshPath, plotMesh, showInfo);

		for (int i = 1; i < size; i++)
		{
			MPI_Send(pair.first.c_str(), pair.first.length()+1, MPI_CHAR, i, 0, PETSC_COMM_WORLD);
			if (i == size-1)
			{
				MPI_Send(&pair.second, 1, MPI_C_BOOL, i, 0, PETSC_COMM_WORLD);
				pair.second = false;
			}
		}
	}
	else
	{
		MPI_Status status;
		MPI_Probe(0, 0, PETSC_COMM_WORLD, &status);
		int count;
		MPI_Get_count(&status, MPI_CHAR, &count);
		char buf[count+1];
		MPI_Recv(&buf, count+1, MPI_CHAR, 0, 0, PETSC_COMM_WORLD, &status);
		pair.first = buf;
		if (rank == size-1)
			MPI_Recv(&pair.second, 1, MPI_C_BOOL, 0, 0, PETSC_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Initial Mesh Generated. Elapsed time: %f\n", elapsed.count());

	readInput(pair.first, pair.second, elementType);
	transferGeometricBoundaryConditions();
	nodalNeighborSearch();
	elementalNeighborSearch();
	if (dimension_ == 2) remesh_ = new TriangularMesher;
	else if (dimension_ == 3) remesh_ = new TetrahedralMesher;

	//This should be later included inside the buildFreeSurface function. Maybe implement a new class called ModelPart, which stores all the mesh information
	for (Node* const& node : nodes_)
    {
		if (node->isIsolated())
		{
			node->setPreviouslyIsolated(true);
			node->setIsolated(false);
		}
		else
		{
			node->setPreviouslyIsolated(false);
		}
    }
	numberOfNodes_ = nodes_.size();
	numberOfIsolatedNodes_ = 0;
	numberOfIsolatedDOFs_ = 0;
	numberOfBlockedNodes_ = numberOfNodes_;
	numberOfBlockedDOFs_ = numberOfDOFs_;
	for (Node* const& node : nodes_)
	{
		unsigned int num_neighbor_nodes = node->getNeighborNodes().size();
		if (num_neighbor_nodes == 1)
		{
			node->setIsolated(true);
			numberOfBlockedNodes_ -= 1;
			numberOfBlockedDOFs_ -= node->getNumberOfDegreesOfFreedom();
			if (!node->isConstrained())
			{
				numberOfIsolatedNodes_++;
				numberOfIsolatedDOFs_ += node->getNumberOfDegreesOfFreedom();
			}
		}
	}

	reorderDOFs();
	domainDecomposition();
	PetscPrintf(PETSC_COMM_WORLD, "...Ending the Pre-processing Procedures...\n");
}

void SolidDomain::solveTransientProblem()
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	setReferenceConfiguration(ReferenceConfiguration::INITIAL);
	if (parameters_->getInitialAccel()) computeInitialAccel();
	
	// Petsc variables
	Mat               tangent;
    Vec               rhs, solution;
    KSP               ksp;
    PC                pc;
	PetscLogDouble    bytes = 0.0;
	PetscBool		  isbjacobi;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	PetscMemoryGetCurrentUsage(&bytes);
	PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes/(1024*1024));

	//Initial position norm
    double initialPositionNorm = getInitialPositionNorm();

	//Array of constrained degrees of freedom
	int numberOfConstrainedDOFs;
	int* constrainedDOFs;
	getConstrainedDOFs(numberOfConstrainedDOFs, constrainedDOFs);

	//Array of external forces
	int ndofsForces; std::vector<DegreeOfFreedom*> dofsForces; double* externalForces;
	getExternalForces(ndofsForces, dofsForces, externalForces);

	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp, KSPFGMRES);
	KSPSetTolerances(ksp, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100);
	KSPGMRESSetRestart(ksp, 100);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCBJACOBI);
	// KSPSetType(ksp, KSPPREONLY);
	// KSPGetPC(ksp, &pc);
	// PCSetType(pc, PCLU);

	createSystemMatrix(tangent);

	//Create PETSc vectors
	VecCreate(PETSC_COMM_WORLD, &rhs);
	VecSetSizes(rhs, PETSC_DECIDE, numberOfBlockedDOFs_);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs, &solution);

	const int numberOfSteps = parameters_->getNumberOfSteps();
	const int maxNonlinearIterations = parameters_->getMaxNonlinearIterations();
	const double nonlinearTolerance = parameters_->getNonlinearTolerance();

	if (rank == 0) exportToParaview(0);
	
	for (int timeStep = 0; timeStep < numberOfSteps; timeStep++)
	{
		PetscPrintf(PETSC_COMM_WORLD, "\n----------------------- TIME STEP = %d, time = %f  -----------------------\n\n", timeStep+1, (double)(timeStep+1)*parameters_->getDeltat());
		parameters_->setCurrentTime(parameters_->getDeltat() * (double)(timeStep+1));
		double modelVolume = 0.0;
    	for (Element* const& el : elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		parameters_->setModelVolume(modelVolume);

		setPastVariables();
		computeCurrentVariables();
		computeIntermediateVariables();
		double positionNorm, pressureNorm;
		for (int iteration = 0; (iteration < maxNonlinearIterations); iteration++)
		{
			applyNeummanConditions(rhs, ndofsForces, dofsForces, externalForces, 1.0);
			assembleTransientLinearSystem(tangent, rhs);
			MatZeroRowsColumns(tangent, numberOfConstrainedDOFs, constrainedDOFs, 1.0, solution, rhs);
			MatView(tangent, PETSC_VIEWER_DRAW_WORLD);
			solveLinearSystem(ksp, tangent, rhs, solution);
			updateVariables(solution, positionNorm, pressureNorm);
			computeCurrentVariables();
			computeIntermediateVariables();

			PetscMemoryGetCurrentUsage(&bytes);
			PetscPrintf(PETSC_COMM_WORLD, "Newton iteration: %d - L2 Position Norm: %E - L2 Pressure Norm: %E\nMemory used by each processor: %f Mb\n", 
						iteration, positionNorm / initialPositionNorm, pressureNorm, bytes/(1024*1024));
	
			MatZeroEntries(tangent);
			VecZeroEntries(rhs);

			if (positionNorm / initialPositionNorm <= nonlinearTolerance)
				break;
		}

		//export results to paraview
		if (rank == 0 && ((timeStep + 1) % parameters_->getExportFrequency() == 0))
		{
			computeCauchyStress();
			exportToParaview(timeStep+1);
		}
	}
	delete[] constrainedDOFs;
	delete[] externalForces;
	KSPDestroy(&ksp);
    VecDestroy(&rhs);
    VecDestroy(&solution);
	MatDestroy(&tangent);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Solid Analysis Done. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::setInitialVelocityX(std::function<double(double, double, double)> function)
{
	for (Node*& node : nodes_)
	{
		double x = node->getDegreeOfFreedom(0)->getInitialValue();
		double y = node->getDegreeOfFreedom(1)->getInitialValue();
		double z = node->getDegreeOfFreedom(2)->getInitialValue();

		double velocity = function(x, y, z);
		node->getDegreeOfFreedom(0)->setInitialFirstTimeDerivative(velocity);
		node->getDegreeOfFreedom(0)->setCurrentFirstTimeDerivative(velocity);
	}
}

void SolidDomain::setInitialVelocityY(std::function<double(double, double, double)> function)
{
	for (Node*& node : nodes_)
	{
		double x = node->getDegreeOfFreedom(0)->getInitialValue();
		double y = node->getDegreeOfFreedom(1)->getInitialValue();
		double z = node->getDegreeOfFreedom(2)->getInitialValue();

		double velocity = function(x, y, z);
		node->getDegreeOfFreedom(1)->setInitialFirstTimeDerivative(velocity);
		node->getDegreeOfFreedom(1)->setCurrentFirstTimeDerivative(velocity);
	}
}

void SolidDomain::setInitialVelocityZ(std::function<double(double, double, double)> function)
{
	for (Node*& node : nodes_)
	{
		double x = node->getDegreeOfFreedom(0)->getInitialValue();
		double y = node->getDegreeOfFreedom(1)->getInitialValue();
		double z = node->getDegreeOfFreedom(2)->getInitialValue();

		double velocity = function(x, y, z);
		node->getDegreeOfFreedom(2)->setInitialFirstTimeDerivative(velocity);
		node->getDegreeOfFreedom(2)->setCurrentFirstTimeDerivative(velocity);
	}
}

// Private methods
Node* SolidDomain::getNode(const int& index)
{
	return nodes_[index];
}

Element* SolidDomain::getElement(const int& index)
{
	return elements_[index];
}

Material* SolidDomain::getMaterial(const int& index)
{
	return materials_[index];
}

const std::vector<DirichletBoundaryCondition*>& SolidDomain::getDirichletBoundaryConditions()
{
	return dirichletBoundaryConditions_;
}

const std::vector<NeumannBoundaryCondition*>& SolidDomain::getNeumannBoundaryConditions()
{
	return neumannBoundaryConditions_;
}

double SolidDomain::getInitialPositionNorm() const
{
    double initialNorm = 0.0;
    for (auto& node : nodes_)
        for (int i = 0; i < dimension_; i++)
		{
			double value = node->getDegreeOfFreedom(i)->getInitialValue();
			initialNorm += value * value;
		}
    return sqrt(initialNorm);
}

void SolidDomain::setPastVariables()
{
	for (Node*& node : nodes_)
	{
		for (int i = 0; i < dimension_; i++)
		{
			DegreeOfFreedom* dof = node->getDegreeOfFreedom(i);
			dof->setPastValue(dof->getCurrentValue());
			dof->setPastFirstTimeDerivative(dof->getCurrentFirstTimeDerivative());
			dof->setPastSecondTimeDerivative(dof->getCurrentSecondTimeDerivative());
		}
	}
}

void SolidDomain::computeCurrentVariables()
{
	double gamma = parameters_->getGamma();
	double beta = parameters_->getBeta();
	double deltat = parameters_->getDeltat();

	for (Node*& node : nodes_)
	{
		for (int i = 0; i < dimension_; i++)
		{
			double vel, accel;
			DegreeOfFreedom* dof = node->getDegreeOfFreedom(i);
			accel = (dof->getCurrentValue() - dof->getPastValue()) / (beta * deltat * deltat) -
					dof->getPastFirstTimeDerivative() / (beta * deltat) - dof->getPastSecondTimeDerivative() * (0.5 / beta - 1.0);
			dof->setCurrentSecondTimeDerivative(accel);
			vel = gamma * deltat * dof->getCurrentSecondTimeDerivative() + dof->getPastFirstTimeDerivative() + 
				  deltat * (1.0 - gamma) * dof->getPastSecondTimeDerivative();
			dof->setCurrentFirstTimeDerivative(vel);
		}
	}
}

void SolidDomain::computeIntermediateVariables()
{
	//Computing alpha-generalized and Newmark parameters
	double alphaM = parameters_->getAlphaM();
	double alphaF = parameters_->getAlphaF();

	for (Node *node : nodes_)
	{
		for (int i = 0; i < dimension_; i++)
		{
			double coordinate, vel, accel;
			DegreeOfFreedom* dof = node->getDegreeOfFreedom(i);
			coordinate = dof->getPastValue() + alphaF * (dof->getCurrentValue() - dof->getPastValue());
			dof->setIntermediateValue(coordinate);
			vel = dof->getPastFirstTimeDerivative() + alphaF * (dof->getCurrentFirstTimeDerivative() - dof->getPastFirstTimeDerivative());
			dof->setIntermediateFirstTimeDerivative(vel);
			accel = dof->getPastSecondTimeDerivative() + alphaM * (dof->getCurrentSecondTimeDerivative() - dof->getPastSecondTimeDerivative());
			dof->setIntermediateSecondTimeDerivative(accel);
		}
	}
}

void SolidDomain::getConstrainedDOFs(int& ndofs, int*& constrainedDOFs)
{
	ndofs = dirichletBoundaryConditions_.size();
	constrainedDOFs = new int[ndofs];
	
	for (int i = 0; i < ndofs; i++)
		constrainedDOFs[i] = dirichletBoundaryConditions_[i]->getDegreeOfFreedom()->getIndex();
}

void SolidDomain::getExternalForces(int& ndofs, std::vector<DegreeOfFreedom*>& dofs, double*& externalForces)
{
	ndofs = 0;
	for (NeumannBoundaryCondition* const& nbc : neumannBoundaryConditions_)
	{
		ndofs += nbc->getNumberOfDOFs();
	}
	dofs.reserve(ndofs);
	externalForces = new double[ndofs];
	int aux = -1;
	for (NeumannBoundaryCondition* const& nbc : neumannBoundaryConditions_)
	{
		int ndof = nbc->getNumberOfDOFs();
		std::vector<DegreeOfFreedom*> nbc_dofs; 
		double *val;
		nbc->getNodalForce(dimension_, nbc_dofs, val);
		for (int i = 0; i < ndof; i++)
		{
			dofs.push_back(nbc_dofs[i]);
			externalForces[++aux] = val[i];
		}
		delete[] val;
	}
}

void SolidDomain::assembleTransientLinearSystem(Mat& mat, Vec& vec)
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	for (Element* const& el : elements_)
	{
		if (el->getRank() == rank && el->isActive())
		{
			int ndofsPosition, ndofsPressure;
			int* indexes;
			double *rhsValues, *hessianValues;
			el->elementContributions(ndofsPosition, ndofsPressure, indexes, rhsValues, hessianValues);
			int ndofs = ndofsPosition + ndofsPressure;
			// dispersing element rhs contribution into global rhs vector
			VecSetValues(vec, ndofs, indexes, rhsValues, ADD_VALUES);

			// dispersing element matrices contribution into global tangent matrix
			for (int i = 0; i < ndofsPosition; i++)
			{
				for (int j = 0; j < ndofsPosition; j++)
				{
					double value = hessianValues[ndofs * i + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[i], 1, &indexes[j], &value, ADD_VALUES);
				}
				for (int j = 0; j < ndofsPressure; j++)
				{
					double value = hessianValues[ndofsPosition + ndofs*i + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[i], 1, &indexes[ndofsPosition + j], &value, ADD_VALUES);
				}
			}
			
			for (int i = 0; i < ndofsPressure; i++)
			{
				for (int j = 0; j < ndofsPosition; j++)
				{
					double value = hessianValues[ndofs*(ndofsPosition+i) + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[ndofsPosition + i], 1, &indexes[j], &value, ADD_VALUES);
				}
				for (int j = 0; j < ndofsPressure; j++)
				{
					double value = hessianValues[ndofs*(ndofsPosition+i)+ndofsPosition + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[ndofsPosition + i], 1, &indexes[ndofsPosition + j], &value, ADD_VALUES);
				}
			}
			delete[] indexes;
			delete[] rhsValues;
			delete[] hessianValues;
		}
	}

	//Assemble matrices and vectors
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);		
    VecAssemblyBegin(vec); 
    VecAssemblyEnd(vec);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Assemble Linear System. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::applyNeummanConditions(Vec& vec, int& ndofs, const std::vector<DegreeOfFreedom*>& dofsForces, double*& externalForces, const double& loadFactor)
{
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	if (rank == 0)
	{
		if (ndofs > 0)
		{
			int* indexes = new int[ndofs];
			double* values = new double[ndofs];
			for (int i = 0; i < ndofs; i++)
			{
				values[i] = externalForces[i] * loadFactor;
				indexes[i] = dofsForces[i]->getIndex();
			}
			VecSetValues(vec, ndofs, indexes, values, ADD_VALUES);
			delete[] indexes;
			delete[] values;
		}
	}
}

void SolidDomain::solveLinearSystem(KSP& ksp, Mat& mat, Vec& rhs, Vec& solution)
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	KSPSetOperators(ksp, mat, mat);
	PC pc;
	KSPGetPC(ksp, &pc);
	PetscBool isbjacobi;
	PetscObjectTypeCompare((PetscObject)pc, PCBJACOBI, &isbjacobi);
	
	if (isbjacobi)
	{
		PetscInt nlocal;
		KSP *subksp;
		PC subpc;
		KSPSetUp(ksp);
		PCBJacobiGetSubKSP(pc, &nlocal, NULL, &subksp);
		for (int i = 0; i < nlocal; i++)
		{
			KSPGetPC(subksp[i], &subpc);
			PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
		}
	}
	else
	{
		PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
	}

	KSPSetFromOptions(ksp);
	KSPSolve(ksp, rhs, solution);
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	if (reason == -3 || reason > 0)
	{
		PetscInt nit;
		KSPGetIterationNumber(ksp, &nit);
		PetscPrintf(PETSC_COMM_WORLD, "Converged with %d iterations. Elapsed time: %f\n", nit, elapsed.count());
	}
	else
	{
		PetscPrintf(PETSC_COMM_WORLD, "Diverged, please check PETSc manual. Reason: %d\n", reason);
		exit(EXIT_FAILURE);
	}
}

void SolidDomain::updateVariables(Vec& solution, double& positionNorm, double& pressureNorm)
{
	Vec All;
	VecScatter ctx;

	//Gathers the solution vector to the master process
	VecScatterCreateToAll(solution, &ctx, &All);
	VecScatterBegin(ctx, solution, All, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(ctx, solution, All, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterDestroy(&ctx);

	//Updates nodal variables
	positionNorm = 0.0;
	pressureNorm = 0.0;

	for (Node* const& node : nodes_)
	{
		int ndof = node->getNumberOfDegreesOfFreedom();
		if (!node->isIsolated())
		{
			for (int j = 0; j < dimension_; j++)
			{
				DegreeOfFreedom* dof = node->getDegreeOfFreedom(j);
				int index = dof->getIndex();
				double val;
				VecGetValues(All, 1, &index, &val);
				positionNorm += val * val;
				dof->incrementCurrentValue(val);
			}
			for (int j = dimension_; j < ndof; j++)
			{
				DegreeOfFreedom* dof = node->getDegreeOfFreedom(j);
				int index = dof->getIndex();
				double val;
				VecGetValues(All, 1, &index, &val);
				pressureNorm += val * val;
				dof->incrementCurrentValue(val);
			}
		}
	}
	
	positionNorm = sqrt(positionNorm);
	pressureNorm = sqrt(pressureNorm);
	VecDestroy(&All);
}

void SolidDomain::computeCauchyStress()
{
	const unsigned int numberOfNodes = nodes_.size();
	int* nodalContribution = new int[numberOfNodes];
	for (int i = 0; i < numberOfNodes; i++)
		nodalContribution[i] = 0;

	const int numberOfStressComponents = dimension_ * (dimension_ + 1) / 2;

	for (Node*& node : nodes_)
		node->clearCauchyStress(numberOfStressComponents);

	for (Element*& el : elements_)
	{
		double** cauchyStress;
		el->getCauchyStress(cauchyStress);
		const std::vector<Node*>& nodes = el->getNodes();
		int nElemNodes = nodes.size();
		for (int i = 0; i < nElemNodes; i++)
		{
			int index = nodes[i]->getIndex();
			nodes_[index]->incrementCauchyStress(numberOfStressComponents, cauchyStress[i]);
			nodalContribution[index]++;
			delete[] cauchyStress[i];
		}
		delete[] cauchyStress;
	}
	
	for (int i = 0; i < numberOfNodes; i++)
	{
		double* cauchy = nodes_[i]->getCauchyStress();
		for (int j = 0; j < numberOfStressComponents; j++)
			cauchy[j] /= nodalContribution[i];
	}
	delete[] nodalContribution;
}

void SolidDomain::computeInitialAccel()
{

}

void SolidDomain::exportToParaview(const int& step)
{
	std::stringstream text;
	text << "results/" << "solidOutput" << step << ".vtu";
	std::ofstream file(text.str());

	unsigned int numberOfElements = 0;
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
			if (elem->getPlot())
				numberOfElements++;
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
			if (elem->getPlot())
				numberOfElements++;
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
			if (elem->getPlot())
				numberOfElements++;
	}

	bool mixed = false;
	if (materials_[0]->getType() == MaterialType::ELASTIC_INCOMPRESSIBLE_SOLID ||
		materials_[0]->getType() == MaterialType::NEWTONIAN_INCOMPRESSIBLE_FLUID ) mixed = true;

	//header
	file << "<?xml version=\"1.0\"?>" << "\n"
         << "<VTKFile type=\"UnstructuredGrid\">" << "\n"
		 << "  <UnstructuredGrid>" << "\n"
         << "  <Piece NumberOfPoints=\"" << nodes_.size()
         << "\"  NumberOfCells=\"" << numberOfElements
         << "\">" << "\n";

	//nodal coordinates
	file << "    <Points>" << "\n"
         << "      <DataArray type=\"Float64\" "
         << "NumberOfComponents=\"" << 3 << "\" format=\"ascii\">" << "\n";
    for (Node*& n : nodes_)
	{
		std::vector<DegreeOfFreedom*> dofs = n->getDegreesOfFreedom();
		file << dofs[0]->getCurrentValue() << " " << dofs[1]->getCurrentValue() << " " << ((dimension_ == 3)? dofs[2]->getCurrentValue() : 0.0) << "\n";
	}
    file << "      </DataArray>" << "\n"
         << "    </Points>" << "\n";

	//element connectivity
	file << "    <Cells>" << "\n"
         << "      <DataArray type=\"Int32\" "
         << "Name=\"connectivity\" format=\"ascii\">" << "\n";
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				const std::vector<int>& vtkConnectivity = parametricElement->getVTKConnectivity();
				const std::vector<Node*>& nodes = elem->getNodes();
				const unsigned int numberOfNodes = nodes.size();
				for (unsigned int i = 0; i < numberOfNodes; i++)
				{
					int nodeIndex = vtkConnectivity[i];
					file << nodes[nodeIndex]->getIndex() << " ";
				}
				file << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				const std::vector<int>& vtkConnectivity = parametricElement->getVTKConnectivity();
				const std::vector<Node*>& nodes = elem->getNodes();
				const unsigned int numberOfNodes = nodes.size();
				for (int i = 0; i < numberOfNodes; i++)
				{
					int nodeIndex = vtkConnectivity[i];
					file << nodes[nodeIndex]->getIndex() << " ";
				}
				file << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				const std::vector<int>& vtkConnectivity = parametricElement->getVTKConnectivity();
				const std::vector<Node*>& nodes = elem->getNodes();
				const unsigned int numberOfNodes = nodes.size();
				for (int i = 0; i < numberOfNodes; i++)
				{
					int nodeIndex = vtkConnectivity[i];
					file << nodes[nodeIndex]->getIndex() << " ";
				}
				file << "\n";
			}
		}
	}
	file << "      </DataArray>" << "\n";

	//offsets
	file << "      <DataArray type=\"Int32\""
		 << " Name=\"offsets\" format=\"ascii\">" << "\n";
	int aux = 0;
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				aux += elem->getNumberOfNodes();
				file << aux << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				aux += elem->getNumberOfNodes();
				file << aux << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				aux += elem->getNumberOfNodes();
				file << aux << "\n";
			}
		}
	}
	file << "      </DataArray>" << "\n";

	//elements type
	file << "      <DataArray type=\"UInt8\" Name=\"types\" "
		 << "format=\"ascii\">" << "\n";
	for (const auto& pair : geometry_->getLines())
	{
		Line* line = pair.second;
		const std::vector<BaseLineElement*>& elements = line->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				VTKCellType vtkType = parametricElement->getVTKCellType();
				file << vtkType << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getSurfaces())
	{
		Surface* surface = pair.second;
		const std::vector<BaseSurfaceElement*>& elements = surface->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				VTKCellType vtkType = parametricElement->getVTKCellType();
				file << vtkType << "\n";
			}
		}
	}
	for (const auto& pair : geometry_->getVolumes())
	{
		Volume* volume = pair.second;
		const std::vector<BaseVolumeElement*>& elements = volume->getBaseElements();
		for (auto& elem : elements)
		{
			if (elem->getPlot())
			{
				ParametricElement* parametricElement = elem->getParametricElement();
				VTKCellType vtkType = parametricElement->getVTKCellType();
				file << vtkType << "\n";
			}
		}
	}
	file << "      </DataArray>" << "\n"
		 << "    </Cells>" << "\n";
	
	//nodal results
	file << "    <PointData>" <<"\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_ <<"\" "
		 << "Name=\"Displacement\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		const std::vector<DegreeOfFreedom*>& dofs = n->getDegreesOfFreedom();
		for (int i = 0; i < dimension_; i++)
			file << dofs[i]->getCurrentValue() - dofs[i]->getInitialValue() << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_ << "\" "
		 << "Name=\"Velocity\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		const std::vector<DegreeOfFreedom*>& dofs = n->getDegreesOfFreedom();
		for (int i = 0; i < dimension_; i++)
			file << dofs[i]->getCurrentFirstTimeDerivative() << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_ << "\" "
		 << "Name=\"Acceleration\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		const std::vector<DegreeOfFreedom*>& dofs = n->getDegreesOfFreedom();
		for (int i = 0; i < dimension_; i++)
			file << dofs[i]->getCurrentSecondTimeDerivative() << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << dimension_*(dimension_+1)/2 << "\" "
		 << "Name=\"CauchyStress\" format=\"ascii\">" << "\n";
	for (Node*& n: nodes_)
	{
		double* cauchyStress = n->getCauchyStress();
		for (int i = 0; i < dimension_*(dimension_+1)/2; i++)
			file << cauchyStress[i] << " ";
		file << "\n";
	}
	file << "      </DataArray> " << "\n";
	if (mixed)
	{
		file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << 1 << "\" "
			<< "Name=\"Pressure\" format=\"ascii\">" << "\n";
		for (Node*& n: nodes_)
		{
			file << n->getDegreeOfFreedom(dimension_)->getCurrentValue() << "\n";
		}
		file << "      </DataArray> " << "\n";
	}
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << 1 <<"\" "
		 << "Name=\"PermutedIndex\" format=\"ascii\">" << "\n";
	for (Node*& node : nodes_)
	{
		file << node->getPermutedIndex() << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "    </PointData>" << "\n";

	//elemental results
	file << "    <CellData>" << "\n";
	file << "      <DataArray type=\"Float64\" NumberOfComponents=\"" << 1 <<"\" "
		 << "Name=\"Rank\" format=\"ascii\">" << "\n";
	for (Element* const& el : elements_)
	{
		file << el->getRank() << "\n";
	}
	file << "      </DataArray> " << "\n";
	file << "    </CellData>" << "\n";

	//footnote
	file << "  </Piece>" << "\n"
		<< "  </UnstructuredGrid>" << "\n"
		<< "</VTKFile>" << "\n";
	file.close();

	// std::stringstream text2;
	// text2 << "error.dat";
	// std::ofstream file2;
	// if (timestep == 0) file2.open(text2.str());
	// else file2.open(text2.str(), std::ofstream::out | std::ofstream::app);
	// file2.precision(16);

	// file2 << timestep << "\n";
	// for (Node* n : nodes_)
	// {
	// 	file2 << n->getCurrentCoordinate()(0) << " " << n->getCurrentCoordinate()(1) << "\n";
	// }

	// std::stringstream text2;
	// text2 << "error.dat";
	// std::ifstream file2;
	// file2.open(text2.str());
	// std::string line;
	// vector<double> reference_position(2 * nodes_.size());
	// for (int i = 0; i < 4000; i++)
	// {
	// 	std::getline(file2, line);
	// 	if (i == 16 * timestep + 15)
	// 	{
	// 		for (int j = 0; j < nodes_.size(); j++)
	// 		{
	// 			file2 >> reference_position(2 * j) >> reference_position(2 * j + 1);
	// 			std::getline(file2, line);
	// 		}
	// 		break;
	// 	}
	// 	else
	// 	{
	// 		for (int j = 0; j < nodes_.size(); j++)
	// 		{
	// 			std::getline(file2, line);
	// 		}
	// 	}
	// }

	// // Local error
	// std::stringstream text3;
	// text3 << "localerror-position" << deltat_ << ".dat";
	// std::ofstream file3;
	// file3.open(text3.str(), std::ofstream::out | std::ofstream::app);
	// double error_position = sqrt((nodes_[1]->getCurrentCoordinate()(0) - reference_position(2)) * (nodes_[1]->getCurrentCoordinate()(0) - reference_position(2)));

	// file3 << deltat_ * (timestep + 1.0) << " " << error_position << "\n";

	// // Norm L2
	// std::stringstream text4;
	// text4 << "globalerror-position" << deltat_ << ".dat";

	// std::ofstream file4;
	// file4.open(text4.str(), std::ofstream::out | std::ofstream::app);
	// error_position = 0.0;
	// for (Element* el : elements_)
	// {
	// 	vector<double> local_reference_position(2*el->getNodes().size());
	// 	for (int i = 0; i < el->getNodes().size(); i++)
	// 	{
	// 		int index = el->getNode(i)->getIndex();
	// 		local_reference_position(2*i) = reference_position(2*index);
	// 		local_reference_position(2*i+1) = reference_position(2*index+1);
	// 	}
	// 	error_position += el->domainIntegration(local_reference_position);
	// }
	// error_position = sqrt(error_position);
	// file4 << deltat_ * (timestep + 1.0) << " " << error_position << "\n";	
}

void SolidDomain::readInput(const std::string& inputFile, const bool& deleteFiles, const PartitionOfUnity elementType)
{
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	auto start_timer = std::chrono::high_resolution_clock::now();

	//defining the maps that are used to store the elements information
	const std::unordered_map<int, ParametricElement*> gmshElements = 
	{
		//lines
		{1, &ParametricLineElement::L2}, {8, &ParametricLineElement::L3}, {26, &ParametricLineElement::L4},
		//triangles
		{2, &ParametricSurfaceElement::T3}, {9, &ParametricSurfaceElement::T6}, {21, &ParametricSurfaceElement::T10},
		//quadrilaterals
		{3, &ParametricSurfaceElement::Q4}, {10, &ParametricSurfaceElement::Q9}, {36, &ParametricSurfaceElement::Q16},
		//tetrahedrons
		{4, &ParametricVolumeElement::TET4}, {11, &ParametricVolumeElement::TET10}, {29, &ParametricVolumeElement::TET20},
		//hexahedrons
		{5, &ParametricVolumeElement::HEX8}, {12, &ParametricVolumeElement::HEX27}, {92, &ParametricVolumeElement::HEX64}
	};

	const std::unordered_map<PartitionOfUnity, unsigned int> dimension =
	{
		//lines
		{L2, 3}, {L3, 3}, {L4, 3},
		//triangles
		{T3, 2}, {T6, 2}, {T10, 2},
		//quadrilaterals
		{Q4, 2}, {Q9, 2}, {Q16, 2},
		//tetrahedrons
		{TET4, 3}, {TET10, 3}, {TET20, 3}
	};

	dimension_ = dimension.at(elementType);
	parameters_->setDimension(dimension_);

	bool mixed = false;
	if (materials_[0]->getType() == MaterialType::ELASTIC_INCOMPRESSIBLE_SOLID ||
		materials_[0]->getType() == MaterialType::NEWTONIAN_INCOMPRESSIBLE_FLUID ) mixed = true;

	int ndofs_per_node = dimension_;
	if (mixed) ndofs_per_node++;

	//opening the .msh file
	std::ifstream file(inputFile);
	std::string line;
	std::getline(file, line); std::getline(file, line); std::getline(file, line); std::getline(file, line);

	//reading physical entities
	unsigned int nEntities;
	file >> nEntities;
	std::getline(file, line);
	std::unordered_map<int, std::string> physicalEntities;
	physicalEntities.reserve(nEntities);
	for (unsigned int i = 0; i < nEntities; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		unsigned int index;
		std::istringstream(tokens[1]) >> index;
		physicalEntities[index] = tokens[2].substr(1, tokens[2].size() - 2);
	}
	std::getline(file, line);
	std::getline(file, line);

	//reading nodes
	file >> numberOfNodes_;
	nodes_.reserve(numberOfNodes_);
	std::getline(file, line);
	for (unsigned int i = 0; i < numberOfNodes_; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		std::vector<DegreeOfFreedom *> degreesOfFreedom;
		degreesOfFreedom.reserve(dimension_);
		for (unsigned int j = 0; j < dimension_; j++)
		{
			double coord;
			std::istringstream(tokens[j + 1]) >> coord;
			degreesOfFreedom.emplace_back(new DegreeOfFreedom(DOFType::POSITION, coord));
			numberOfDOFs_++;
		}
		nodes_.emplace_back(new Node(i, degreesOfFreedom));
	}
	std::getline(file, line);
	std::getline(file, line);

	// Adding pressure degrees of freedom for mixed formulation
	if (mixed)
		for (Node *&node : nodes_)
		{
			DegreeOfFreedom *dof = new DegreeOfFreedom(DOFType::PRESSURE, 0.0);
			node->addDegreeOfFreedom(dof);
			numberOfDOFs_++;
		}

	//Pre allocating elements
	unsigned int nElements;
	file >> nElements;
	std::getline(file, line);
	const std::streampos position_to_rewind = file.tellg();
	unsigned int nLineElements = 0;
	unsigned int nSurfaceElements = 0;
	unsigned int nVolumeElements = 0;
	for (int i = 0; i < nElements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		
		std::string name = physicalEntities[values[3]];
		if (name[0] == 'l')
		{
			ParametricElement* type = gmshElements.at(values[1]);
			int numberOfNodes = type->getNumberOfNodes();
			nLineElements++;
			Line *object = geometry_->getLine(name);
			object->numberOfElements_++;
			object->numberOfNodes_ += numberOfNodes;
		}
		else if (name[0] == 's')
		{
			ParametricElement* type = gmshElements.at(values[1]);
			int numberOfNodes = type->getNumberOfNodes();
			nSurfaceElements++;
			Surface *object = geometry_->getSurface(name);
			object->numberOfElements_++;
			object->numberOfNodes_ += numberOfNodes;
		}
		else if (name[0] == 'v')
		{
			ParametricElement* type = gmshElements.at(values[1]);
			int numberOfNodes = type->getNumberOfNodes();
			nVolumeElements++;
			Volume *object = geometry_->getVolume(name);
			object->numberOfElements_++;
			object->numberOfNodes_ += numberOfNodes;
		}
	}

	for (auto& pair : geometry_->lines_)
	{
		auto& line = pair.second;
		line->baseElements_.reserve(line->numberOfElements_);
		line->elements_.reserve(line->numberOfElements_);
		line->nodes_.reserve(line->numberOfNodes_);
	}
	for (auto& pair : geometry_->surfaces_)
	{
		auto& surface = pair.second;
		surface->baseElements_.reserve(surface->numberOfElements_);
		surface->elements_.reserve(surface->numberOfElements_);
		surface->nodes_.reserve(surface->numberOfNodes_);
	}
	for (auto& pair : geometry_->volumes_)
	{
		auto& volume = pair.second;
		volume->baseElements_.reserve(volume->numberOfElements_);
		volume->elements_.reserve(volume->numberOfElements_);
		volume->nodes_.reserve(volume->numberOfNodes_);
	}

	if (nVolumeElements)
		elements_.reserve(nVolumeElements);
	else if (nSurfaceElements)
		elements_.reserve(nSurfaceElements);
	else
		elements_.reserve(nLineElements);

	//reading elements
	file.seekg(position_to_rewind);
	unsigned int nPoints = nElements - nLineElements - nSurfaceElements - nVolumeElements;
	for (int i = 0; i < nPoints; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		Node* node = nodes_[values[5] - 1];
		std::string name = physicalEntities[values[3]];
		Point* object = geometry_->getPoint(name);
		object->addNode(node);
	}
	for (int i = 0; i < nLineElements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		ParametricLineElement* type = static_cast<ParametricLineElement*>(gmshElements.at(values[1]));
		int numberOfNodes = type->getNumberOfNodes();
		std::vector<Node *> elementNodes;
		elementNodes.reserve(numberOfNodes);
		for (size_t j = 5; j < values.size(); j++)
			elementNodes.push_back(nodes_[values[j] - 1]);
		std::string name = physicalEntities[values[3]];
		Line* object = geometry_->getLine(name);
		BaseLineElement *base_elem = new BaseLineElement(i, *type, elementNodes);
		object->addBaseElement(base_elem);
		object->addNodes(elementNodes);
		switch (elementType)
		{
		case L2:
		case L3:
		case L4:
			base_elem->setPlot(true);
			Material *mat = object->getMaterial();
			int ndofs = ndofs_per_node * elementNodes.size();
			std::vector<DegreeOfFreedom*> dofs;
			dofs.reserve(ndofs);
			for (auto& node : elementNodes)
				for (int j = 0; j < dimension_; j++)
					dofs.push_back(node->getDegreeOfFreedom(j));
			if (mixed)
				for (auto& node : elementNodes)
					dofs.push_back(node->getDegreeOfFreedom(dimension_));
			Element* element = new LineElement(i, dofs, mat, base_elem, parameters_);
			elements_.push_back(element);
			object->addElement(element);
			break;
		}
	}
	for (int i = 0; i < nSurfaceElements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		ParametricSurfaceElement* type = static_cast<ParametricSurfaceElement*>(gmshElements.at(values[1]));
		int numberOfNodes = type->getNumberOfNodes();
		std::vector<Node*> elementNodes;
		elementNodes.reserve(numberOfNodes);
		for (size_t j = 5; j < values.size(); j++)
			elementNodes.push_back(nodes_[values[j] - 1]);
		std::string name = physicalEntities[values[3]];
		Surface* object = geometry_->getSurface(name);
		BaseSurfaceElement* base_elem = new BaseSurfaceElement(i, *type, elementNodes);
		object->addBaseElement(base_elem);
		object->addNodes(elementNodes);
		switch (elementType)
		{
		case T3:
		case T6:
		case T10:
		case Q4:
		case Q9:
		case Q16:
			base_elem->setPlot(true);
			Material* mat = object->getMaterial();
			int ndofs = ndofs_per_node * elementNodes.size();
			std::vector<DegreeOfFreedom*> dofs;
			dofs.reserve(ndofs);
			for (auto& node : elementNodes)
				for (int j = 0; j < dimension_; j++)
					dofs.push_back(node->getDegreeOfFreedom(j));
			if (mixed)
				for (auto& node : elementNodes)
					dofs.push_back(node->getDegreeOfFreedom(dimension_));
			Element* element = new PlaneElement(i, dofs, mat, base_elem, parameters_);
			elements_.push_back(element);
			object->addElement(element);
			break;
		}
	}
	for (int i = 0; i < nVolumeElements; i++)
	{
		std::getline(file, line);
		std::vector<std::string> tokens = split(line);
		std::vector<int> values(tokens.size(), 0);
		for (size_t j = 0; j < tokens.size(); j++)
			std::istringstream(tokens[j]) >> values[j];
		ParametricVolumeElement* type = static_cast<ParametricVolumeElement*>(gmshElements.at(values[1]));
		int numberOfNodes = type->getNumberOfNodes();
		std::vector<Node*> elementNodes;
		elementNodes.reserve(numberOfNodes);
		for (size_t j = 5; j < values.size(); j++)
			elementNodes.push_back(nodes_[values[j] - 1]);
		std::string name = physicalEntities[values[3]];
		Volume* object = geometry_->getVolume(name);
		BaseVolumeElement* base_elem = new BaseVolumeElement(i, *type, elementNodes);
		object->addBaseElement(base_elem);
		object->addNodes(elementNodes);
		switch (elementType)
		{
		case TET4:
		case TET10:
		case TET20:
		case HEX8:
		case HEX27:
		case HEX64:
			base_elem->setPlot(true);
			Material* mat = object->getMaterial();
			int ndofs = ndofs_per_node * elementNodes.size();
			std::vector<DegreeOfFreedom*> dofs;
			dofs.reserve(ndofs);
			for (auto& node : elementNodes)
				for (int j = 0; j < dimension_; j++)
					dofs.push_back(node->getDegreeOfFreedom(j));
			if (mixed)
				for (auto& node : elementNodes)
					dofs.push_back(node->getDegreeOfFreedom(dimension_));
			Element* element = new VolumeElement(i, dofs, mat, base_elem, parameters_);
			elements_.push_back(element);
			object->addElement(element);
			break;
		}
	}
	for (auto& pair : geometry_->lines_)
	{
		auto& line = pair.second;
		line->elements_.shrink_to_fit();
		line->nodes_.shrink_to_fit();
	}
	for (auto& pair : geometry_->surfaces_)
	{
		auto& surface = pair.second;
		surface->elements_.shrink_to_fit();
		surface->nodes_.shrink_to_fit();
	}
	for (auto& pair : geometry_->volumes_)
	{
		auto& volume = pair.second;
		volume->elements_.shrink_to_fit();
		volume->nodes_.shrink_to_fit();
	}

	//Closing the file
	file.close();
	if (deleteFiles)
		system(("rm " + inputFile).c_str());

	//Set node materials
	for (auto& pair : geometry_->volumes_)
	{
		auto& volume = pair.second;
		Material* mat = volume->getMaterial();
		if (mat)
			for (auto& node : volume->nodes_)
				node->setMaterial(mat);
	}
	for (auto& pair : geometry_->surfaces_)
	{
		auto& surface = pair.second;
		Material* mat = surface->getMaterial();
		if (mat)
			for (auto &node : surface->nodes_)
				node->setMaterial(mat);
	}
	for (auto &pair : geometry_->lines_)
	{
		auto &line = pair.second;
		Material *mat = line->getMaterial();
		if (mat)
			for (auto &node : line->nodes_)
				node->setMaterial(mat);
	}

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Data reading executed. Elapsed time: %f\n", elapsed.count());
	PetscPrintf(PETSC_COMM_WORLD, "Model information: \n%d Solid Nodes \n%d Solid Elements \n%d Solid DOFs \n", 
								   nodes_.size(), elements_.size(), numberOfDOFs_);
}

void SolidDomain::transferGeometricBoundaryConditions()
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	// Transfering neumann conditions
	int index = 0;
	for (auto& nbc : geometry_->getNeumannBoundaryConditions())
	{
		// Transfering point loads applied to points
		if (nbc->getPoint())
		{
			Point* p = nbc->getPoint();
			Node* n = p->getNode();
			double x = nbc->getValueX();
			double y = nbc->getValueY();
			double z = nbc->getValueZ();
			neumannBoundaryConditions_.emplace_back(new PointLoad(index, dimension_, n, x, y, z));
			index++;
		}
		// Transfering distributed loads applied to line
		else if (nbc->getLine())
		{
			Line* l = nbc->getLine();
			const std::vector<BaseLineElement*>& elements = l->getBaseElements();
			double x = nbc->getValueX();
			double y = nbc->getValueY();
			double z = nbc->getValueZ();
			for (auto& el : elements)
			{
				int ndofs = dimension_ * el->getNumberOfNodes();
				neumannBoundaryConditions_.emplace_back(new LineLoad(index, ndofs, el, x, y, z));
				index++;
			}
		}
		// Transfering distributed loads applied to surface
		else if (nbc->getSurface())
		{
			Surface* s = nbc->getSurface();
			const std::vector<BaseSurfaceElement*>& elements = s->getBaseElements();
			double x = nbc->getValueX();
			double y = nbc->getValueY();
			double z = nbc->getValueZ();
			for (auto& el : elements)
			{
				int ndofs = dimension_ * el->getNumberOfNodes();
				neumannBoundaryConditions_.emplace_back(new SurfaceLoad(index, ndofs, el, x, y, z));
				index++;
			}
		}
	}

	// Transfering dirichlet conditions
	for (auto& dbc : geometry_->getDirichletBoundaryConditions())
	{
		// Transfering nodal displacements applied to points
		if (dbc->getPoint())
		{
			Point* p = dbc->getPoint();
			Node* n = p->getNode();
			n->setConstrain(true);
			n->setBoundary(true);
			double val = dbc->getValue();
			ConstrainedDOF constrainedDOF = dbc->getDegreeOfFreedom();
			switch (constrainedDOF)
			{
				case X:
				{
					DegreeOfFreedom* dof = n->getDegreeOfFreedom(0);
					dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
					break;
				}
				case Y:
				{
					DegreeOfFreedom* dof = n->getDegreeOfFreedom(1);
					dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
					break;
				}
				case Z:
				{
					DegreeOfFreedom* dof = n->getDegreeOfFreedom(2);
					dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
					break;
				}
				case ALL:
				{
					for (int i = 0; i < dimension_; i++)
					{
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(i);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
					}
					break;
				}
			}
		}
		// Transfering nodal displacements applied to lines
		else if (dbc->getLine())
		{
			Line* l = dbc->getLine();
			double val = dbc->getValue();
			ConstrainedDOF constrainedDOF = dbc->getDegreeOfFreedom();
			const std::vector<Node*>& nodes = l->getNodes();
			switch (constrainedDOF)
			{
				case X:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(0);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case Y:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(1);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case Z:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(2);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case ALL:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						for (int i = 0; i < dimension_; i++)
						{
							DegreeOfFreedom* dof = n->getDegreeOfFreedom(i);
							dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
						}
					}
					break;
				}
			}
		}
		// Transfering nodal displacements applied to surfaces
		else if(dbc->getSurface())
		{
			Surface* s = dbc->getSurface();
			double val = dbc->getValue();
			ConstrainedDOF constrainedDOF = dbc->getDegreeOfFreedom();
			const std::vector<Node*>& nodes = s->getNodes();
			switch (constrainedDOF)
			{
				case X:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(0);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case Y:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(1);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case Z:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(2);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case ALL:
				{
					for (auto& n : nodes)
					{
						n->setConstrain(true);
						n->setBoundary(true);
						for (int i = 0; i < dimension_; i++)
						{
							DegreeOfFreedom* dof = n->getDegreeOfFreedom(i);
							dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
						}
					}
					break;
				}
			}
		}
		// Transfering nodal displacements applied to volumes
		else if(dbc->getVolume())
		{
			Volume* v = dbc->getVolume();
			double val = dbc->getValue();
			ConstrainedDOF constrainedDOF = dbc->getDegreeOfFreedom();
			const std::vector<Node*>& nodes = v->getNodes();
			switch (constrainedDOF)
			{
				case X:
				{
					for (auto& n : nodes)
					{
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(0);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case Y:
				{
					for (auto& n : nodes)
					{
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(1);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case Z:
				{
					for (auto& n : nodes)
					{
						DegreeOfFreedom* dof = n->getDegreeOfFreedom(2);
						dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));	
					}
					break;
				}
				case ALL:
				{
					for (auto& n : nodes)
					{
						for (int i = 0; i < dimension_; i++)
						{
							DegreeOfFreedom* dof = n->getDegreeOfFreedom(i);
							dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(n, dof, val));
						}
					}
					break;
				}
			}
		}
	}

	// Transfering interface conditions
	int numberOfInterfaceNodes = 0;
	for (auto& ibc : geometry_->getInterfaceBoundaryConditions())
	{
		if (ibc->getPoint()) numberOfInterfaceNodes++;
		else if (ibc->getLine()) numberOfInterfaceNodes += ibc->getLine()->getNumberOfNodes();
		else if (ibc->getSurface()) numberOfInterfaceNodes += ibc->getSurface()->getNumberOfNodes();
		else if (ibc->getVolume()) numberOfInterfaceNodes += ibc->getVolume()->getNumberOfNodes();
	}
	interfaceNodes_.reserve(numberOfInterfaceNodes);

	for (auto& ibc : geometry_->getInterfaceBoundaryConditions())
	{
		// Transfering interface points
		if (ibc->getPoint())
		{
			Point* p = ibc->getPoint();
			Node* n = p->getNode();
			if (!n->isInterface())
			{
				n->setInterface(true);
				interfaceNodes_.push_back(n);
			}
		}
		// Transfering interface lines
		else if (ibc->getLine())
		{
			Line* l = ibc->getLine();
			const std::vector<Node*>& nodes = l->getNodes();
			for (auto& n : nodes)
			{
				if (!n->isInterface())
				{
					n->setInterface(true);
					interfaceNodes_.push_back(n);
				}
			}
		}
		// Transfering interface surfaces
		else if(ibc->getSurface())
		{
			Surface* s = ibc->getSurface();
			const std::vector<Node*>& nodes = s->getNodes();
			for (auto& n : nodes)
			{
				if (!n->isInterface())
				{
					n->setInterface(true);
					interfaceNodes_.push_back(n);
				}
			}
		}
	}

	interfaceNodes_.shrink_to_fit();

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Transfering Boundary Conditions from Entities to Model. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::domainDecompositionMETIS(const PartitionOfUnity& elementType)
{
	auto start_timer = std::chrono::high_resolution_clock::now();

    const std::unordered_map<PartitionOfUnity, int> numNodes = 
	{
		//lines
		{L2, 2}, {L3, 3}, {L4, 4},
		//triangles
		{T3, 3}, {T6, 6}, {T10, 10},
		//quadrilaterals
		{Q4, 4}, {Q9, 9}, {Q16, 16},
		//tetrahedrons
		{TET4, 4}, {TET10, 10}, {TET20, 20}
	};
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = elements_.size();
    idx_t numNd = nodes_.size();
    idx_t ssize = size;
    idx_t one = 1;
	idx_t n = numNodes.at(elementType);
	idx_t* elem_start = new idx_t[numEl+1]; idx_t* elem_connec = new idx_t[n*numEl];
    elementPartition_ = new idx_t[numEl];
    nodePartition_ = new idx_t[numNd];
    for (idx_t i = 0; i < numEl+1; i++)
	{
        elem_start[i]=n*i;
    }
    for (idx_t jel = 0; jel < numEl; jel++)
	{
        for (idx_t i=0; i < elements_[jel]->getNodes().size(); i++)
		{	
			int nodeIndex = elements_[jel]->getNodes()[i]->getIndex();
        	elem_connec[n*jel+i] = nodeIndex;
        }
    }

    //Performs the domain decomposition
    METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec, \
                              NULL, NULL, &one, &ssize, NULL, NULL,    \
                              &objval, elementPartition_, nodePartition_);
	delete[] elem_start; delete[] elem_connec;
   
	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "6/6 - Partitioning the domain. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::nodalNeighborSearch() const
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	//clearing the neighborhood information of the previous mesh
	for (Node* const& node : nodes_)
	{
		node->clearNeighborNodes();
		node->clearNeighborElements();
	}

	//searching the neighbor elements of each node
	for (Element* const& el : elements_)
	{
		const std::vector<Node*>& nodes = el->getNodes();
		for (Node* const& node : nodes)
		{
			node->addNeighborElement(el);
		}
	}

	//searching the neighbor nodes of each node
	for (Node* const& node : nodes_)
	{
		node->addNeighborNode(node);
		const std::vector<Element*>& neighborElements = node->getNeighborElements();
		for (auto& el : neighborElements)
		{
			const std::vector<Node*>& elemNodes = el->getNodes();
			for (auto& elemNode : elemNodes)
			{
				bool alreadyAdded = false;
				for (auto& n : node->getNeighborNodes())
				{
					if (elemNode->getIndex() == n->getIndex())
					{
						alreadyAdded = true;
						break;
					}
				}
				if (!alreadyAdded)
				{
					node->addNeighborNode(elemNode);
				}
			}
		}
	}

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;
}

void SolidDomain::elementalNeighborSearch() const
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	//Reseting the boundary flags for nodes and elements, except for the nodes that are constrained
    for (Node* const& node : nodes_)
        if (!node->isConstrained())
            node->setBoundary(false);
    
    for (Element* const& el : elements_)
        el->setBoundary(false);

	//clearing the neighborhood information of the previous mesh
	for (Node* const& node : nodes_)
	{
		node->clearNeighborElements();
	}
	for (Element* const& element : elements_)
	{
		element->clearNeighborElements();
	}

	//searching the neighbor elements of each node
	for (Element* const& el : elements_)
	{
		const std::vector<Node*>& nodes = el->getNodes();
		for (Node* const& node : nodes)
		{
			node->addNeighborElement(el);
		}
	}

	//searching the neighbor elements of each element
	for (Element* const& el : elements_)
	{
		ParametricElement* parametric = el->getBaseElement()->getParametricElement();
		int numberOfFaces = parametric->getNumberOfFaces();
		const std::vector<Node*>& elementNodes = el->getNodes();

		for (int i = 0; i < numberOfFaces; i++)
		{
			bool isBoundary = true;
			const std::vector<int>& faceVerticesIndex = parametric->getFaceVertices(i);
			const std::vector<int>& faceNodesIndex = parametric->getFaceNodes(i);
			unsigned int numberOfVerticesPerFace = faceVerticesIndex.size();
			unsigned int numberOfNodesPerFace = faceNodesIndex.size();
			std::vector<Node*> faceVertices(numberOfVerticesPerFace);
			for (int j = 0; j < numberOfVerticesPerFace; j++)
			{
				faceVertices[j] = elementNodes[faceVerticesIndex[j]];
			}
			std::vector<Node*> faceNodes(numberOfNodesPerFace);
			for (int j = 0; j < numberOfNodesPerFace; j++)
			{
				faceNodes[j] = elementNodes[faceNodesIndex[j]];
			}

			const std::vector<Element*>& neighborElements = faceVertices[0]->getNeighborElements();

			for (Element* const& nElem : neighborElements)
			{
				if (nElem->getIndex() == el->getIndex()) continue;
				int numberOfCoincidentVertices = 1; //it is 1 because face vertex 0 is always coincident as we are searching on its neighbors
				const std::vector<Node*>& nNodes = nElem->getNodes();
				for (Node* const& node : nNodes)
				{
					for (int j = 1; j < numberOfVerticesPerFace; j++)
					{
						if (node->getIndex() == faceVertices[j]->getIndex())
						{
							numberOfCoincidentVertices++;
							break;
						}
					}
				}
				if (numberOfCoincidentVertices == numberOfVerticesPerFace)
				{
					el->addNeighborElement(nElem);
					isBoundary = false;
					break;
				}	
			}
			if (isBoundary) //no neighbor element was found for face i, so we add the element itself to the neighbor elements container
			{
				el->addNeighborElement(el);
				el->setBoundary(true);
				for (Node* const& node : faceNodes)
					node->setBoundary(true);
			}
		}
	}

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;
}

void SolidDomain::createSystemMatrix(Mat& mat)
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	//getting processor rank and total number of processors
	int rank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

	int N = numberOfBlockedDOFs_;
	int nb = numberOfBlockedNodes_;

	//defining the ownership range of each processor
    int start[size], end[size], n = PETSC_DECIDE;
    PetscSplitOwnership(PETSC_COMM_WORLD, &n, &N);
    MPI_Scan(&n, &end[rank], 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    start[rank] = end[rank] - n;

    for (int i = 0; i < size; i++)
    {
        MPI_Bcast(&start[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
        MPI_Bcast(&end[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
    }

	//Number of nonzero entries in diagonal and off diagonal parts of the matrix
	int* d_nnz_all = new int[N];
    int* o_nnz_all = new int[N];

	if (rank == 0)
	{
		//Counting the number of nonzeroes in the global matrix
		int nnz = 0;
		for (Node*& node : nodes_)
		{
			if (!node->isIsolated())
			{
				int ndof = node->getNumberOfDegreesOfFreedom();
				const std::vector<Node*>& neighborNodes = node->getNeighborNodes();
				for (auto& neighborNode : neighborNodes)
				{
					nnz += ndof * neighborNode->getNumberOfDegreesOfFreedom();
				}
			}
		}

		int* col = new int[nnz]; //stores the number of the column of each nonzero term
		int* ptr = new int[N+1]; //accumulates the number of nonzeroes per row in a CSR format
		ptr[0] = 0;

		//Iterating over the nodes in the permuted order. This block of code needs to acces the dofs in a crescent and consecutive order
		//As only the dofs of blocked nodes are contributing to the matrix, we access only the firsts nblocked nodes
		int sum = -1;
		for (unsigned int k = 0; k < nb; k++)
		{
			Node*& node = nodes_[perm_[k]]; 
			const std::vector<DegreeOfFreedom*>& i_dofs = node->getDegreesOfFreedom();
			for (DegreeOfFreedom* const& i_dof : i_dofs)
			{
				int i = i_dof->getIndex();
				const std::vector<Node*>& neighborNodes = node->getNeighborNodes();
				for (auto& neighborNode : neighborNodes)
				{
					const std::vector<DegreeOfFreedom*>& j_dofs = neighborNode->getDegreesOfFreedom();
					for (DegreeOfFreedom* const& j_dof : j_dofs)
					{
						int j = j_dof->getIndex();
						col[++sum] = j;
					}
				}
				ptr[i+1] = sum+1;
			}
		}

		// number of non-zeros by row
        for (int k = 0; k < size; k++)
        {
            int low = start[k];
            int high = end[k];
            for (int i = low; i < high; i++)
            {
                int inf = ptr[i];
                int sup = ptr[i + 1];
                d_nnz_all[i] = 0;
                o_nnz_all[i] = 0;
                for (int j = inf; j < sup; j++)
                {
                    if (col[j] >= low && col[j] < high)
						++d_nnz_all[i];
                    else
						++o_nnz_all[i];
                }
            }
        }
		delete[] ptr;
        delete[] col;
	}

	MPI_Bcast(d_nnz_all, N, MPI_INT, 0, PETSC_COMM_WORLD);
    MPI_Bcast(o_nnz_all, N, MPI_INT, 0, PETSC_COMM_WORLD);

    int* d_nnz = &d_nnz_all[start[rank]];
    int* o_nnz = &o_nnz_all[start[rank]];

    // Create PETSc sparse matrix
	MatCreateAIJ(PETSC_COMM_WORLD, n, n, N, N, 0, d_nnz, 0, o_nnz, &mat);

	delete[] d_nnz_all;
    delete[] o_nnz_all;

	MatSetFromOptions(mat);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	//PetscPrintf(PETSC_COMM_WORLD, "Allocating memory for the sparse matrix. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::reorderDOFs()
{
	//getting processor rank and total number of processors
	int rank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

	auto start_timer = std::chrono::high_resolution_clock::now();

	int n = numberOfNodes_;
	int nb = numberOfBlockedNodes_;
	int m = 10 * n;
	
    if (!perm_) perm_ = new idx_t[n]; //array of index permutation according to METIS. perm[i] is the ith node in container

	if (rank == 0)
	{
		idx_t* xadj = new idx_t[n+1];
		xadj[0] = 0;
		std::vector<int> adjacency_aux;
		adjacency_aux.reserve(m);

		for (Node*& node : nodes_)
		{
			int index = node->getIndex();
			const std::vector<Node*>& neighborNodes = node->getNeighborNodes();
			unsigned int num_neighbor_nodes = neighborNodes.size() - 1;
			xadj[index+1] = xadj[index] + num_neighbor_nodes;
			for (unsigned int i = 0; i < num_neighbor_nodes; i++)
			{
				adjacency_aux.push_back(neighborNodes[i+1]->getIndex());
			}
		}

		m = adjacency_aux.size(); //2x number of edges
		idx_t* adjncy = new idx_t[m];

		for (int i = 0; i < m; i++)
		{
			adjncy[i] = adjacency_aux[i];
		}

		adjacency_aux.clear();
		adjacency_aux.shrink_to_fit();

		// iperm[i] is the permuted index of initial node i.
		// perm[i] is the initial index of permuted node i.
		// IMPORTANT! The nodes container will not be sorted in relation to the permuted indexes, it will remain fixed.
		// To loop over the nodes in the permuted order, the elements should be accessed via nodes[perm[i]].
		
		idx_t* perm = new idx_t[n];
		idx_t* iperm = new idx_t[n];
		
    	METIS_NodeND(&n, xadj, adjncy, NULL, NULL, perm, iperm);

		//Placing the blocked nodes first and the isolated ones after them
		int countBlocked = -1;
		int countIsolated = -1;
		for (int i = 0; i < n; i++)
		{
			int initial_index = perm[i];
			if (!nodes_[initial_index]->isIsolated())
			{
				perm_[++countBlocked] = initial_index;
			}
			else
			{
				countIsolated++;
				perm_[nb + countIsolated] = initial_index;
			}
		}

		delete[] xadj;
		delete[] adjncy;
		delete[] perm;
		delete[] iperm;
	}

	//Share METIS node renumbering to all processors
    MPI_Bcast(&perm_[0],n,MPI_INT,0,PETSC_COMM_WORLD);

	//Setting the new index for nodes and their DOFs
	int dof_index = -1;
	for (int i = 0; i < n; i++) //i is the new index of node perm[i]
	{
		int initial_index = perm_[i];
		Node* node = nodes_[initial_index];
		node->setPermutedIndex(i);
		std::vector<DegreeOfFreedom*> dofs = node->getDegreesOfFreedom();
		for (DegreeOfFreedom*& dof : dofs)
		{
			dof->setIndex(++dof_index);
		}
	}

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	//PetscPrintf(PETSC_COMM_WORLD, "4/6 - Reordering DOFs. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::domainDecomposition()
{
	/*	The domain decomposition is performed in the following way:
		- The blocked nodes are distributed along the processors according to their dofs.
		- Then, the elements are decomposed based on their nodes ranks.
		- After that, a rank is assigned to the isolated nodes.
	*/
	auto start_timer = std::chrono::high_resolution_clock::now();

	int rank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

	unsigned int n_nodes = nodes_.size();
    unsigned int n_elem = elements_.size();

	elementPartition_ = new idx_t[n_elem];
    nodePartition_ = new idx_t[n_nodes];
    
	//defining the ownership range of blocked DOFs
	int N = numberOfBlockedDOFs_;
    int start[size], end[size], n = PETSC_DECIDE;
    PetscSplitOwnership(PETSC_COMM_WORLD, &n, &N);
    MPI_Scan(&n, &end[rank], 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    start[rank] = end[rank] - n;

    for (int i = 0; i < size; i++)
    {
        MPI_Bcast(&start[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
        MPI_Bcast(&end[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
    }

	if (rank == 0)
	{
        //Performs the domain decomposition
        for (unsigned int i = 0; i < numberOfBlockedNodes_; i++) //Looping over the permuted order of nodes
		{
			Node* node = nodes_[perm_[i]];
			int ndof = node->getNumberOfDegreesOfFreedom();
            for (int j = 0; j < size; j++)
			{
                if (((ndof*i) >= start[j]) && ((ndof*i) <= end[j]))
				{
					node->setRank(j);
					nodePartition_[perm_[i]] = j;
				}
            }
        }

        for (unsigned int i = 0; i < n_elem; i++)
		{
            Element* elem = elements_[i];
			Node* last_node = elem->getNodes().back(); //could be any nodes of the element
			unsigned int last_node_index = last_node->getPermutedIndex();
			int ndof = last_node->getNumberOfDegreesOfFreedom();

            for (int j = 0; j < size; j++)
			{
                if (((ndof*last_node_index) >= start[j]) && ((ndof*last_node_index) <= end[j]))
				{
                    elem->setRank(j);
					elementPartition_[i] = j;
                    break;
                }
            }
        }
    }

	//defining the ownership range of isolated nodes
	N = numberOfNodes_ - numberOfBlockedNodes_;
	if  (N > 0)
	{
		n = PETSC_DECIDE;
		PetscSplitOwnership(PETSC_COMM_WORLD, &n, &N);
    	MPI_Scan(&n, &end[rank], 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD);
    	start[rank] = end[rank] - n;

		for (int i = 0; i < size; i++)
		{
			MPI_Bcast(&start[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
			MPI_Bcast(&end[i], 1, MPI_INT, i, PETSC_COMM_WORLD);
		}
	}

	if (rank == 0)
	{
		for (unsigned int i = numberOfBlockedNodes_; i < n_nodes; i++)
		{
			Node* node = nodes_[perm_[i]];
			int index = i - numberOfBlockedNodes_;
			for (int j = 0; j < size; j++)
			{
				if (index >= start[j] && index <= end[j])
				{
					node->setRank(j);
					nodePartition_[perm_[i]] = j;
				}
			}
		}
	}

	MPI_Bcast(elementPartition_,n_elem,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(nodePartition_,n_nodes,MPI_INT,0,PETSC_COMM_WORLD);

	if (rank != 0)
	{
		for (unsigned int i = 0; i < n_nodes; i++)
		{
			nodes_[i]->setRank(nodePartition_[i]);
		}
		for (unsigned int i = 0; i < n_elem; i++)
		{
			elements_[i]->setRank(elementPartition_[i]);
		}
	}

	MPI_Barrier(PETSC_COMM_WORLD);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	//PetscPrintf(PETSC_COMM_WORLD, "6/6 - Partitioning the domain. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::solveStaticProblem()
{
	auto start_timer = std::chrono::high_resolution_clock::now();
	
	setReferenceConfiguration(ReferenceConfiguration::INITIAL);
	parameters_->setStaticAnalysis(true);

	// Petsc variables
	Mat               tangent;
    Vec               rhs, solution;
    KSP               ksp;
    PC                pc;
	PetscLogDouble    bytes = 0.0;
	PetscBool		  isbjacobi;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	PetscMemoryGetCurrentUsage(&bytes);
	PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes/(1024*1024));

	//Initial position norm
    double initialPositionNorm = getInitialPositionNorm();

	//Array of constrained degrees of freedom
	int numberOfConstrainedDOFs; int* constrainedDOFs;
	getConstrainedDOFs(numberOfConstrainedDOFs, constrainedDOFs);
	
	//Array of external forces
	int ndofsForces; std::vector<DegreeOfFreedom*> dofsForces; double* externalForces;
	getExternalForces(ndofsForces, dofsForces, externalForces);

	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp, KSPFGMRES);
	KSPSetTolerances(ksp, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 100);
	KSPGMRESSetRestart(ksp, 100);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCBJACOBI);

	createSystemMatrix(tangent);

	//Create PETSc vectors
	VecCreate(PETSC_COMM_WORLD, &rhs);
	VecSetSizes(rhs, PETSC_DECIDE, numberOfBlockedDOFs_);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs, &solution);

	const int numberOfSteps = parameters_->getNumberOfSteps();
	const int maxNonlinearIterations = parameters_->getMaxNonlinearIterations();
	const double nonlinearTolerance = parameters_->getNonlinearTolerance();
	
	if (rank == 0) exportToParaview(0);

	for (int step = 0; step < numberOfSteps; step++)
	{
		double loadFactor = ((double)step + 1.0) / (double)numberOfSteps;

		PetscPrintf(PETSC_COMM_WORLD, "\n----------------------- STEP = %d, Loadfactor = %f  -----------------------\n\n", step+1, loadFactor);
		parameters_->setCurrentTime(loadFactor);
		double modelVolume = 0.0;
    	for (Element* const& el : elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		parameters_->setModelVolume(modelVolume);

		double positionNorm, pressureNorm;
		for (int iteration = 0; iteration < maxNonlinearIterations; iteration++)
		{
			setPastVariables();
			computeIntermediateVariables();
			applyNeummanConditions(rhs, ndofsForces, dofsForces, externalForces, loadFactor);
			assembleStaticLinearSystem(tangent, rhs);
			MatZeroRowsColumns(tangent, numberOfConstrainedDOFs, constrainedDOFs, 1.0, solution, rhs);
			MatView(tangent, PETSC_VIEWER_DRAW_WORLD);
			solveLinearSystem(ksp, tangent, rhs, solution);
			updateVariables(solution, positionNorm, pressureNorm);
			setPastVariables();
			computeIntermediateVariables();
			
			PetscMemoryGetCurrentUsage(&bytes);
			PetscPrintf(PETSC_COMM_WORLD, "Newton iteration: %d - L2 Position Norm: %E - L2 Pressure Norm: %E\nMemory used by each processor: %f Mb\n", 
						iteration, positionNorm / initialPositionNorm, pressureNorm, bytes/(1024*1024));
	
			MatZeroEntries(tangent);
			VecZeroEntries(rhs);

			if (positionNorm / initialPositionNorm <= nonlinearTolerance)
				break;
		}

		//export results to paraview
		if (rank == 0)
		{
			computeCauchyStress();
			exportToParaview(step+1);
		}
	}
	delete[] constrainedDOFs;
	delete[] externalForces;
	MatDestroy(&tangent);
	KSPDestroy(&ksp);
    VecDestroy(&rhs);
    VecDestroy(&solution);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Solid Analysis Done. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::assembleStaticLinearSystem(Mat& mat, Vec& vec)
{
	auto start_timer = std::chrono::high_resolution_clock::now();
	
	int rank, size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	for (Element* const& el : elements_)
	{
		if (el->getRank() == rank && el->isActive())
		{
			int ndofsPosition, ndofsPressure;
			int* indexes;
			double *rhsValues, *hessianValues;
			el->elementContributions(ndofsPosition, ndofsPressure, indexes, rhsValues, hessianValues);
			int ndofs = ndofsPosition + ndofsPressure;
			// dispersing element rhs contribution into global rhs vector
			VecSetValues(vec, ndofs, indexes, rhsValues, ADD_VALUES);

			// dispersing element matrices contribution into global tangent matrix
			for (int i = 0; i < ndofsPosition; i++)
			{
				for (int j = 0; j < ndofsPosition; j++)
				{
					double value = hessianValues[ndofs * i + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[i], 1, &indexes[j], &value, ADD_VALUES);
				}
				for (int j = 0; j < ndofsPressure; j++)
				{
					double value = hessianValues[ndofsPosition + ndofs*i + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[i], 1, &indexes[ndofsPosition + j], &value, ADD_VALUES);
				}
			}
			
			for (int i = 0; i < ndofsPressure; i++)
			{
				for (int j = 0; j < ndofsPosition; j++)
				{
					double value = hessianValues[ndofs*(ndofsPosition+i) + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[ndofsPosition + i], 1, &indexes[j], &value, ADD_VALUES);
				}
				for (int j = 0; j < ndofsPressure; j++)
				{
					double value = hessianValues[ndofs*(ndofsPosition+i)+ndofsPosition + j]; //value = 0.0;
					MatSetValues(mat, 1, &indexes[ndofsPosition + i], 1, &indexes[ndofsPosition + j], &value, ADD_VALUES);
				}
			}
			delete[] indexes;
			delete[] rhsValues;
			delete[] hessianValues;
		}
	}

	//Assemble matrices and vectors
    MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);		
    VecAssemblyBegin(vec); 
    VecAssemblyEnd(vec);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Assemble Linear System. Elapsed time: %f\n", elapsed.count());
}

void SolidDomain::solveStaggeredProblem(int& ndofsInterfaceForces,
                                        std::vector<DegreeOfFreedom*>& dofsInterfaceForces,
                                        double*& interfaceForces)
{
	// Petsc variables
	Mat               tangent;
    Vec               rhs, solution;
    KSP               ksp;
    PC                pc;
	PetscLogDouble    bytes = 0.0;
	
	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	//Initial position norm
    double initialPositionNorm = getInitialPositionNorm();

	//Array of constrained degrees of freedom
	int numberOfConstrainedDOFs;
	int* constrainedDOFs;
	getConstrainedDOFs(numberOfConstrainedDOFs, constrainedDOFs);

	//Array of external forces
	int ndofsForces; std::vector<DegreeOfFreedom*> dofsForces; double* externalForces;
	getExternalForces(ndofsForces, dofsForces, externalForces);

	computeCurrentVariables();
	computeIntermediateVariables();

	createSystemMatrix(tangent);

	//Create PETSc vectors
	VecCreate(PETSC_COMM_WORLD, &rhs);
	VecSetSizes(rhs, PETSC_DECIDE, numberOfDOFs_);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs, &solution);
	
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetType(ksp, KSPPREONLY);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCLU);
	PCFactorSetMatSolverType(pc, MATSOLVERMUMPS);
	KSPSetFromOptions(ksp);

	const int maxNonlinearIterations = parameters_->getMaxNonlinearIterations();
	const double nonlinearTolerance = parameters_->getNonlinearTolerance();
	
	double positionNorm, pressureNorm;
	for (int iteration = 0; (iteration < maxNonlinearIterations); iteration++)
	{
		applyNeummanConditions(rhs, ndofsForces, dofsForces, externalForces, 1.0);
		applyNeummanConditions(rhs, ndofsInterfaceForces, dofsInterfaceForces, interfaceForces, 1.0);
		assembleTransientLinearSystem(tangent, rhs);
		MatZeroRowsColumns(tangent, numberOfConstrainedDOFs, constrainedDOFs, 1.0, solution, rhs);
		//MatView(tangent, PETSC_VIEWER_DRAW_WORLD);
		solveLinearSystem(ksp, tangent, rhs, solution);
		updateVariables(solution, positionNorm, pressureNorm);
		computeCurrentVariables();
		computeIntermediateVariables();

		PetscMemoryGetCurrentUsage(&bytes);
		PetscPrintf(PETSC_COMM_WORLD, "Solid iteration: %d - L2 Position Norm: %E - L2 Pressure Norm: %E\n", 
						iteration, positionNorm / initialPositionNorm, pressureNorm);
	
		MatZeroEntries(tangent);
		VecZeroEntries(rhs);

		if (positionNorm / initialPositionNorm <= nonlinearTolerance)
			break;
	}

	delete[] constrainedDOFs;
	delete[] externalForces;
	KSPDestroy(&ksp);
    VecDestroy(&rhs);
    VecDestroy(&solution);
	MatDestroy(&tangent);
}