#include "CoupledDomain.h"

CoupledDomain::CoupledDomain(FluidDomain* fluid, SolidDomain* solid)
	: fluid_(fluid),
      solid_(solid), 
	  dimension_(2),
	  numberOfDOFs_(0),
	  numberOfBlockedDOFs_(0),
	  numberOfIsolatedDOFs_(0),
	  numberOfCoupledDOFs_(0),
	  numberOfNodes_(0),
	  numberOfBlockedNodes_(0),
	  numberOfIsolatedNodes_(0),
      numberOfCoupledNodes_(0),
      numberOfCoupledBlockedDOFs_(0),
	  pressureCoupled_(false)
{
	int fail = system("mkdir -p ./results");
	fail = system("rm ./results/*.vtu 2> /dev/null");
	fail = system("rm ./results/*.txt 2> /dev/null");

    if (solid_->dimension_ == fluid_->dimension_) 
        dimension_ = solid_->dimension_;
    else
    {
        std::cout << "Solid and Fluid domains are not compatible.\n";
        exit(EXIT_FAILURE);
    }
}

CoupledDomain::~CoupledDomain() {}

void CoupledDomain::setNumberOfSteps(const int numberOfSteps)
{
    fluid_->setNumberOfSteps(numberOfSteps);
    solid_->setNumberOfSteps(numberOfSteps);
}

void CoupledDomain::setMaxNonlinearIterations(const int maxNonlinearIterations)
{
    fluid_->setMaxNonlinearIterations(maxNonlinearIterations);
    solid_->setMaxNonlinearIterations(maxNonlinearIterations);
}

void CoupledDomain::setNonlinearTolerance(const double nonlinearTolerance)
{
    fluid_->setNonlinearTolerance(nonlinearTolerance);
    solid_->setNonlinearTolerance(nonlinearTolerance);
}

void CoupledDomain::setDeltat(const double deltat)
{  
    fluid_->setDeltat(deltat);
    solid_->setDeltat(deltat);
}
    
void CoupledDomain::setSpectralRadius(const double rhoInf)
{
    fluid_->setSpectralRadius(rhoInf);
    solid_->setSpectralRadius(rhoInf);
}

void CoupledDomain::setGeneralizedAlphas(const double alphaM, const double alphaF)
{
    fluid_->setGeneralizedAlphas(alphaM, alphaF);
    solid_->setGeneralizedAlphas(alphaM, alphaF);
}

void CoupledDomain::setInitialAccel(const bool initialAccel)
{
    fluid_->setInitialAccel(initialAccel);
    solid_->setInitialAccel(initialAccel);
}

void CoupledDomain::setExportFrequency(const int& freq)
{
	fluid_->setExportFrequency(freq);
    solid_->setExportFrequency(freq);
}

void CoupledDomain::solveMonolithicCoupledProblem()
{
    auto start_timer = std::chrono::high_resolution_clock::now();

    solid_->setReferenceConfiguration(ReferenceConfiguration::INITIAL);

    identifyMonolithicCoupledDOFs();
    coupleInterfaceDOFs();

    numberOfNodes_ = solid_->numberOfNodes_ + fluid_->numberOfNodes_;
    numberOfBlockedNodes_ = solid_->numberOfBlockedNodes_ + fluid_->numberOfBlockedNodes_;
    numberOfCoupledBlockedDOFs_ = numberOfCoupledDOFs_;
	for (Node* const& fluidNode : fluid_->interfaceNodes_)
	{
		if (fluidNode->isIsolated())
		{
			numberOfCoupledBlockedDOFs_ -= dimension_;
		}
	}
    numberOfBlockedDOFs_ = solid_->numberOfBlockedDOFs_ + fluid_->numberOfBlockedDOFs_ - numberOfCoupledBlockedDOFs_;
    numberOfIsolatedNodes_ = fluid_->numberOfIsolatedNodes_;
    numberOfIsolatedDOFs_ = fluid_->numberOfIsolatedDOFs_;

    domainDecomposition();

    if (fluid_->parameters_->getInitialAccel())
    {
        solid_->computeInitialAccel();
        fluid_->computeInitialAccel();
    }

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
	PetscPrintf(PETSC_COMM_SELF,
                "Memory used by each processor to store problem data: %f Mb\n",
                bytes/(1024*1024));

    //Initial position norm
    //double initialPositionNorm = getInitialPositionNorm();

	//Array of constrained degrees of freedom
	int numberOfConstrainedDOFs;
	int* constrainedDOFs;
	getConstrainedDOFs(numberOfConstrainedDOFs, constrainedDOFs);

	//Array of external forces
	int ndofsForces; std::vector<DegreeOfFreedom*> dofsForces; double* externalForces;
	getExternalForces(ndofsForces, dofsForces, externalForces);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetFromOptions(ksp);
	KSPSetType(ksp, KSPFGMRES);
	KSPSetTolerances(ksp, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
	KSPGMRESSetRestart(ksp, 500);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCBJACOBI);

	createSystemMatrix(tangent);

	//Create PETSc vectors
	VecCreate(PETSC_COMM_WORLD, &rhs);
	VecSetSizes(rhs, PETSC_DECIDE, numberOfBlockedDOFs_);
	VecSetFromOptions(rhs);
	VecDuplicate(rhs, &solution);

	const int numberOfSteps = fluid_->parameters_->getNumberOfSteps();
	const int maxNonlinearIterations = fluid_->parameters_->getMaxNonlinearIterations();
	const double nonlinearTolerance = fluid_->parameters_->getNonlinearTolerance();

    if (rank == 0)
    {
        solid_->exportToParaview(0);
        fluid_->exportToParaview(0);
    }
	
	for (int timeStep = 0; timeStep < numberOfSteps; timeStep++)
	{
		PetscPrintf(PETSC_COMM_WORLD,
                    "\n----------------------- TIME STEP = %d, time = %f  -----------------------\n\n",
                    timeStep+1, (double)(timeStep+1)*fluid_->parameters_->getDeltat());
		solid_->parameters_->setCurrentTime(solid_->parameters_->getDeltat() * (double)(timeStep+1));
        fluid_->parameters_->setCurrentTime(fluid_->parameters_->getDeltat() * (double)(timeStep+1));

        double modelVolume = 0.0;
    	for (Element* const& el : solid_->elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		solid_->parameters_->setModelVolume(modelVolume);
        modelVolume = 0.0;
    	for (Element* const& el : fluid_->elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		fluid_->parameters_->setModelVolume(modelVolume);

		solid_->setPastVariables();
		solid_->computeCurrentVariables();
		solid_->computeIntermediateVariables();
        fluid_->setPastVariables();
        fluid_->computeCurrentVariables();
        fluid_->computeIntermediateVariables();

		double positionNorm, pressureNorm, initialPositionNorm, initialPressureNorm;
		for (int iteration = 0; (iteration < maxNonlinearIterations); iteration++)
		{
			if (iteration == 0)
			{
				initialPositionNorm = computePositionNorm();
				initialPressureNorm = computePressureNorm();
			}
			applyNeummanConditions(rhs, ndofsForces, dofsForces, externalForces, 1.0);
			assembleMonolithicLinearSystem(tangent, rhs);
			MatZeroRowsColumns(tangent, numberOfConstrainedDOFs, constrainedDOFs, 1.0, solution, rhs);
            MatView(tangent, PETSC_VIEWER_DRAW_WORLD);
			solveLinearSystem(ksp, tangent, rhs, solution);
			updateVariables(solution, positionNorm, pressureNorm);
			positionNorm /= initialPositionNorm;
			pressureNorm /= initialPressureNorm;
			solid_->computeCurrentVariables();
		    solid_->computeIntermediateVariables();
            fluid_->computeCurrentVariables();
            fluid_->computeIntermediateVariables();

			PetscMemoryGetCurrentUsage(&bytes);
			PetscPrintf(PETSC_COMM_WORLD,
                        "Newton iteration: %d - L2 Position Norm: %E - L2 Pressure Norm: %E\nMemory used by each processor: %f Mb\n", 
						iteration, positionNorm, pressureNorm, bytes/(1024*1024));
	
			MatZeroEntries(tangent);
			VecZeroEntries(rhs);

			if (positionNorm <= nonlinearTolerance && pressureNorm <= nonlinearTolerance)
				break;
		}

		//export results to paraview
		if (rank == 0 && ((timeStep + 1) % solid_->parameters_->getExportFrequency() == 0))
		{
			solid_->computeCauchyStress();
            fluid_->computeCauchyStress();
			solid_->exportToParaview(timeStep+1);
            fluid_->exportToParaview(timeStep+1);
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

	PetscPrintf(PETSC_COMM_WORLD, "Monolithic Coupled Analysis Done. Elapsed time: %f\n", elapsed.count());
}

void CoupledDomain::solveMonolithicCoupledPFEMProblem()
{
    auto start_timer = std::chrono::high_resolution_clock::now();

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    solid_->setReferenceConfiguration(ReferenceConfiguration::INITIAL);

	identifyMonolithicCoupledDOFs();

    numberOfNodes_ = solid_->numberOfNodes_ + fluid_->numberOfNodes_;
    numberOfBlockedNodes_ = solid_->numberOfBlockedNodes_ + fluid_->numberOfBlockedNodes_;
	numberOfCoupledBlockedDOFs_ = numberOfCoupledDOFs_;
	for (Node* const& fluidNode : fluid_->interfaceNodes_)
	{
		if (fluidNode->isIsolated())
		{
			numberOfCoupledBlockedDOFs_ -= dimension_;
		}
	}
    numberOfBlockedDOFs_ = solid_->numberOfBlockedDOFs_ + fluid_->numberOfBlockedDOFs_ - numberOfCoupledBlockedDOFs_;
    numberOfIsolatedNodes_ = fluid_->numberOfIsolatedNodes_;
    numberOfIsolatedDOFs_ = fluid_->numberOfIsolatedDOFs_;

    if (fluid_->parameters_->getInitialAccel())
    {
        solid_->computeInitialAccel();
        fluid_->computeInitialAccel();
    }
	
	// Petsc variables
	Mat               tangent;
    Vec               rhs, solution;
    KSP               ksp;
    PC                pc;
	PetscLogDouble    bytes = 0.0;
	PetscBool		  isbjacobi;

	PetscMemoryGetCurrentUsage(&bytes);
	PetscPrintf(PETSC_COMM_WORLD, "Memory used by each processor to store problem data: %f Mb\n", bytes/(1024*1024));

    //Initial position norm
    //double initialPositionNorm = getInitialPositionNorm();

	//Array of external forces
	int ndofsForces; std::vector<DegreeOfFreedom*> dofsForces; double* externalForces;
	getExternalForces(ndofsForces, dofsForces, externalForces);
	
	KSPCreate(PETSC_COMM_WORLD, &ksp);
	KSPSetFromOptions(ksp);
	KSPSetType(ksp, KSPFGMRES);
	KSPSetTolerances(ksp, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 1000);
	KSPGMRESSetRestart(ksp, 500);
	KSPGetPC(ksp, &pc);
	PCSetType(pc, PCBJACOBI);
	// KSPSetType(ksp, KSPPREONLY);
	// KSPGetPC(ksp, &pc);
	// PCSetType(pc, PCLU);

	const int numberOfSteps = fluid_->parameters_->getNumberOfSteps();
	const int maxNonlinearIterations = fluid_->parameters_->getMaxNonlinearIterations();
	const double nonlinearTolerance = fluid_->parameters_->getNonlinearTolerance();

	if (rank == 0)
    {
        solid_->exportToParaview(0);
        fluid_->exportToParaview(0);
    }
	
	for (int timeStep = 0; timeStep < numberOfSteps; timeStep++)
	{
		PetscPrintf(PETSC_COMM_WORLD,
                    "\n----------------------- TIME STEP = %d, time = %f  -----------------------\n\n",
                    timeStep+1, (double)(timeStep+1)*fluid_->parameters_->getDeltat());
		solid_->parameters_->setCurrentTime(solid_->parameters_->getDeltat() * (double)(timeStep+1));
        fluid_->parameters_->setCurrentTime(fluid_->parameters_->getDeltat() * (double)(timeStep+1));
		
        double modelVolume = 0.0;
    	for (Element* const& el : solid_->elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		solid_->parameters_->setModelVolume(modelVolume);
        modelVolume = 0.0;
    	for (Element* const& el : fluid_->elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		fluid_->parameters_->setModelVolume(modelVolume);

		fluid_->executeMesh();
		fluid_->nodalNeighborSearch();

		modelVolume = 0.0;
    	for (Element* const& el : fluid_->elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		
		if (rank == 0)
		{
			std::ofstream massConservation("results/mass.txt", std::ofstream::app);
			massConservation << (timeStep + 1) << " ";
			massConservation.precision(16);
			massConservation << std::fixed << modelVolume << " ";
			massConservation.close();
		}

		//Identification of isolated particles in the fluid domain
		for (Node* const& node : fluid_->nodes_)
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
		fluid_->numberOfNodes_ = fluid_->nodes_.size();
		fluid_->numberOfDOFs_ = 0;
		for (Node* const& node : fluid_->nodes_)
			fluid_->numberOfDOFs_ += node->getNumberOfDegreesOfFreedom();
		fluid_->numberOfIsolatedNodes_ = 0;
		fluid_->numberOfIsolatedDOFs_ = 0;
		fluid_->numberOfBlockedNodes_ = fluid_->numberOfNodes_;
		fluid_->numberOfBlockedDOFs_ = fluid_->numberOfDOFs_;
		for (Node* const& node : fluid_->nodes_)
		{
			node->setNewEntity(false);
			unsigned int num_neighbor_nodes = node->getNeighborNodes().size();
			if (num_neighbor_nodes == 1)
			{
				node->setIsolated(true);
				fluid_->numberOfBlockedNodes_ -= 1;
				fluid_->numberOfBlockedDOFs_ -= node->getNumberOfDegreesOfFreedom();
				node->getDegreeOfFreedom(dimension_)->setCurrentValue(0.0);
				node->getDegreeOfFreedom(dimension_)->setCurrentFirstTimeDerivative(0.0);
				node->getDegreeOfFreedom(dimension_)->setCurrentSecondTimeDerivative(0.0);

				if (!node->isConstrained() && !node->isInterface())
				{
					fluid_->numberOfIsolatedNodes_++;
					fluid_->numberOfIsolatedDOFs_ += node->getNumberOfDegreesOfFreedom();
				}
			}
		}

        numberOfNodes_ = solid_->numberOfNodes_ + fluid_->numberOfNodes_;
        numberOfBlockedNodes_ = solid_->numberOfBlockedNodes_ + fluid_->numberOfBlockedNodes_;
        numberOfDOFs_ = solid_->numberOfDOFs_ + fluid_->numberOfDOFs_ - numberOfCoupledDOFs_;
        numberOfCoupledBlockedDOFs_ = numberOfCoupledDOFs_;
		for (Node* const& fluidNode : fluid_->interfaceNodes_)
		{
			if (fluidNode->isIsolated())
			{
				numberOfCoupledBlockedDOFs_ -= dimension_;
			}
		}
    	numberOfBlockedDOFs_ = solid_->numberOfBlockedDOFs_ + fluid_->numberOfBlockedDOFs_ - numberOfCoupledBlockedDOFs_;
        numberOfIsolatedNodes_ = fluid_->numberOfIsolatedNodes_;
        numberOfIsolatedDOFs_ = fluid_->numberOfIsolatedDOFs_;

		fluid_->reorderDOFs();

		coupleInterfaceDOFs();
		domainDecomposition();
		fluid_->deactivateSlivers();

		//Array of constrained degrees of freedom
		int numberOfConstrainedDOFs; int* constrainedDOFs;
		getConstrainedDOFs(numberOfConstrainedDOFs, constrainedDOFs);
		
		solid_->setPastVariables();
		solid_->computeCurrentVariables();
		solid_->computeIntermediateVariables();
        fluid_->setPastVariables();
        fluid_->computeCurrentVariables();
        fluid_->computeIntermediateVariables();

		//Create PETSc vectors
		VecCreate(PETSC_COMM_WORLD, &rhs);
		VecSetSizes(rhs, PETSC_DECIDE, numberOfBlockedDOFs_);
		VecSetFromOptions(rhs);
		VecDuplicate(rhs, &solution);

		createSystemMatrix(tangent);

		double positionNorm, pressureNorm, initialPositionNorm, initialPressureNorm;
		for (int iteration = 0; (iteration < maxNonlinearIterations); iteration++)
		{
			if (iteration == 0)
			{
				initialPositionNorm = computePositionNorm();
				initialPressureNorm = computePressureNorm();
			}
			applyNeummanConditions(rhs, ndofsForces, dofsForces, externalForces, 1.0);
			assembleMonolithicLinearSystem(tangent, rhs);
            MatView(tangent, PETSC_VIEWER_DRAW_WORLD);
			MatZeroRowsColumns(tangent, numberOfConstrainedDOFs, constrainedDOFs, 1.0, solution, rhs);
			solveLinearSystem(ksp, tangent, rhs, solution);
			updateVariables(solution, positionNorm, pressureNorm);
			positionNorm /= initialPositionNorm;
			pressureNorm /= initialPressureNorm;
			solid_->computeCurrentVariables();
		    solid_->computeIntermediateVariables();
            fluid_->computeCurrentVariables();
            fluid_->computeIntermediateVariables();

			PetscMemoryGetCurrentUsage(&bytes);
			PetscPrintf(PETSC_COMM_WORLD, "Newton iteration: %d - L2 Position Norm: %E - L2 Pressure Norm: %E\nMemory used by each processor: %f Mb\n", 
						iteration, positionNorm, pressureNorm, bytes/(1024*1024));

			if (positionNorm <= nonlinearTolerance && pressureNorm <= nonlinearTolerance)
				break;
			
			MatZeroEntries(tangent);
			VecZeroEntries(rhs);
		}
		MatDestroy(&tangent);
		VecDestroy(&rhs);
    	VecDestroy(&solution);

        //Compute the isolated particles motion
		for (Node* const& node : fluid_->nodes_)
		{
			if (node->isIsolated() && !node->isConstrained() && !node->isInterface())
			{
				double* gravity = fluid_->parameters_->getGravity();
				for (int j = 0; j < dimension_; j++)
				{
					DegreeOfFreedom* dof = node->getDegreeOfFreedom(j);
					double acel = gravity[j];
					dof->setCurrentSecondTimeDerivative(acel);
					double vel = dof->getPastFirstTimeDerivative() + acel * fluid_->parameters_->getDeltat();
					dof->setCurrentFirstTimeDerivative(vel);
					double coord = dof->getPastValue() + vel * fluid_->parameters_->getDeltat() + 0.5 * fluid_->parameters_->getDeltat() * fluid_->parameters_->getDeltat();
					dof->setCurrentValue(coord);
				}
			}
		}

		modelVolume = 0.0;
    	for (Element* const& el : fluid_->elements_)
		{
        	modelVolume += el->getBaseElement()->getJacobianIntegration();
		}
		if (rank == 0)
		{
			std::ofstream massConservation("results/mass.txt", std::ofstream::app);
			massConservation.precision(16);
			massConservation << std::fixed << modelVolume << "\n";
			massConservation.close();
		}

		if (rank == 0)
		{
			std::ofstream interface("results/interface.txt");
			interface.precision(16);
			for (Node* const& fluidNode : fluid_->interfaceNodes_)
			{
				Node* solidNode = fluidNode->getInterfaceNode();
				for (int i = 0; i < dimension_; i++)
				{
					interface << solidNode->getDegreeOfFreedom(i)->getIndex() << "\n";
					interface << "Position: ";
					interface << std::fixed << solidNode->getDegreeOfFreedom(i)->getCurrentValue() << " " << fluidNode->getDegreeOfFreedom(i)->getCurrentValue() << "\n";
					interface << "Velocity: ";
					interface << std::fixed << solidNode->getDegreeOfFreedom(i)->getCurrentFirstTimeDerivative() << " " << fluidNode->getDegreeOfFreedom(i)->getCurrentFirstTimeDerivative() << "\n";
					interface << "Acceleration: ";
					interface << std::fixed << solidNode->getDegreeOfFreedom(i)->getCurrentSecondTimeDerivative() << " " << fluidNode->getDegreeOfFreedom(i)->getCurrentSecondTimeDerivative() << "\n";
				}
				interface << "\n";
			}
			interface << "\n";
			interface.close();
		}

		//export results to paraview
		if (rank == 0 && ((timeStep + 1) % solid_->parameters_->getExportFrequency() == 0))
		{
			solid_->computeCauchyStress();
            //fluid_->computeCauchyStress();
			solid_->exportToParaview(timeStep+1);
            fluid_->exportToParaview(timeStep+1);
		}
		delete[] constrainedDOFs;
	}
	delete[] externalForces;
	KSPDestroy(&ksp);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	PetscPrintf(PETSC_COMM_WORLD, "Monolithic Coupled PFEM Analysis Done. Elapsed time: %f\n", elapsed.count());
}

void CoupledDomain::solvePartitionedCoupledProblem()
{

}

double CoupledDomain::getInitialPositionNorm() const
{
    double initialNorm = 0.0;
    for (Node*& node : solid_->nodes_)
        for (int i = 0; i < dimension_; i++)
        {
            double value = node->getDegreeOfFreedom(i)->getInitialValue();
            initialNorm += value * value;
        }
    for (Node*& node : fluid_->nodes_)
        for (int i = 0; i < dimension_; i++)
            if (!node->isIsolated() && !node->isInterface())
            {
                double value = node->getDegreeOfFreedom(i)->getInitialValue();
                initialNorm += value * value;
            }
    return sqrt(initialNorm);
}

double CoupledDomain::computePositionNorm() const
{
    double norm = 0.0;
    for (Node*& node : solid_->nodes_)
        for (int i = 0; i < dimension_; i++)
        {
            double value = node->getDegreeOfFreedom(i)->getCurrentValue();
            norm += value * value;
        }
    for (Node*& node : fluid_->nodes_)
        for (int i = 0; i < dimension_; i++)
            if (!node->isIsolated() && !node->isInterface())
            {
                double value = node->getDegreeOfFreedom(i)->getCurrentValue();
                norm += value * value;
            }
    norm = sqrt(norm);

	if (norm < 1.0e-12)
		norm = 1.0;
	
	return norm;
}

double CoupledDomain::computePressureNorm() const
{
    double norm = 0.0;
        
    for (Node*& node : fluid_->nodes_)
    {
        double value = node->getDegreeOfFreedom(dimension_)->getCurrentValue();
        norm += value * value;
    }
    norm = sqrt(norm);

	if (norm < 1.0e-12)
		norm = 1.0;
	
	return norm;
}

void CoupledDomain::getConstrainedDOFs(int& ndofs, int*& constrainedDOFs)
{
    int nSolidDOFs; int* constrainedSolidDOFs;
	solid_->getConstrainedDOFs(nSolidDOFs, constrainedSolidDOFs);
    
    int nFluidDOFs; int* constrainedFluidDOFs;
	fluid_->getConstrainedDOFs(nFluidDOFs, constrainedFluidDOFs);

    ndofs = nSolidDOFs + nFluidDOFs;
    constrainedDOFs = new int[ndofs];
    for (int i = 0; i < nSolidDOFs; i++)
        constrainedDOFs[i] = constrainedSolidDOFs[i];
    for (int i = 0; i < nFluidDOFs; i++)
        constrainedDOFs[nSolidDOFs + i] = constrainedFluidDOFs[i];
    
    delete[] constrainedSolidDOFs; 
    delete[] constrainedFluidDOFs;
}

void CoupledDomain::getExternalForces(int& ndofs, std::vector<DegreeOfFreedom*>& dofsForces, double*& externalForces)
{
    int nSolidDOFs; std::vector<DegreeOfFreedom*> solidDofsForces; double* solidExternalForces;
	solid_->getExternalForces(nSolidDOFs, solidDofsForces, solidExternalForces);

    int nFluidDOFs; std::vector<DegreeOfFreedom*> fluidDofsForces; double* fluidExternalForces;
	fluid_->getExternalForces(nFluidDOFs, fluidDofsForces, fluidExternalForces);

    ndofs = nSolidDOFs + nFluidDOFs;
    dofsForces.reserve(ndofs);
    externalForces = new double[ndofs];
    for (int i = 0; i < nSolidDOFs; i++)
    {
        dofsForces.push_back(solidDofsForces[i]);
        externalForces[i] = solidExternalForces[i];
    }
    for (int i = 0; i < nFluidDOFs; i++)
    {
        dofsForces.push_back(fluidDofsForces[i]);
        externalForces[nSolidDOFs + i] = fluidExternalForces[i];
    }

    delete[] solidExternalForces;
    delete[] fluidExternalForces;
}

void CoupledDomain::applyNeummanConditions(Vec& vec, int& ndofs, const std::vector<DegreeOfFreedom*>& dofsForces, double*& externalForces, const double& loadFactor)
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

void CoupledDomain::assembleMonolithicLinearSystem(Mat& mat, Vec& vec)
{
    auto start_timer = std::chrono::high_resolution_clock::now();

    int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);


    //Computing solid elements contribution
    for (Element* const& el : solid_->elements_)
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

	//Computing fluid elements contribution
    for (Element* const& el : fluid_->elements_)
    {
		const std::vector<Node*>& nodes = el->getNodes();
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

	PetscPrintf(PETSC_COMM_WORLD, "Assemble Monolithic Linear System. Elapsed time: %f\n", elapsed.count());
}

void CoupledDomain::solveLinearSystem(KSP& ksp, Mat& mat, Vec& rhs, Vec& solution)
{
	auto start_timer = std::chrono::high_resolution_clock::now();

	KSPReset(ksp);
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
			//PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
		}
	}
	else
	{
		//PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
	}
	KSPSolve(ksp, rhs, solution);
	KSPConvergedReason reason;
	KSPGetConvergedReason(ksp,&reason);
	PetscInt nit;
	KSPGetIterationNumber(ksp, &nit);
	PetscReal norm;
	KSPGetResidualNorm(ksp, &norm);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	if (reason > 0)
	{
		PetscPrintf(PETSC_COMM_WORLD, "Solver converged within %d iterations. Elapsed time: %f\n", nit, elapsed.count());
	}
	else
	{
		if (reason == -3)
			PetscPrintf(PETSC_COMM_WORLD, "Solver convergence is very slow. Modifying the solver in order to improve the convergence...\n");
		else
			PetscPrintf(PETSC_COMM_WORLD, "Solver diverged, reason %d. Modifying the solver in order to improve the convergence...\n", reason);
		KSP ksp2;
		KSPCreate(PETSC_COMM_WORLD, &ksp2);
		KSPSetType(ksp2, KSPPREONLY);
		KSPSetTolerances(ksp2, 1.0e-8, PETSC_DEFAULT, PETSC_DEFAULT, 5000);
		KSPGMRESSetRestart(ksp2, 30);
		PetscInt restart;
		KSPGMRESGetRestart(ksp2, &restart);
		PC pc2;
		KSPGetPC(ksp2, &pc2);
		PCSetType(pc2, PCLU);
		KSPSetOperators(ksp2, mat, mat);
		PetscObjectTypeCompare((PetscObject)pc2, PCBJACOBI, &isbjacobi);
		if (isbjacobi)
		{
			PetscInt nlocal;
			KSP *subksp;
			PC subpc;
			KSPSetUp(ksp2);
			PCBJacobiGetSubKSP(pc2, &nlocal, NULL, &subksp);
			for (int i = 0; i < nlocal; i++)
			{
				KSPGetPC(subksp[i], &subpc);
				//PCFactorSetShiftType(subpc, MAT_SHIFT_NONZERO);
			}
		}
		else
		{
			//PCFactorSetShiftType(pc2, MAT_SHIFT_NONZERO);
		}
		KSPSolve(ksp2, rhs, solution);
		KSPGetConvergedReason(ksp2, &reason);
		KSPGetIterationNumber(ksp2, &nit);
		if (reason > 0)
		{
			PetscPrintf(PETSC_COMM_WORLD, "Solver converged within %d iterations.\n", nit);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "Changing the solver did not improve the convergence.\n");
		}
		KSPDestroy(&ksp2);
	}
}

void CoupledDomain::updateVariables(Vec& solution, double& positionNorm, double& pressureNorm)
{
	Vec All;
	VecScatter ctx;

	//Gathers the solution vector to the master process
	VecScatterCreateToAll(solution, &ctx, &All);
	VecScatterBegin(ctx, solution, All, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(ctx, solution, All, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterDestroy(&ctx);

	positionNorm = 0.0;
	pressureNorm = 0.0;
	//Updates solid nodal variables
    for (Node* const& node : solid_->nodes_)
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

    //Updates fluid nodal variables
    for (Node* const& node : fluid_->nodes_)
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
		else if (node->isIsolated() && node->isInterface())
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
		}
	}
    positionNorm = sqrt(positionNorm);
    pressureNorm = sqrt(pressureNorm);
	VecDestroy(&All);
}

void CoupledDomain::identifyMonolithicCoupledDOFs()
{
    for (Node*& solidNode : solid_->interfaceNodes_)
    {
        for (Node*& fluidNode : fluid_->interfaceNodes_)
        {
            bool coincidentNodes = true;
            for (int i = 0; i < dimension_; i++)
            {
                if (abs(solidNode->getDegreeOfFreedom(i)->getInitialValue() - fluidNode->getDegreeOfFreedom(i)->getInitialValue()) > 1.0e-10)
                    coincidentNodes = false;
            }
            if (coincidentNodes)
            {
                numberOfCoupledDOFs_ += solidNode->getNumberOfDegreesOfFreedom();
                numberOfCoupledNodes_++;
                solidNode->setInterfaceNode(fluidNode);
                fluidNode->setInterfaceNode(solidNode);
				for (int i = 0; i < dimension_; i++) //this is to avoid round errors between the solid and fluid nodes on the interface
            	{
					double coord = solidNode->getDegreeOfFreedom(i)->getInitialValue();
					fluidNode->getDegreeOfFreedom(i)->setInitialValue(coord);
					fluidNode->getDegreeOfFreedom(i)->setCurrentValue(coord);
            	}
                break;
            }
        }
    }

    pressureCoupled_ = numberOfCoupledDOFs_ - dimension_ * numberOfCoupledNodes_;

    numberOfDOFs_ = solid_->numberOfDOFs_ + fluid_->numberOfDOFs_ - numberOfCoupledDOFs_;

    PetscPrintf(PETSC_COMM_WORLD,
                "Solving FSI problem in a Monolithic fashion.\nSolid DOFs: %d\nFluid DOFs: %d\nCoupled DOFs: %d\nTotal DOFs: %d\n",
                solid_->numberOfDOFs_, fluid_->numberOfDOFs_, numberOfCoupledDOFs_, numberOfDOFs_);    
}

void CoupledDomain::identifyPartitionedCoupledDOFs()
{
    int numberOfCoincidentNodes = 0;
    for (Node*& solidNode : solid_->interfaceNodes_)
    {
        for (Node*& fluidNode : fluid_->interfaceNodes_)
        {
            bool coincidentNodes = true;
            for (int i = 0; i < dimension_; i++)
            {
                if (abs(solidNode->getDegreeOfFreedom(i)->getInitialValue() - fluidNode->getDegreeOfFreedom(i)->getInitialValue()) > 1.0e-10)
                    coincidentNodes = false;
            }
            if (coincidentNodes)
            {
                numberOfCoincidentNodes++;
                solidNode->setInterfaceNode(fluidNode);
                fluidNode->setInterfaceNode(solidNode);
            }
        }
    }
    numberOfCoupledDOFs_ = numberOfCoincidentNodes * dimension_;
    numberOfDOFs_ = solid_->numberOfDOFs_ + fluid_->numberOfDOFs_;

    for (Node*& fluidNode : fluid_->interfaceNodes_)
    {
        fluidNode->setConstrain(true);
        for (int i = 0; i < dimension_; i++)
		{
			DegreeOfFreedom* dof = fluidNode->getDegreeOfFreedom(i);
			fluid_->dirichletBoundaryConditions_.emplace_back(new DirichletBoundaryCondition(fluidNode, dof, 0.0));
		}
    }

    int nconditions = 0;
    for (auto& ibc : solid_->geometry_->getInterfaceBoundaryConditions())
	{
		if (ibc->getPoint()) nconditions++;
		else if (ibc->getLine()) nconditions += ibc->getLine()->getNumberOfBaseElements();
		else if (ibc->getSurface()) nconditions += ibc->getSurface()->getNumberOfBaseElements();
		else if (ibc->getVolume()) nconditions += ibc->getVolume()->getNumberOfBaseElements();
	}

    interfaceNeumannConditions_.reserve(nconditions);
    int index = -1;
    for (auto& ibc : solid_->geometry_->getInterfaceBoundaryConditions())
	{
        if (ibc->getLine())
		{
            Line* l = ibc->getLine();
			const std::vector<BaseLineElement*>& elements = l->getBaseElements();
			for (auto& el : elements)
			{
				int ndofs = dimension_ * el->getNumberOfNodes();
				interfaceNeumannConditions_.emplace_back(new InterfaceLineLoad(++index, ndofs, el));
			}
        }
        if (ibc->getSurface())
		{
            Surface* s = ibc->getSurface();
			const std::vector<BaseSurfaceElement*>& elements = s->getBaseElements();
			for (auto& el : elements)
			{
				int ndofs = dimension_ * el->getNumberOfNodes();
				interfaceNeumannConditions_.emplace_back(new InterfaceSurfaceLoad(++index, ndofs, el));
			}
        }
    }
    
    PetscPrintf(PETSC_COMM_WORLD,
                "Solving FSI problem in a Partitioned fashion.\nSolid DOFs: %d\nFluid DOFs: %d\nCoupled DOFs: %d\nTotal DOFs: %d\n",
                solid_->numberOfDOFs_, fluid_->numberOfDOFs_, numberOfCoupledDOFs_, numberOfDOFs_);    
}

void CoupledDomain::solvePartitionedCoupledPFEMProblem()
{
    // solid_->setReferenceConfiguration(ReferenceConfiguration::INITIAL);

    // identifyPartitionedCoupledDOFs();

    // if (fluid_->parameters_->getInitialAccel())
    // {
    //     solid_->computeInitialAccel();
    //     fluid_->computeInitialAccel();
    // }
	
	// int rank;
	// MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    // if (rank == 0)
    // {
    //     solid_->exportToParaview(0);
    //     fluid_->exportToParaview(0);
    // }

    // PetscLogDouble    bytes;
	// PetscMemoryGetCurrentUsage(&bytes);
	// PetscPrintf(PETSC_COMM_SELF,
    //             "Memory used by each processor to store problem data: %f Mb\n",
    //             bytes/(1024*1024));

    // fluid_->identifyFreeSurfaces();
	// fluid_->computeNodesCloudArea();

    // const int numberOfSteps = fluid_->parameters_->getNumberOfSteps();
	// const int maxNonlinearIterations = fluid_->parameters_->getMaxNonlinearIterations();
	// const double nonlinearTolerance = fluid_->parameters_->getNonlinearTolerance();
	
    // double* interfaceResidual = new double[numberOfCoupledDOFs_];

	// for (int timeStep = 0; timeStep < numberOfSteps; timeStep++)
	// {
	// 	PetscPrintf(PETSC_COMM_WORLD,
    //                 "\n----------------------- TIME STEP = %d, time = %f  -----------------------\n\n",
    //                 timeStep+1, (double)(timeStep+1)*fluid_->parameters_->getDeltat());
		
	// 	solid_->setPastVariables();
    //     fluid_->setPastVariables();

	// 	bool converged = false;   int iteration = 0; 
    //     while (!converged && iteration < 3)
    //     {  
    //         PetscPrintf(PETSC_COMM_WORLD, "Gauss-Seidel iteration: %d\n", iteration);
    //         transferSolidDisplacements();
    //         fluid_->solveStaggeredPFEMProblem();
    //         int ndofsInterfaceForces; std::vector<DegreeOfFreedom*> dofsForces; double* interfaceForces;
    //         transferFluidForces(ndofsInterfaceForces, dofsForces, interfaceForces);
    //         solid_->solveStaggeredProblem(ndofsInterfaceForces, dofsForces, interfaceForces);
    //         delete[] interfaceForces;

           
    //         int index = -1;
    //         for (auto& solidNode : solid_->interfaceNodes_)
    //         {
    //             index++;
    //             Node* fluidNode = solidNode->getInterfaceNode();
    //             for (int i = 0; i < dimension_; i++)
    //             {
    //                 interfaceResidual[dimension_ * index + i] = solidNode->getDegreeOfFreedom(i)->getCurrentValue() - 
    //                                                             fluidNode->getDegreeOfFreedom(i)->getCurrentValue();
    //             }
    //         }
    //         double interfaceNorm = 0.0;
    //         for (int i = 0; i < numberOfCoupledDOFs_; i++)
    //             interfaceNorm += interfaceResidual[i] * interfaceResidual[i];
    //         interfaceNorm = sqrt(interfaceNorm);
    //         if (interfaceNorm < 1.0e-6)
    //             converged = true;
    //         PetscPrintf(PETSC_COMM_WORLD, "Gauss-Seidel norm: %E\n", interfaceNorm);
    //         iteration++;
    //     }

	// 	//export results to paraview
	// 	if (rank == 0)
	// 	{
	// 		solid_->computeCauchyStress();
    //         fluid_->computeCauchyStress();
	// 		solid_->exportToParaview(timeStep+1);
    //         fluid_->exportToParaview(timeStep+1);
	// 	}
    //     fluid_->checkMeshQuality();
	// 	fluid_->computeNodesCloudArea();
	// 	fluid_->remesh();
	// 	fluid_->alphaShape();
	// 	fluid_->nodalNeighborSearch();
	// 	fluid_->reorderDOFs();
	// 	fluid_->identifyFreeSurfaces();
	// }
    // delete[] interfaceResidual;
}

void CoupledDomain::transferSolidDisplacements()
{
    for (Node*& solidNode : solid_->interfaceNodes_)
    {
        Node* fluidNode = solidNode->getInterfaceNode();
        for (int i = 0; i < dimension_; i++)
        {
            double val = solidNode->getDegreeOfFreedom(i)->getCurrentValue();
            fluidNode->getDegreeOfFreedom(i)->setCurrentValue(val);
        }
    }
}

void CoupledDomain::transferFluidForces(int& ndofs,
                                        std::vector<DegreeOfFreedom*>& dofs,
                                        double*& interfaceForces)
{
    ndofs = 0;
	for (InterfaceNeumannBoundaryCondition* const& ibc : interfaceNeumannConditions_)
	{
		ndofs += ibc->getNumberOfDOFs();
	}
	dofs.reserve(ndofs);
	interfaceForces = new double[ndofs];
	int aux = -1;
	for (InterfaceNeumannBoundaryCondition* const& ibc : interfaceNeumannConditions_)
	{
		int ndof = ibc->getNumberOfDOFs();
		std::vector<DegreeOfFreedom*> ibc_dofs; 
		double *val;
		ibc->getNodalForce(dimension_, ibc_dofs, val);
		for (int i = 0; i < ndof; i++)
		{
			dofs.push_back(ibc_dofs[i]);
			interfaceForces[++aux] = val[i];
		}
		delete[] val;
	}
}

void CoupledDomain::createSystemMatrix(Mat& mat)
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
		for (Node* const& solidNode : solid_->nodes_)
		{
			if (!solidNode->isIsolated())
			{
				int ndof = solidNode->getNumberOfDegreesOfFreedom();
				const std::vector<Node*>& neighborNodes = solidNode->getNeighborNodes();
				for (auto& neighborNode : neighborNodes)
				{
					nnz += ndof * neighborNode->getNumberOfDegreesOfFreedom();
				}
			}
		}
        for (Node* const& fluidNode : fluid_->nodes_)
		{
			if (!fluidNode->isIsolated())
			{
				int ndof = fluidNode->getNumberOfDegreesOfFreedom();
				const std::vector<Node*>& neighborNodes = fluidNode->getNeighborNodes();
				for (auto& neighborNode : neighborNodes)
				{
					nnz += ndof * neighborNode->getNumberOfDegreesOfFreedom();
				}
			}
		}
        //Exclude the duplicated non zeroes at the interface DOFs
        for (Node* const& fluidNode : fluid_->interfaceNodes_)
        {
			if (!fluidNode->isIsolated())
			{
				const std::vector<Node*>& neighborNodes = fluidNode->getNeighborNodes();
				for (auto& neighborNode : neighborNodes)
				{
					if (neighborNode->isInterface())
						nnz -= dimension_ * dimension_; //assuming that only positions are coupled
				}
			}
        }

		int* col = new int[nnz]; //stores the number of the column of each nonzero term
		int* ptr = new int[N+1]; //accumulates the number of nonzeroes per row in a CSR format
		ptr[0] = 0;

		//Iterating over the nodes in the permuted order. This block of code needs to acces the dofs in a crescent and consecutive order
		//As only the dofs of blocked nodes are contributing to the matrix, we access only the firsts nblocked nodes
		int sum = -1;
		for (unsigned int k = 0; k < solid_->numberOfBlockedNodes_; k++)
		{
			Node*& solidNode = solid_->nodes_[solid_->perm_[k]];
			const std::vector<Node*>& neighborNodes = solidNode->getNeighborNodes();
			const std::vector<DegreeOfFreedom*>& i_dofs = solidNode->getDegreesOfFreedom();
            Node* fluidNode = solidNode->getInterfaceNode();
			for (DegreeOfFreedom* const& i_dof : i_dofs)
			{
				int i = i_dof->getIndex();
				for (auto& neighborNode : neighborNodes)
				{
					const std::vector<DegreeOfFreedom*>& j_dofs = neighborNode->getDegreesOfFreedom();
					for (DegreeOfFreedom* const& j_dof : j_dofs)
					{
						int j = j_dof->getIndex();
						col[++sum] = j;
					}
				}
                if (solidNode->isInterface() && !fluidNode->isIsolated())
                {
                    const std::vector<Node*>& neighborNodesFluid = fluidNode->getNeighborNodes();
                    for (auto& neighborNode : neighborNodesFluid)
                    {
                        const std::vector<DegreeOfFreedom*>& j_dofs = neighborNode->getDegreesOfFreedom();
                        if (neighborNode->isInterface())
                        {
                            int j = j_dofs.back()->getIndex();
                            col[++sum] = j;
                        }
                        else
                        {
                            for (DegreeOfFreedom* const& j_dof : j_dofs)
                            {
                                int j = j_dof->getIndex();
                                col[++sum] = j;
                            }
                        }
                    }
                }
				ptr[i+1] = sum+1;
			}
		}

        for (unsigned int k = 0; k < fluid_->numberOfBlockedNodes_; k++)
		{
			Node*& fluidNode = fluid_->nodes_[fluid_->perm_[k]];
            const std::vector<Node*>& neighborNodes = fluidNode->getNeighborNodes();
			const std::vector<DegreeOfFreedom*>& i_dofs = fluidNode->getDegreesOfFreedom();
            Node* solidNode = fluidNode->getInterfaceNode();

            if (fluidNode->isInterface()) //only pressure degree of freedom needs to be computed
            {
                int i = i_dofs.back()->getIndex();
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
            else
            {
                for (DegreeOfFreedom* const& i_dof : i_dofs)
			    {
                    int i = i_dof->getIndex();
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
	//MatCreateAIJ(PETSC_COMM_WORLD, n, n, N, N, 0, d_nnz, 0, o_nnz, &mat);
	MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, N, N, 200, PETSC_NULL, 200, PETSC_NULL, &mat);

	delete[] d_nnz_all;
    delete[] o_nnz_all;

	MatSetFromOptions(mat);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;

	//PetscPrintf(PETSC_COMM_WORLD, "Allocating memory for the sparse matrix. Elapsed time: %f\n", elapsed.count());
}

void CoupledDomain::coupleInterfaceDOFs()
{
    /*
        Fluid interface DOFs receive the index of Solid interface DOFs.
        This way, we guarantee that the interface is always included in the system.
    */
    
    int dofCount = solid_->numberOfBlockedDOFs_ - 1;

    int nFluidNodes = fluid_->nodes_.size();

    for (int i = 0; i < nFluidNodes; i++)
    {
        Node* fluidNode = fluid_->nodes_[fluid_->perm_[i]];
        const std::vector<DegreeOfFreedom*>& fluidDOFs = fluidNode->getDegreesOfFreedom();

        if (fluidNode->isInterface())
        {
            Node* solidNode = fluidNode->getInterfaceNode();
            const std::vector<DegreeOfFreedom*>& solidDOFs = solidNode->getDegreesOfFreedom();
            for (int j = 0; j < dimension_; j++)
            {
                int id = solidDOFs[j]->getIndex();
                fluidDOFs[j]->setIndex(id);
            }
            fluidDOFs[dimension_]->setIndex(++dofCount);
        }
        else
        {
            for (DegreeOfFreedom* const& dof : fluidDOFs)
            {
                dof->setIndex(++dofCount);
            }

        }
    }
}

void CoupledDomain::domainDecomposition()
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

	unsigned int n_nodes_fluid = fluid_->nodes_.size();
    unsigned int n_elem_fluid = fluid_->elements_.size();

    if (fluid_->elementPartition_) delete[] fluid_->elementPartition_;
    if (fluid_->nodePartition_) delete[] fluid_->nodePartition_;
	fluid_->elementPartition_ = new idx_t[n_elem_fluid];
    fluid_->nodePartition_ = new idx_t[n_nodes_fluid];
    
    unsigned int n_nodes_solid = solid_->nodes_.size();
    unsigned int n_elem_solid = solid_->elements_.size();

    if (solid_->elementPartition_) delete[] solid_->elementPartition_;
    if (solid_->nodePartition_) delete[] solid_->nodePartition_;
	solid_->elementPartition_ = new idx_t[n_elem_solid];
    solid_->nodePartition_ = new idx_t[n_nodes_solid];

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
        for (unsigned int i = 0; i < solid_->numberOfBlockedNodes_; i++) //Looping over the permuted order of nodes
		{
			Node* solidNode = solid_->nodes_[solid_->perm_[i]];
			int ndof = solidNode->getNumberOfDegreesOfFreedom();
            for (int j = 0; j < size; j++)
			{
                if (((ndof*i) >= start[j]) && ((ndof*i) <= end[j]))
				{
					solidNode->setRank(j);
					solid_->nodePartition_[solid_->perm_[i]] = j;
				}
            }
        }
        for (unsigned int i = 0; i < fluid_->numberOfBlockedNodes_; i++) //Looping over the permuted order of nodes
		{
			Node* fluidNode = fluid_->nodes_[fluid_->perm_[i]];
            if (fluidNode->isInterface())
            {
                int rank = fluidNode->getInterfaceNode()->getRank();
                fluidNode->setRank(rank);
				fluid_->nodePartition_[fluid_->perm_[i]] = rank;
            }
            else
            {
                int ndof = fluidNode->getNumberOfDegreesOfFreedom();
                for (int j = 0; j < size; j++)
                {
                    if (((ndof*i + solid_->numberOfBlockedDOFs_) >= start[j]) && ((ndof*i + solid_->numberOfBlockedDOFs_) <= end[j]))
                    {
                        fluidNode->setRank(j);
                        fluid_->nodePartition_[fluid_->perm_[i]] = j;
                    }
                }

            }
        }

        for (unsigned int i = 0; i < n_elem_solid; i++)
		{
            Element* elem = solid_->elements_[i];
			Node* last_node = elem->getNodes().back(); //could be any nodes of the element
            int last_dof_index = last_node->getDegreesOfFreedom().back()->getIndex();

            for (int j = 0; j < size; j++)
			{
                if ((last_dof_index >= start[j]) && (last_dof_index <= end[j]))
				{
                    elem->setRank(j);
					solid_->elementPartition_[i] = j;
                    break;
                }
            }
        }

        for (unsigned int i = 0; i < n_elem_fluid; i++)
		{
            Element* elem = fluid_->elements_[i];
			Node* last_node = elem->getNodes().back(); //could be any nodes of the element
            int last_dof_index = last_node->getDegreesOfFreedom().back()->getIndex();

            for (int j = 0; j < size; j++)
			{
                if ((last_dof_index >= start[j]) && (last_dof_index <= end[j]))
				{
                    elem->setRank(j);
					fluid_->elementPartition_[i] = j;
                    break;
                }
            }
        }
    }

	//defining the ownership range of isolated nodes
	N = numberOfNodes_ - numberOfBlockedNodes_;
    if (N > 0)
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
		for (unsigned int i = solid_->numberOfBlockedNodes_; i < n_nodes_solid; i++)
		{
			Node* solidNode = solid_->nodes_[solid_->perm_[i]];
			int index = i - solid_->numberOfBlockedNodes_;
			for (int j = 0; j < size; j++)
			{
				if (index >= start[j] && index <= end[j])
				{
					solidNode->setRank(j);
					solid_->nodePartition_[solid_->perm_[i]] = j;
				}
			}
		}

        for (unsigned int i = fluid_->numberOfBlockedNodes_; i < n_nodes_fluid; i++)
		{
			Node* fluidNode = fluid_->nodes_[fluid_->perm_[i]];
			int index = i - fluid_->numberOfBlockedNodes_ + n_nodes_solid - solid_->numberOfBlockedNodes_;
			for (int j = 0; j < size; j++)
			{
				if (index >= start[j] && index <= end[j])
				{
					fluidNode->setRank(j);
					fluid_->nodePartition_[fluid_->perm_[i]] = j;
				}
			}
		}
	}

	MPI_Bcast(solid_->elementPartition_,n_elem_solid,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(solid_->nodePartition_,n_nodes_solid,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(fluid_->elementPartition_,n_elem_fluid,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(fluid_->nodePartition_,n_nodes_fluid,MPI_INT,0,PETSC_COMM_WORLD);

	if (rank != 0)
	{
		for (unsigned int i = 0; i < n_nodes_solid; i++)
		{
			solid_->nodes_[i]->setRank(solid_->nodePartition_[i]);
		}
		for (unsigned int i = 0; i < n_elem_solid; i++)
		{
			solid_->elements_[i]->setRank(solid_->elementPartition_[i]);
		}
        for (unsigned int i = 0; i < n_nodes_fluid; i++)
		{
			fluid_->nodes_[i]->setRank(fluid_->nodePartition_[i]);
		}
		for (unsigned int i = 0; i < n_elem_fluid; i++)
		{
			fluid_->elements_[i]->setRank(fluid_->elementPartition_[i]);
		}
	}

	MPI_Barrier(PETSC_COMM_WORLD);

	auto end_timer = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end_timer - start_timer;
}