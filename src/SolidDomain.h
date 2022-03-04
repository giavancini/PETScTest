#pragma once

#include "mesh_interface/Geometry.h"
#include "mesh_interface/Mesh.h"
#include "AnalysisParameters.h"
#include "DirichletBoundaryCondition.h"
#include "NeumannBoundaryCondition.h"
#include "InterfaceBoundaryCondition.h"
#include "LineElement.h"
#include "PlaneElement.h"
#include "VolumeElement.h"
#include "TriangularMesher.h"
#include "TetrahedralMesher.h"
#include <unordered_map>
#include <petscksp.h>
#include <metis.h>
#include <chrono>

class SolidDomain
{
public:
	SolidDomain(Geometry* geometry, const int& index = 0);

	~SolidDomain();

	void setNumberOfSteps(const int numberOfSteps);

	void setMaxNonlinearIterations(const int maxNonlinearIterations);

	void setNonlinearTolerance(const double nonlinearTolerance);

	void setDeltat(const double deltat);

	void setGravity(const double gravity_x,
                    const double gravity_y,
                    const double gravity_z = 0.0);

	void setSpectralRadius(const double rhoInf);

	void setGeneralizedAlphas(const double alphaM, const double alphaF);

    void setAlpha(const double alpha);

    void setMeshLength(const double h);

    void setInitialAccel(const bool initialAccel);

	void setExportFrequency(const int& freq);

	void setLumpedMass(const bool& useLumpedMass);

	void setReferenceConfiguration(const ReferenceConfiguration reference);

	void applyMaterial(const std::vector<Line*> lines, Material*& material);

	void applyMaterial(const std::vector<Surface*> surfaces, Material*& material);

	void applyMaterial(const std::vector<Volume*> volumes, Material*& material);

	void generateMesh(const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm = AUTO, std::string geofile = std::string(), 
					  const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);

	void solveStaticProblem();

	void solveTransientProblem();

	void solveStaggeredProblem(int& ndofsInterfaceForces,
                               std::vector<DegreeOfFreedom*>& dofsInterfaceForces,
                               double*& interfaceForces);

	void setInitialVelocityX(std::function<double(double, double, double)> function);

	void setInitialVelocityY(std::function<double(double, double, double)> function);

	void setInitialVelocityZ(std::function<double(double, double, double)> function);

private:
	Node* getNode(const int& index);

	Element* getElement(const int& index);

	Material* getMaterial(const int& index);

	const std::vector<DirichletBoundaryCondition*>& getDirichletBoundaryConditions();

	const std::vector<NeumannBoundaryCondition*>& getNeumannBoundaryConditions();

	double getInitialPositionNorm() const;

	void setPastVariables();

	void computeCurrentVariables();

	void computeIntermediateVariables();

	void getConstrainedDOFs(int& ndofs, int*& constrainedDOFs);

	void getExternalForces(int& ndofs, std::vector<DegreeOfFreedom*>& dofs, double*& externalForces);

	void assembleStaticLinearSystem(Mat& mat, Vec& vec);

	void assembleTransientLinearSystem(Mat& mat, Vec& vec);

	void applyNeummanConditions(Vec& vec, int& ndofs, const std::vector<DegreeOfFreedom*>& dofsForces, double*& externalForces, const double& loadFactor);

	void solveLinearSystem(KSP& ksp, Mat& mat, Vec& rhs, Vec& solution);

	void updateVariables(Vec& solution, double& positionNorm, double& pressureNorm);

	void computeCauchyStress();

	void computeInitialAccel();

	void exportToParaview(const int& step);

	void readInput(const std::string& inputFile, const bool& deleteFiles, const PartitionOfUnity elementType);

	void transferGeometricBoundaryConditions();

	void domainDecompositionMETIS(const PartitionOfUnity& elementType);

	void nodalNeighborSearch() const;

	void elementalNeighborSearch() const;

	void createSystemMatrix(Mat& mat);

	void reorderDOFs();

	void domainDecomposition();
private:
	int index_;
	int dimension_;
	int numberOfDOFs_;
	int numberOfBlockedDOFs_;
	int numberOfIsolatedDOFs_;
	int numberOfNodes_;
	int numberOfBlockedNodes_;
	int numberOfIsolatedNodes_;
	AnalysisParameters* parameters_;
	Geometry * geometry_;
	Mesher* remesh_;
	std::vector<Node*> nodes_;
	std::vector<Node*> interfaceNodes_;
	std::vector<Element*> elements_;
	std::vector<Material*> materials_;
	std::vector<DirichletBoundaryCondition*> dirichletBoundaryConditions_;
	std::vector<NeumannBoundaryCondition*> neumannBoundaryConditions_;
	idx_t* elementPartition_;
	idx_t* nodePartition_;
	idx_t* perm_;

public:
	friend class CoupledDomain;
};

