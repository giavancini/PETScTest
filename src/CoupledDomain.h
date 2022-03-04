#pragma once

#include "FluidDomain.h"
#include "SolidDomain.h"

class CoupledDomain
{
    public:
	CoupledDomain(FluidDomain* fluid, SolidDomain* solid);

	~CoupledDomain();

	void setNumberOfSteps(const int numberOfSteps);

    void setMaxNonlinearIterations(const int maxNonlinearIterations);

    void setNonlinearTolerance(const double nonlinearTolerance);

    void setDeltat(const double deltat);
    
    void setSpectralRadius(const double rhoInf);

    void setGeneralizedAlphas(const double alphaM, const double alphaF);

    void setInitialAccel(const bool initialAccel);

    void setExportFrequency(const int& freq);

	void solveMonolithicCoupledProblem();

	void solveMonolithicCoupledPFEMProblem();

	void solvePartitionedCoupledProblem();

    void solvePartitionedCoupledPFEMProblem();

    private:

    double getInitialPositionNorm() const;

    double computePositionNorm() const;

    double computePressureNorm() const;

    void getConstrainedDOFs(int& ndofs, int*& constrainedDOFs);

    void getExternalForces(int& ndofs, std::vector<DegreeOfFreedom*>& dofs, double*& externalForces);

    void applyNeummanConditions(Vec& vec, int& ndofs, const std::vector<DegreeOfFreedom*>& dofsForces, double*& externalForces, const double& loadFactor);

    void assembleMonolithicLinearSystem(Mat& mat, Vec& vec);

	void assembleMonolithicPFEMLinearSystem(Mat& mat, Vec& vec);

    void solveLinearSystem(KSP& ksp, Mat& mat, Vec& rhs, Vec& solution);

    void updateVariables(Vec& solution, double& positionNorm, double& pressureNorm);

	void identifyMonolithicCoupledDOFs();

    void identifyPartitionedCoupledDOFs();

    void transferSolidDisplacements();

    void transferFluidForces(int& ndofs,
                             std::vector<DegreeOfFreedom*>& dofs,
                             double*& interfaceForces);
    
    void createSystemMatrix(Mat& mat);

	void coupleInterfaceDOFs();

	void domainDecomposition();

    private:
    FluidDomain* fluid_;
    SolidDomain* solid_;
	int dimension_;
	int numberOfDOFs_;
	int numberOfBlockedDOFs_;
	int numberOfIsolatedDOFs_;
	int numberOfCoupledDOFs_;
	int numberOfNodes_;
	int numberOfBlockedNodes_;
	int numberOfIsolatedNodes_;
    int numberOfCoupledNodes_;
    int numberOfCoupledBlockedDOFs_;
	bool pressureCoupled_;
    std::vector<InterfaceNeumannBoundaryCondition*> interfaceNeumannConditions_;
};



