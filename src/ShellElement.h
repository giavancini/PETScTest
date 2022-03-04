#pragma once

#include "Element.h"
#include "BaseSurfaceElement.h"
#include "Material.h"
#include "AnalysisParameters.h"

class ShellElement : public Element
{
    public:

    ShellElement(const int& index,
                 const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                 const double& thickness,
                 Material* material,
                 BaseSurfaceElement* base,
                 AnalysisParameters* params);

    ~ShellElement() override;

    void setMaterial(Material* material);

    void setBaseElement(BaseSurfaceElement* base);

    void setAnalysisParameters(AnalysisParameters* params);

    Material* getMaterial() const;

    BaseSurfaceElement* getBaseElement() const;
    
    const std::vector<Node*>& getNodes() const override;
        
    void getDOFIndexes(unsigned int& ndof,
                       int*& indexes) const override;

    void getCauchyStress(double**& nodalCauchyStress) const override;

    void elementContributions(int& ndofs1,
                              int& ndofs2,
                              int*& indexes,
                              double*& rhsValues,
                              double*& hessianValues) const override;

    void incrementNormalVector();

    inline void getSurfaceDependentVariables(double* phi,
                                             double** dphi_dxsi,
                                             double n[3],
                                             double g[3],
                                             double& t0,
                                             double& t,
                                             double dx_dxsi[3][2],
                                             double dy_dxsi[3][2],
                                             double dn_dxsi[3][2],
                                             double dg_dxsi[3][2],
                                             double dt0_dxsi[2],
                                             double dt_dxsi[2]) const;

    inline void getReferenceJacobianMatrix(const double t0,
                                           const double h0,
                                           const double hm,
                                           const double n[3],
                                           const double dx_dxsi[3][2],
                                           const double dn_dxsi[3][2],
                                           const double dt0_dxsi[2],
                                           double A0[3][3]) const;

    inline void getCurrentJacobianMatrix(const double t,
                                         const double h0,
                                         const double hm,
                                         const double g[3],
                                         const double dy_dxsi[3][2],
                                         const double dg_dxsi[3][2],
                                         const double dt_dxsi[2],
                                         double A1[3][3]) const;
    
    inline double getMatrixDeterminant(const double matrix[3][3]) const;
        
    inline void getInverseMatrix(const double matrix[3][3],
                                 const double& determinant,
                                 double inverse[3][3]) const;

    inline void getDeformationGradient(const double A1[3][3],
                                       const double A0I[3][3],
                                       double A[3][3]) const;

    inline void getStrainTensor(const double A[3][3],
                                double E[6]) const;

    inline void getFirstDerivatives(double* phi,
                                    double** dphi_dxsi,
                                    const double hm,
                                    const double h0,
                                    const double t,
                                    const double g[3],
                                    const double dt_dxsi[2],
                                    const double A0I[3][3],
                                    const double A[3][3],
                                    const double dg_dxsi[3][2],
                                    double dA_dy[][3][3],
                                    double dE_dy[][6]) const;

    inline void getStrainTensorSecondDerivative(double* phi,
                                                double** dphi_dxsi,
                                                const double h0,
                                                const double hm,
                                                const int i,
                                                const int j,
                                                const double A0I[3][3],
                                                const double A[3][3],
                                                const double dA_dy[][3][3],
                                                double d2E_dydy[6]) const;

    inline double doubleContraction(const double v1[6],
                                    const double v2[6]) const;

    private:
    double thickness_;
    Material* material_;
    BaseSurfaceElement* base_;
    AnalysisParameters* parameters_;
    int numberOfThicknessIntegrationPoints_;
};