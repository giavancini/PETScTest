#pragma once

#include "Element.h"
#include "Material.h"
#include "AnalysisParameters.h"

class PlaneElement : public Element
{
    public:
    PlaneElement(const int& index,
                 const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                 Material* material,
                 BaseSurfaceElement* base,
                 AnalysisParameters* params);

    ~PlaneElement() override;

    void setMaterial(Material* material);

    void setBaseElement(BaseSurfaceElement* base);

    void setAnalysisParameters(AnalysisParameters* params);

    Material* getMaterial() const;

    BaseSurfaceElement* getBaseElement() const override;

    ParametricElement* getParametricElement() const override;
    
    const std::vector<Node*>& getNodes() const override;
        
    void getDOFIndexes(unsigned int& ndof,
                       int*& indexes) const override;

    void getCauchyStress(double**& nodalCauchyStress) const override;

    void elementContributions(int& ndofs1,
                              int& ndofs2,
                              int*& indexes,
                              double*& rhsValues,
                              double*& hessianValues) const override;

    void clearNeighborElements() override;

    inline void getReferenceJacobianMatrix(double** dphi_dxsi,
                                           double A0[2][2]) const;

    inline void getCurrentJacobianMatrix(double** dphi_dxsi,
                                         double A1[2][2]) const;

    inline void getCurrentJacobianMatrixTimeDerivative(double** dphi_dxsi,
                                                       double dA1_dt[2][2]) const;

    inline double getMatrixDeterminant(const double matrix[2][2]) const;
        
    inline void getInverseMatrix(const double matrix[2][2],
                                 const double& determinant,
                                 double inverse[2][2]) const;

    inline void getReferenceShapeFunctionsGradient(double** dphi_dxsi,
                                                   double dxsi_dx[2][2], double dphi_dx[][2]) const;

    inline void getCurrentShapeFunctionsGradient(double** dphi_dxsi,
                                                 double dxsi_dy[2][2], double dphi_dy[][2]) const;

    inline void getDeformationGradient(double dphi_dx[][2],
                                       double A[2][2]) const;

    inline void getDeformationGradientTimeDerivative(double dphi_dx[][2],
                                                     double dA_dt[2][2]) const;

    inline void getStrainTensor(const double A[2][2],
                                double E[3]) const;

    inline void getStrainTensorTimeDerivative(const double A[2][2],
                                              const double dA_dt[2][2],
                                              double dE_dt[3]) const;

    inline void getInverseRightCauchyTensor(const double AI[2][2],
                                            double CI[3]) const;

    inline void getStrainTensorFirstDerivative(const double A[2][2],
                                               double dphi_dx[][2],
                                               double dE_dy[][3]) const;

    inline void getStrainTensorTimeDerivativeFirstDerivative(const double A[2][2],
                                                             const double dA_dt[2][2],
                                                             double dphi_dx[][2],
                                                             double dE_dtdy[][3]) const;

    inline void getStrainTensorSecondDerivative(const int i,
                                                const int j,
                                                const double dphi_dx[][2],
                                                double d2E_dydy[3]) const;

    inline double doubleContraction(const double v1[3],
                                    const double v2[3]) const;

    inline double computeStabilizationParameter() const;

    private:
    Material* material_;
    BaseSurfaceElement* base_;
    AnalysisParameters* parameters_;
};