#pragma once

#include "Element.h"
#include "Material.h"
#include "AnalysisParameters.h"

class VolumeElement : public Element
{
    public:
    VolumeElement(const int& index,
                 const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                 Material* material,
                 BaseVolumeElement* base,
                 AnalysisParameters* params);

    ~VolumeElement() override;

    void setMaterial(Material* material);

    void setBaseElement(BaseVolumeElement* base);

    void setAnalysisParameters(AnalysisParameters* params);

    Material* getMaterial() const;

    BaseVolumeElement* getBaseElement() const override;

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
                                           double A0[3][3]) const;

    inline void getCurrentJacobianMatrix(double** dphi_dxsi,
                                         double A1[3][3]) const;

    inline void getCurrentJacobianMatrixTimeDerivative(double** dphi_dxsi,
                                                       double dA1_dt[3][3]) const;

    inline double getMatrixDeterminant(const double matrix[3][3]) const;
        
    inline void getInverseMatrix(const double matrix[3][3],
                                 const double& determinant,
                                 double inverse[3][3]) const;

    inline void getReferenceShapeFunctionsGradient(double** dphi_dxsi,
                                                   double dxsi_dx[3][3], double dphi_dx[][3]) const;

    inline void getCurrentShapeFunctionsGradient(double** dphi_dxsi,
                                                 double dxsi_dy[3][3], double dphi_dy[][3]) const;

    inline void getDeformationGradient(double dphi_dx[][3],
                                       double A[3][3]) const;

    inline void getDeformationGradientTimeDerivative(double dphi_dx[][3],
                                                     double dA_dt[3][3]) const;

    inline void getStrainTensor(const double A[3][3],
                                double E[6]) const;

    inline void getStrainTensorTimeDerivative(const double A[3][3],
                                              const double dA_dt[3][3],
                                              double dE_dt[6]) const;

    inline void getInverseRightCauchyTensor(const double AI[3][3],
                                            double CI[6]) const;

    inline void getStrainTensorFirstDerivative(const double A[3][3],
                                               double dphi_dx[][3],
                                               double dE_dy[][6]) const;

    inline void getStrainTensorTimeDerivativeFirstDerivative(const double A[3][3],
                                                             const double dA_dt[3][3],
                                                             double dphi_dx[][3],
                                                             double dE_dtdy[][6]) const;

    inline void getStrainTensorSecondDerivative(const int i,
                                                const int j,
                                                const double dphi_dx[][3],
                                                double d2E_dydy[6]) const;

    inline double doubleContraction(const double v1[6],
                                    const double v2[6]) const;

    private:
    Material* material_;
    BaseVolumeElement* base_;
    AnalysisParameters* parameters_;
};