#pragma once
#include "BaseElement.h"

class BaseSurfaceElement : public BaseElement
{
    public:
    BaseSurfaceElement(const int index,
                       ParametricSurfaceElement& parametricElement,     
                       const std::vector<Node*>& nodes);
    
    ~BaseSurfaceElement() override;

    double getRadius() const override;

    double getJacobianIntegration() const override;

    void getInitialNormalVector(double* xsi,
                                double normal[3]) const;
        
    inline void getCurrentJacobianMatrix(double** dphi_dxsi,
                                         double dy_dxsi[2][2]) const;

    inline double getMatrixDeterminant(const double jacobianMatrix[2][2]) const;
};