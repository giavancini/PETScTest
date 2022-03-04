#pragma once
#include "BaseElement.h"

class BaseVolumeElement : public BaseElement
{
    public:
    BaseVolumeElement(const int index,
                       ParametricVolumeElement& parametricElement,     
                       const std::vector<Node*>& nodes);
    
    ~BaseVolumeElement() override;

    double getRadius() const override;

    double getJacobianIntegration() const override;

    inline void getCurrentJacobianMatrix(double** dphi_dxsi,
                                         double dy_dxsi[3][3]) const;

    inline double getMatrixDeterminant(const double jacobianMatrix[3][3]) const;        
};