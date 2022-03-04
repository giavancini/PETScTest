#pragma once
#include "BaseElement.h"

class BaseLineElement : public BaseElement
{
    public:
    BaseLineElement(const int index,
                    ParametricLineElement& parametricElement,     
                    const std::vector<Node*>& nodes);

    ~BaseLineElement() override;

    double getRadius() const override;

    double getJacobianIntegration() const override;
};