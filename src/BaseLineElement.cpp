#include "BaseLineElement.h"

BaseLineElement::BaseLineElement(const int index,
                                 ParametricLineElement& parametricElement,     
                                 const std::vector<Node*>& nodes)
        : BaseElement(index, parametricElement, nodes) {}

BaseLineElement::~BaseLineElement() {}

double BaseLineElement::getRadius() const
{
        return 0.0;
}

double BaseLineElement::getJacobianIntegration() const
{
        return 0.0;
}