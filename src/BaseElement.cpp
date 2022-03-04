#include "BaseElement.h"

BaseElement::BaseElement(const int index,
                         ParametricElement& parametricElement,     
                         const std::vector<Node*>& nodes)
        : index_(index),
          parametricElement_(&parametricElement),
          nodes_(nodes),
          plot_(false) {}

int BaseElement::getIndex() const
{
    return index_;
}

int BaseElement::getNumberOfNodes() const
{
    return nodes_.size();
}

ParametricElement* BaseElement::getParametricElement() const
{
    return parametricElement_;
}

const std::vector<Node*>& BaseElement::getNodes() const
{
    return nodes_;
}

bool BaseElement::getPlot() const
{
    return plot_;
}

void BaseElement::setPlot(const bool plot)
{
    plot_ = plot;
}

void BaseElement::setIndex(const int index)
{
    index_ = index;
}