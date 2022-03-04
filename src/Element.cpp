#include "Element.h"

Element::Element(const int index,
                 const std::vector<DegreeOfFreedom*>& degreesOfFreedom)
        : index_(index),
          rank_(0),
          isBoundary_(false),
          isActive_(true),
          referenceConfiguration_(ReferenceConfiguration::PAST),
          degreesOfFreedom_(degreesOfFreedom) {}

void Element::setIndex(const int& index)
{
    index_ = index;
}

void Element::setRank(const int& rank)
{
    rank_ = rank;
}

void Element::setBoundary(const bool& isBoundary)
{
    isBoundary_ = isBoundary;
}

void Element::setActive(const bool& isActive)
{
    isActive_ = isActive;
}

void Element::setReferenceConfiguration(const ReferenceConfiguration reference)
{
    referenceConfiguration_ = reference;
}

void Element::setDegreesOfFreedom(const std::vector<DegreeOfFreedom*>& degreesOfFreedom)
{
    degreesOfFreedom_ = degreesOfFreedom;
}

int Element::getIndex() const
{
    return index_;
}

int Element::getRank() const
{
    return rank_;
}

bool Element::isBoundary() const
{
    return isBoundary_;
}

bool Element::isActive() const
{
    return isActive_;
}

ReferenceConfiguration Element::getReferenceConfiguration() const
{
    return referenceConfiguration_;
}

const std::vector<DegreeOfFreedom*>& Element::getDegreesOfFreedom() const
{
    return degreesOfFreedom_;
}

const unsigned int Element::getNumberOfDOFs() const
{
    return degreesOfFreedom_.size();
}

void Element::addNeighborElement(Element* el)
{
    neighborElements_.push_back(el);
}

const std::vector<Element*>& Element::getNeighborElements() const
{
    return neighborElements_;
}

