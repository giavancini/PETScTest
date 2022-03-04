#include "DegreeOfFreedom.h"

DegreeOfFreedom::DegreeOfFreedom(const DOFType& type, const double value)
    : index_(0),
      type_(type),
      initialValue_(value),
      currentValue_(value),
      pastValue_(value),
      intermediateValue_(value),
      initialFirstTimeDerivative_(0.0),
      currentFirstTimeDerivative_(0.0),
      pastFirstTimeDerivative_(0.0),
      intermediateFirstTimeDerivative_(0.0),
      initialSecondTimeDerivative_(0.0),
      currentSecondTimeDerivative_(0.0),
      pastSecondTimeDerivative_(0.0),
      intermediateSecondTimeDerivative_(0.0) {}

int DegreeOfFreedom::getIndex() const
{ 
    return index_;
}

DOFType DegreeOfFreedom::getType() const
{
    return type_;
}

double DegreeOfFreedom::getInitialValue() const
{
    return initialValue_;
}

double DegreeOfFreedom::getCurrentValue() const
{
    return currentValue_;
}

double DegreeOfFreedom::getPastValue() const
{
    return pastValue_;
}

double DegreeOfFreedom::getIntermediateValue() const
{
    return intermediateValue_;
}

double DegreeOfFreedom::getInitialFirstTimeDerivative() const
{
    return initialFirstTimeDerivative_;
}

double DegreeOfFreedom::getCurrentFirstTimeDerivative() const
{ 
    return currentFirstTimeDerivative_;
}

double DegreeOfFreedom::getPastFirstTimeDerivative() const
{
    return pastFirstTimeDerivative_;
}

double DegreeOfFreedom::getIntermediateFirstTimeDerivative() const
{
    return intermediateFirstTimeDerivative_;
}

double DegreeOfFreedom::getInitialSecondTimeDerivative() const
{ 
    return initialSecondTimeDerivative_;
}

double DegreeOfFreedom::getCurrentSecondTimeDerivative() const
{
    return currentSecondTimeDerivative_;
}

double DegreeOfFreedom::getPastSecondTimeDerivative() const
{
    return pastSecondTimeDerivative_;
}

double DegreeOfFreedom::getIntermediateSecondTimeDerivative() const
{
    return intermediateSecondTimeDerivative_;
}

void DegreeOfFreedom::setIndex(const int& index)
{
    index_ = index;
}

void DegreeOfFreedom::setInitialValue(const double& value)
{
    initialValue_ = value;
}

void DegreeOfFreedom::setCurrentValue(const double& value)
{
    currentValue_ = value;
}

void DegreeOfFreedom::setPastValue(const double& value)
{
    pastValue_ = value;
}

void DegreeOfFreedom::setIntermediateValue(const double& value)
{
    intermediateValue_ = value;
}

void DegreeOfFreedom::incrementCurrentValue(const double& value)
{
    currentValue_ += value;
}

void DegreeOfFreedom::setInitialFirstTimeDerivative(const double& value)
{
    initialFirstTimeDerivative_ = value;
}

void DegreeOfFreedom::setCurrentFirstTimeDerivative(const double& value)
{
    currentFirstTimeDerivative_ = value;
}

void DegreeOfFreedom::setPastFirstTimeDerivative(const double& value)
{
    pastFirstTimeDerivative_ = value;
}

void DegreeOfFreedom::setIntermediateFirstTimeDerivative(const double& value)
{
    intermediateFirstTimeDerivative_ = value;
}

void DegreeOfFreedom::setInitialSecondTimeDerivative(const double& value)
{
    initialSecondTimeDerivative_ = value;
}

void DegreeOfFreedom::setCurrentSecondTimeDerivative(const double& value)
{
    currentSecondTimeDerivative_ = value;
}

void DegreeOfFreedom::setPastSecondTimeDerivative(const double& value)
{
    pastSecondTimeDerivative_ = value;
}

void DegreeOfFreedom::setIntermediateSecondTimeDerivative(const double& value)
{
    intermediateSecondTimeDerivative_ = value;
}