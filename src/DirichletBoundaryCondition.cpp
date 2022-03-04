#include "DirichletBoundaryCondition.h"

DirichletBoundaryCondition::DirichletBoundaryCondition(Node* const node, DegreeOfFreedom* const dof, const double value)
    : node_(node), dof_(dof), value_(value) {}

DirichletBoundaryCondition::~DirichletBoundaryCondition() {}

Node* DirichletBoundaryCondition::getNode() const
{
    return node_;
}

DegreeOfFreedom* DirichletBoundaryCondition::getDegreeOfFreedom() const
{
    return dof_;
}

double DirichletBoundaryCondition::getValue() const
{
    return value_;
}