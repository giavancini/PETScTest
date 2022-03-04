#pragma once
#include "Node.h"

class DirichletBoundaryCondition
{
    private:
        Node* node_;
        DegreeOfFreedom* dof_;
        double value_;

    public:

        DirichletBoundaryCondition(Node* const node, DegreeOfFreedom* const dof, const double value);

        ~DirichletBoundaryCondition();

        Node* getNode() const;

        DegreeOfFreedom* getDegreeOfFreedom() const;

        double getValue() const;
};