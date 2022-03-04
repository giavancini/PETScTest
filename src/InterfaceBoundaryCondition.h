#pragma once
#include "Node.h"
#include "BaseLineElement.h"
#include "BaseSurfaceElement.h"

class InterfaceNeumannBoundaryCondition
{
public:
    InterfaceNeumannBoundaryCondition(const int index,
                             const int ndofs);

    virtual ~InterfaceNeumannBoundaryCondition() = default;

    int getNumberOfDOFs() const;

    virtual void getNodalForce(const int& dimension,
                               std::vector<DegreeOfFreedom*>& dofs,
                               double*& values) const = 0;

protected:
    int index_;
    int ndofs_;
};

class InterfaceLineLoad : public InterfaceNeumannBoundaryCondition
{
public:
    InterfaceLineLoad(const int index,
                      const int ndofs,
                      BaseLineElement* const element);

    ~InterfaceLineLoad() override;

    void getNodalForce(const int& dimension,
                       std::vector<DegreeOfFreedom*>& dofs,
                       double*& values) const override;

private:
    BaseLineElement* element_;
};

class InterfaceSurfaceLoad : public InterfaceNeumannBoundaryCondition
{
public:
    InterfaceSurfaceLoad(const int index,
                      const int ndofs,
                      BaseSurfaceElement* const element);

    ~InterfaceSurfaceLoad() override;

    void getNodalForce(const int& dimension,
                       std::vector<DegreeOfFreedom*>& dofs,
                       double*& values) const override;

private:
    BaseSurfaceElement* element_;
};