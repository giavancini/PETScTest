#pragma once
#include "Node.h"
#include "BaseLineElement.h"
#include "BaseSurfaceElement.h"

class NeumannBoundaryCondition
{
    protected:
        int index_;
        int ndofs_;
        double forces_[3];

    public:

        NeumannBoundaryCondition(const int index, const int ndofs, const double valueX, const double valueY, const double valueZ);

        virtual ~NeumannBoundaryCondition() = 0;

        int getNumberOfDOFs() const;

        double getForce(const int& dof) const;

        virtual void getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const = 0;
};

class PointLoad : public NeumannBoundaryCondition
{
    private:
        Node* node_;
    
    public:

        PointLoad(const int index, const int ndofs, Node* const node, const double valueX, const double valueY, const double valueZ);

        ~PointLoad() override;

        void getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const override;
};

class LineLoad : public NeumannBoundaryCondition
{
    private:
        BaseLineElement* element_;

    public:

        LineLoad(const int index, const int ndofs, BaseLineElement* const element, const double valueX, const double valueY, const double valueZ);

        ~LineLoad() override;

        void getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const override;
};

class SurfaceLoad : public NeumannBoundaryCondition
{
    private:
        BaseSurfaceElement* element_;

    public:

        SurfaceLoad(const int index, const int ndofs, BaseSurfaceElement* const element, const double valueX, const double valueY, const double valueZ);

        ~SurfaceLoad() override;
        
        void getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const override;
};