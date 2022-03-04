#include "NeumannBoundaryCondition.h"
#include <cmath>
#include <iostream>

NeumannBoundaryCondition::NeumannBoundaryCondition(const int index, const int ndofs, const double valueX, const double valueY, const double valueZ)
    : index_(index), ndofs_(ndofs), forces_{valueX, valueY, valueZ}  {}

NeumannBoundaryCondition::~NeumannBoundaryCondition() {}

int NeumannBoundaryCondition::getNumberOfDOFs() const
{
    return ndofs_;
}

double NeumannBoundaryCondition::getForce(const int& dof) const
{
    return forces_[dof];
}

PointLoad::PointLoad(const int index, const int ndofs, Node* const node, const double valueX, const double valueY, const double valueZ)
    : NeumannBoundaryCondition(index, ndofs, valueX, valueY, valueZ), node_(node) {}

PointLoad::~PointLoad() {}

void PointLoad::getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const
{
    dofs.reserve(ndofs_);
    values = new double[ndofs_];
    for (int i = 0; i < dimension; i++)
    {
        dofs.push_back(node_->getDegreeOfFreedom(i));
        values[i] = forces_[i];
    }
}

LineLoad::LineLoad(const int index, const int ndofs, BaseLineElement* const element, const double valueX, const double valueY, const double valueZ)
    : NeumannBoundaryCondition(index, ndofs, valueX, valueY, valueZ), element_(element) {}

LineLoad::~LineLoad() {}

void LineLoad::getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const
{
    dofs.reserve(ndofs_);
    values = new double[ndofs_];

    for (int i = 0; i < ndofs_; i++)
        values[i] = 0.0;
    
    const int numberOfNodes = element_->getNumberOfNodes();
    const std::vector<Node*>& nodes = element_->getNodes();
    
    std::vector<QuadraturePoint*> quadraturePoints = element_->getParametricElement()->getQuadraturePoints();

    for (auto& qp : quadraturePoints)
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dx_dxsi[dimension];
        double interpolatedForce[dimension];
        for (int i = 0; i < dimension; i++)
        {
            dx_dxsi[i] = 0.0;
            interpolatedForce[i] = 0.0;
        }
        
        for (int i = 0; i < numberOfNodes; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                double coord = nodes[i]->getDegreeOfFreedom(j)->getInitialValue();
                dx_dxsi[j] += dphi_dxsi[0][i] * coord;
                interpolatedForce[j] += phi[i] * forces_[j];
            }
        }
        double j0 = 0.0;
        for (int i = 0; i < dimension; i++)
        {
            j0 += dx_dxsi[i] * dx_dxsi[i];
        }
        j0 = sqrt(j0);

        for (int i = 0; i < numberOfNodes; i++)
        {
            for (int j = 0; j < dimension; j++)
            {
                values[dimension*i+j] += interpolatedForce[j] * phi[i] * weight * j0;
            }
        }
    }

    for (Node* const& node : nodes)
    {
        for (int i = 0; i < dimension; i++)
        {
            dofs.push_back(node->getDegreeOfFreedom(i));
        }
    }
}

SurfaceLoad::SurfaceLoad(const int index, const int ndofs, BaseSurfaceElement* const element, const double valueX, const double valueY, const double valueZ)
    : NeumannBoundaryCondition(index, ndofs, valueX, valueY, valueZ), element_(element) {}

SurfaceLoad::~SurfaceLoad() {}

void SurfaceLoad::getNodalForce(const int& dimension, std::vector<DegreeOfFreedom*>& dofs, double*& values) const
{
    dofs.reserve(ndofs_);
    values = new double[ndofs_];

    for (int i = 0; i < ndofs_; i++)
        values[i] = 0.0;
    
    const int numberOfNodes = element_->getNumberOfNodes();
    const std::vector<Node*>& nodes = element_->getNodes();
    
    std::vector<QuadraturePoint*> quadraturePoints = element_->getParametricElement()->getQuadraturePoints();

    for (auto& qp : quadraturePoints)
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dx_dxsi1[3], dx_dxsi2[3];
        dx_dxsi1[0] = 0.0; dx_dxsi1[1] = 0.0; dx_dxsi1[2] = 0.0;
        dx_dxsi2[0] = 0.0; dx_dxsi2[1] = 0.0; dx_dxsi2[2] = 0.0;
        double interpolatedForce[3];
        interpolatedForce[0] = 0.0; interpolatedForce[1] = 0.0; interpolatedForce[2] = 0.0;

        for (int i = 0; i < numberOfNodes; i++)
        {
            double coords[3];
            coords[0] = nodes[i]->getDegreeOfFreedom(0)->getInitialValue();
            coords[1] = nodes[i]->getDegreeOfFreedom(1)->getInitialValue();
            coords[2] = nodes[i]->getDegreeOfFreedom(2)->getInitialValue();

            dx_dxsi1[0] += dphi_dxsi[0][i] * coords[0];
            dx_dxsi1[1] += dphi_dxsi[0][i] * coords[1];
            dx_dxsi1[2] += dphi_dxsi[0][i] * coords[2];
            
            dx_dxsi2[0] += dphi_dxsi[1][i] * coords[0];
            dx_dxsi2[1] += dphi_dxsi[1][i] * coords[1];
            dx_dxsi2[2] += dphi_dxsi[1][i] * coords[2];

            interpolatedForce[0] += phi[i] * forces_[0];
            interpolatedForce[1] += phi[i] * forces_[1];
            interpolatedForce[2] += phi[i] * forces_[2];
        }
        double j0 = (dx_dxsi1[0]*dx_dxsi1[0] + dx_dxsi1[1]*dx_dxsi1[1] + dx_dxsi1[2]*dx_dxsi1[2]) * 
                    (dx_dxsi2[0]*dx_dxsi2[0] + dx_dxsi2[1]*dx_dxsi2[1] + dx_dxsi2[2]*dx_dxsi2[2]) - 
					(dx_dxsi1[0]*dx_dxsi2[0] + dx_dxsi1[1]*dx_dxsi2[1] + dx_dxsi1[2]*dx_dxsi2[2]) * 
                    (dx_dxsi2[0]*dx_dxsi1[0] + dx_dxsi2[1]*dx_dxsi1[1] + dx_dxsi2[2]*dx_dxsi1[2]);
		j0 = sqrt(j0);

        double factor = weight * j0;

        for (int i = 0; i < numberOfNodes; i++)
        {
            values[3*i+0] += interpolatedForce[0] * phi[i] * factor;
            values[3*i+1] += interpolatedForce[1] * phi[i] * factor;
            values[3*i+2] += interpolatedForce[2] * phi[i] * factor;
        }
    }

    for (Node* const& node : nodes)
    {
        for (int i = 0; i < dimension; i++)
        {
            dofs.push_back(node->getDegreeOfFreedom(i));
        }
    }
}