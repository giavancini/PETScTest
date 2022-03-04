#include "InterfaceBoundaryCondition.h"
#include <cmath>
#include <iostream>

InterfaceNeumannBoundaryCondition::InterfaceNeumannBoundaryCondition(const int index,
                                                                     const int ndofs)
    : index_(index),
      ndofs_(ndofs) {}

int InterfaceNeumannBoundaryCondition::getNumberOfDOFs() const
{
    return ndofs_;
}

InterfaceLineLoad::InterfaceLineLoad(const int index,
                                     const int ndofs,
                                     BaseLineElement* const element)
    : InterfaceNeumannBoundaryCondition(index, ndofs),
      element_(element) {}

InterfaceLineLoad::~InterfaceLineLoad() {}

void InterfaceLineLoad::getNodalForce(const int& dimension,
                                      std::vector<DegreeOfFreedom*>& dofs,
                                      double*& values) const
{
    dofs.reserve(ndofs_);
    values = new double[ndofs_];

    for (int i = 0; i < ndofs_; i++)
        values[i] = 0.0;
    
    const int numberOfNodes = element_->getNumberOfNodes();
    const std::vector<Node*>& nodes = element_->getNodes();
    
    const std::vector<QuadraturePoint*>& quadraturePoints = element_->getParametricElement()->getQuadraturePoints();

    for (auto& qp : quadraturePoints)
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dx_dxsi[2];
        double interpolatedForce[2];
        double interpolatedStress[3];
        for (int i = 0; i < 2; i++)
        {
            dx_dxsi[i] = 0.0;
        }
        for (int i = 0; i < 3; i++)
        {
            interpolatedStress[i] = 0.0;
        }

        for (int i = 0; i < numberOfNodes; i++)
        {
            Node* fluidNode = nodes[i]->getInterfaceNode();
            double* stress = fluidNode->getCauchyStress();
            for (int j = 0; j < 2; j++)
            {
                double coord = nodes[i]->getDegreeOfFreedom(j)->getCurrentValue();
                dx_dxsi[j] += dphi_dxsi[0][i] * coord;
                stress[i] += fluidNode->getDegreeOfFreedom(2)->getCurrentValue(); //Add pressure to the diagonal of stress tensor
            }
            for (int j = 0; j < 3; j++)
            {
                interpolatedStress[j] += stress[j] * phi[i];
            }
        }

        // double jac = 0.0;
        // for (int i = 0; i < 2; i++)
        // {
        //     jac += dx_dxsi[i] * dx_dxsi[i];
        // }
        // jac = sqrt(jac);
        
        double normal[2];
        normal[0] = -dx_dxsi[1];
        normal[1] = dx_dxsi[0];

        interpolatedForce[0] = interpolatedStress[0] * normal[0] + interpolatedStress[2] * normal[1];
        interpolatedForce[1] = interpolatedStress[2] * normal[0] + interpolatedStress[1] * normal[1];

        for (int i = 0; i < numberOfNodes; i++)
        {
            double factor = phi[i] * weight;
            for (int j = 0; j < 2; j++)
            {
                values[dimension*i+j] += interpolatedForce[j] * factor;
            }
        }
    }

    for (Node* const& node : nodes)
    {
        for (int i = 0; i < 2; i++)
        {
            dofs.push_back(node->getDegreeOfFreedom(i));
        }
    }
}

InterfaceSurfaceLoad::InterfaceSurfaceLoad(const int index,
                                           const int ndofs,
                                           BaseSurfaceElement* const element)
    : InterfaceNeumannBoundaryCondition(index, ndofs),
      element_(element) {}

InterfaceSurfaceLoad::~InterfaceSurfaceLoad() {}

void InterfaceSurfaceLoad::getNodalForce(const int& dimension,
                                         std::vector<DegreeOfFreedom*>& dofs,
                                         double*& values) const
{
    dofs.reserve(ndofs_);
    values = new double[ndofs_];

    for (int i = 0; i < ndofs_; i++)
        values[i] = 0.0;
    
    const int numberOfNodes = element_->getNumberOfNodes();
    const std::vector<Node*>& nodes = element_->getNodes();
    
    const std::vector<QuadraturePoint*>& quadraturePoints = element_->getParametricElement()->getQuadraturePoints();

    for (auto& qp : quadraturePoints)
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dx_dxsi1[3], dx_dxsi2[3];
        dx_dxsi1[0] = 0.0; dx_dxsi1[1] = 0.0; dx_dxsi1[2] = 0.0;
        dx_dxsi2[0] = 0.0; dx_dxsi2[1] = 0.0; dx_dxsi2[2] = 0.0;
        double interpolatedForce[3];
        double interpolatedStress[6];
        interpolatedStress[0] = 0.0; interpolatedStress[1] = 0.0; interpolatedStress[2] = 0.0;
        interpolatedStress[3] = 0.0; interpolatedStress[4] = 0.0; interpolatedStress[5] = 0.0;

        for (int i = 0; i < numberOfNodes; i++)
        {
            double coords[3];
            coords[0] = nodes[i]->getDegreeOfFreedom(0)->getCurrentValue();
            coords[1] = nodes[i]->getDegreeOfFreedom(1)->getCurrentValue();
            coords[2] = nodes[i]->getDegreeOfFreedom(2)->getCurrentValue();

            dx_dxsi1[0] += dphi_dxsi[0][i] * coords[0];
            dx_dxsi1[1] += dphi_dxsi[0][i] * coords[1];
            dx_dxsi1[2] += dphi_dxsi[0][i] * coords[2];
            
            dx_dxsi2[0] += dphi_dxsi[1][i] * coords[0];
            dx_dxsi2[1] += dphi_dxsi[1][i] * coords[1];
            dx_dxsi2[2] += dphi_dxsi[1][i] * coords[2];

            Node* fluidNode = nodes[i]->getInterfaceNode();
            double* stress = fluidNode->getCauchyStress();
            stress[0] += fluidNode->getDegreeOfFreedom(3)->getCurrentValue();
            stress[1] += fluidNode->getDegreeOfFreedom(3)->getCurrentValue();
            stress[2] += fluidNode->getDegreeOfFreedom(3)->getCurrentValue();
            for (int j = 0; j < 6; j++)
            {
                interpolatedStress[j] += stress[j] * phi[i];
            }
        }

        double normal[3];
        normal[0] = dx_dxsi1[1] * dx_dxsi2[2] - dx_dxsi1[2] * dx_dxsi2[1];
        normal[1] = dx_dxsi1[2] * dx_dxsi2[0] - dx_dxsi1[0] * dx_dxsi2[2];
        normal[2] = dx_dxsi1[0] * dx_dxsi2[1] - dx_dxsi1[1] * dx_dxsi2[0];

        interpolatedForce[0] = interpolatedStress[0] * normal[0] + interpolatedStress[3] * normal[1] + interpolatedStress[4] * normal[2];
        interpolatedForce[1] = interpolatedStress[3] * normal[0] + interpolatedStress[1] * normal[1] + interpolatedStress[5] * normal[2];
        interpolatedForce[2] = interpolatedStress[4] * normal[0] + interpolatedStress[5] * normal[1] + interpolatedStress[2] * normal[2];

        for (int i = 0; i < numberOfNodes; i++)
        {
            double factor = phi[i] * weight;
            values[3*i+0] += interpolatedForce[0] * factor;
            values[3*i+1] += interpolatedForce[1] * factor;
            values[3*i+2] += interpolatedForce[2] * factor;
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
