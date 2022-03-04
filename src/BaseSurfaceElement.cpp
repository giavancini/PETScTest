#include "BaseSurfaceElement.h"
#include <algorithm>

BaseSurfaceElement::BaseSurfaceElement(const int index,
                                       ParametricSurfaceElement& parametricElement,     
                                       const std::vector<Node*>& nodes)
        : BaseElement(index, parametricElement, nodes) {}

BaseSurfaceElement::~BaseSurfaceElement() {}

double BaseSurfaceElement::getRadius() const
{
    // This function works only for triangles, as it is not guaranteed that a quadrilateral has a circumcentre
    
    // Use coordinates relative to point a of the triangle
    double x_ba = nodes_[1]->getDegreeOfFreedom(0)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(0)->getCurrentValue();
    double y_ba = nodes_[1]->getDegreeOfFreedom(1)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(1)->getCurrentValue();
    double x_ca = nodes_[2]->getDegreeOfFreedom(0)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(0)->getCurrentValue();
    double y_ca = nodes_[2]->getDegreeOfFreedom(1)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(1)->getCurrentValue();

    // Squares of lengths of the edges incident to point a
    double length_ba = x_ba * x_ba + y_ba * y_ba;
    double length_ca = x_ca * x_ca + y_ca * y_ca;

    // Calculate the denominator of the formula
    double denominator = 0.5 / (x_ba * y_ca - y_ba * x_ca);
    
    // Calculate circumcentre relative to point a
    double x_circumcentre = (y_ca * length_ba - y_ba * length_ca) * denominator;
    double y_circumcentre = (x_ba * length_ca - x_ca * length_ba) * denominator;
    
    // Circumradius is the norm of the circumcentre vector
    return sqrt((x_circumcentre * x_circumcentre) + (y_circumcentre * y_circumcentre));
}

double BaseSurfaceElement::getJacobianIntegration() const
{
    const std::vector<QuadraturePoint*>& quadraturePoints = parametricElement_->getQuadraturePoints();

    double area = 0.0;
    for (auto& qp : quadraturePoints)
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dy_dxsi[2][2];
        getCurrentJacobianMatrix(dphi_dxsi, dy_dxsi);
        double j = getMatrixDeterminant(dy_dxsi);

        area += j * weight;
    }
    return area;
}

void BaseSurfaceElement::getInitialNormalVector(double* xsi,
                                                double normal[3]) const
{
    double** dphi_dxsi;
    parametricElement_->getShapeFunctionsDerivatives(xsi, dphi_dxsi);
    const int nnodes = parametricElement_->getNumberOfNodes();

    double dx_dxsi[2][3];
    dx_dxsi[0][0] = 0.0; dx_dxsi[0][1] = 0.0; dx_dxsi[0][2] = 0.0;
    dx_dxsi[1][0] = 0.0; dx_dxsi[1][1] = 0.0; dx_dxsi[1][2] = 0.0;
    
    for (int i = 0; i < nnodes; i++)
    {
        double x[3];
        x[0] = nodes_[i]->getDegreeOfFreedom(0)->getInitialValue();
        x[1] = nodes_[i]->getDegreeOfFreedom(1)->getInitialValue();
        x[2] = nodes_[i]->getDegreeOfFreedom(2)->getInitialValue();

        dx_dxsi[0][0] += dphi_dxsi[0][i] * x[0];
        dx_dxsi[0][1] += dphi_dxsi[0][i] * x[1];
        dx_dxsi[0][2] += dphi_dxsi[0][i] * x[2];
        dx_dxsi[1][0] += dphi_dxsi[1][i] * x[0];
        dx_dxsi[1][1] += dphi_dxsi[1][i] * x[1];
        dx_dxsi[1][2] += dphi_dxsi[1][i] * x[2];
    }

    normal[0] = dx_dxsi[0][1] * dx_dxsi[1][2] - dx_dxsi[0][2] * dx_dxsi[1][1];
    normal[1] = dx_dxsi[0][2] * dx_dxsi[1][0] - dx_dxsi[0][0] * dx_dxsi[1][2];
    normal[2] = dx_dxsi[0][0] * dx_dxsi[1][1] - dx_dxsi[0][1] * dx_dxsi[1][0];

    double norm = std::sqrt(pow(normal[0], 2) + pow(normal[1], 2) + pow(normal[2], 2));
    norm = 1.0 / norm;

    normal[0] *= norm;
    normal[1] *= norm;
    normal[2] *= norm;
}

inline void BaseSurfaceElement::getCurrentJacobianMatrix(double** dphi_dxsi,
                                                         double dy_dxsi[2][2]) const
{
    dy_dxsi[0][0] = 0.0; dy_dxsi[0][1] = 0.0;
    dy_dxsi[1][0] = 0.0; dy_dxsi[1][1] = 0.0;
    
    const unsigned int numberOfNodes = nodes_.size();
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        double y[2];
        y[0] = nodes_[i]->getDegreeOfFreedom(0)->getCurrentValue();
        y[1] = nodes_[i]->getDegreeOfFreedom(1)->getCurrentValue();
        dy_dxsi[0][0] += dphi_dxsi[0][i] * y[0];
        dy_dxsi[0][1] += dphi_dxsi[1][i] * y[0];
        dy_dxsi[1][0] += dphi_dxsi[0][i] * y[1];
        dy_dxsi[1][1] += dphi_dxsi[1][i] * y[1];
    }
}

inline double BaseSurfaceElement::getMatrixDeterminant(const double jacobianMatrix[2][2]) const
{
    return jacobianMatrix[0][0] * jacobianMatrix[1][1] -
           jacobianMatrix[0][1] * jacobianMatrix[1][0];
        
}

