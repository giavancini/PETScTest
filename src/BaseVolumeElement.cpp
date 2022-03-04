#include "BaseVolumeElement.h"

BaseVolumeElement::BaseVolumeElement(const int index,
                                     ParametricVolumeElement& parametricElement,     
                                     const std::vector<Node*>& nodes)
        : BaseElement(index, parametricElement, nodes) {}

BaseVolumeElement::~BaseVolumeElement() {}

double BaseVolumeElement::getRadius() const
{
    // This function works only for tetrahedrons, as it is not guaranteed that a hexaedro has a circumcentre

    // Use coordinates relative to point a of the tetrahedron
    double x_ba = nodes_[1]->getDegreeOfFreedom(0)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(0)->getCurrentValue();
    double y_ba = nodes_[1]->getDegreeOfFreedom(1)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(1)->getCurrentValue();
    double z_ba = nodes_[1]->getDegreeOfFreedom(2)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(2)->getCurrentValue();
    double x_ca = nodes_[2]->getDegreeOfFreedom(0)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(0)->getCurrentValue();
    double y_ca = nodes_[2]->getDegreeOfFreedom(1)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(1)->getCurrentValue();
    double z_ca = nodes_[2]->getDegreeOfFreedom(2)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(2)->getCurrentValue();
    double x_da = nodes_[3]->getDegreeOfFreedom(0)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(0)->getCurrentValue();
    double y_da = nodes_[3]->getDegreeOfFreedom(1)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(1)->getCurrentValue();
    double z_da = nodes_[3]->getDegreeOfFreedom(2)->getCurrentValue() - nodes_[0]->getDegreeOfFreedom(2)->getCurrentValue();

    // Squares of lengths of the edges incident to point a
    double length_ba = x_ba * x_ba + y_ba * y_ba + z_ba * z_ba;
    double length_ca = x_ca * x_ca + y_ca * y_ca + z_ca * z_ca;
    double length_da = x_da * x_da + y_da * y_da + z_da * z_da;

    // Cross products of these edges
    double x_cross_cd = y_ca * z_da - y_da * z_ca;
    double y_cross_cd = z_ca * x_da - z_da * x_ca;
    double z_cross_cd = x_ca * y_da - x_da * y_ca;
    double x_cross_db = y_da * z_ba - y_ba * z_da;
    double y_cross_db = z_da * x_ba - z_ba * x_da;
    double z_cross_db = x_da * y_ba - x_ba * y_da;
    double x_cross_bc = y_ba * z_ca - y_ca * z_ba;
    double y_cross_bc = z_ba * x_ca - z_ca * x_ba;
    double z_cross_bc = x_ba * y_ca - x_ca * y_ba;

    // Calculate the denominator of the formula
    double denominator = 0.5 / (x_ba * x_cross_cd + y_ba * y_cross_cd + z_ba * z_cross_cd);

    // Calculate circumcentre relative to point a
    double x_circumcentre = (length_ba * x_cross_cd + length_ca * x_cross_db + length_da * x_cross_bc) *
            denominator;
    double y_circumcentre = (length_ba * y_cross_cd + length_ca * y_cross_db + length_da * y_cross_bc) *
            denominator;
    double z_circumcentre = (length_ba * z_cross_cd + length_ca * z_cross_db + length_da * z_cross_bc) *
            denominator;
    
    // Circumradius is the norm of the circumcentre vector
    return sqrt((x_circumcentre * x_circumcentre) + (y_circumcentre * y_circumcentre) + (z_circumcentre * z_circumcentre));
}

double BaseVolumeElement::getJacobianIntegration() const
{
        const std::vector<QuadraturePoint*>& quadraturePoints = parametricElement_->getQuadraturePoints();

        double volume = 0.0;
        for (auto& qp : quadraturePoints)
        {
            double* phi = qp->getShapeFunctionsValues();
            double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
            double weight = qp->getWeight();

            double dy_dxsi[3][3];
            getCurrentJacobianMatrix(dphi_dxsi, dy_dxsi);
            double j = getMatrixDeterminant(dy_dxsi);

            volume += j * weight;
        }
        return volume;
}

inline void BaseVolumeElement::getCurrentJacobianMatrix(double** dphi_dxsi,
                                                        double dy_dxsi[3][3]) const
{
    dy_dxsi[0][0] = 0.0; dy_dxsi[0][1] = 0.0; dy_dxsi[0][2] = 0.0;
    dy_dxsi[1][0] = 0.0; dy_dxsi[1][1] = 0.0; dy_dxsi[1][2] = 0.0;
    dy_dxsi[2][0] = 0.0; dy_dxsi[2][1] = 0.0; dy_dxsi[2][2] = 0.0;
    
    const unsigned int numberOfNodes = nodes_.size();
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        double y[3];
        y[0] = nodes_[i]->getDegreeOfFreedom(0)->getCurrentValue();
        y[1] = nodes_[i]->getDegreeOfFreedom(1)->getCurrentValue();
        y[2] = nodes_[i]->getDegreeOfFreedom(2)->getCurrentValue();
        dy_dxsi[0][0] += dphi_dxsi[0][i] * y[0];
        dy_dxsi[0][1] += dphi_dxsi[1][i] * y[0];
        dy_dxsi[0][2] += dphi_dxsi[2][i] * y[0];
        dy_dxsi[1][0] += dphi_dxsi[0][i] * y[1];
        dy_dxsi[1][1] += dphi_dxsi[1][i] * y[1];
        dy_dxsi[1][2] += dphi_dxsi[2][i] * y[1];
        dy_dxsi[2][0] += dphi_dxsi[0][i] * y[2];
        dy_dxsi[2][1] += dphi_dxsi[1][i] * y[2];
        dy_dxsi[2][2] += dphi_dxsi[2][i] * y[2];
    }
}

inline double BaseVolumeElement::getMatrixDeterminant(const double jacobianMatrix[3][3]) const
{
    return jacobianMatrix[0][0] * jacobianMatrix[1][1] * jacobianMatrix[2][2] +
           jacobianMatrix[0][1] * jacobianMatrix[1][2] * jacobianMatrix[2][0] +
           jacobianMatrix[0][2] * jacobianMatrix[1][0] * jacobianMatrix[2][1] -
           jacobianMatrix[0][2] * jacobianMatrix[1][1] * jacobianMatrix[2][0] -
           jacobianMatrix[0][1] * jacobianMatrix[1][0] * jacobianMatrix[2][2] -
           jacobianMatrix[0][0] * jacobianMatrix[1][2] * jacobianMatrix[2][1];
        
}
