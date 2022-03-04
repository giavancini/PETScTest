#include "QuadraturePoint.h"

QuadraturePoint::QuadraturePoint(double* xsi,
                                 double* phi,
                                 double** dphi_dxsi,
                                 const double weight,
                                 const int dimension,
                                 const int numberOfNodes)
                : xsi_(xsi),
                  phi_(phi),
                  dphi_dxsi_(dphi_dxsi),
                  weight_(weight),
                  dimension_(dimension),
                  numberOfNodes_(numberOfNodes) {}

QuadraturePoint::~QuadraturePoint()
{
    delete[] xsi_;
    delete[] phi_;
    for (int i = 0; i < dimension_; i++)
    {
        delete[] dphi_dxsi_[i];
    }
    delete[] dphi_dxsi_;
}

double* QuadraturePoint::getParametricCoordinates() const
{
    return xsi_;
}

double* QuadraturePoint::getShapeFunctionsValues() const
{
    return phi_;
}

double** QuadraturePoint::getShapeFunctionsDerivativesValues() const
{
    return dphi_dxsi_;
}

double QuadraturePoint::getWeight() const
{
    return weight_;
}

int QuadraturePoint::getDimension() const
{
    return dimension_;
}

int QuadraturePoint::getNumberOfNodes() const
{
    return numberOfNodes_;
}