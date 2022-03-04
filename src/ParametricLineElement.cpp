#include "ParametricLineElement.h"

ParametricLineElement ParametricLineElement::L2(PartitionOfUnity::L2);
ParametricLineElement ParametricLineElement::L3(PartitionOfUnity::L3);
ParametricLineElement ParametricLineElement::L4(PartitionOfUnity::L4);

ParametricLineElement::ParametricLineElement(const PartitionOfUnity elementType)
    :ParametricElement(elementType)
{
    vtkCellType_ = VTK_LAGRANGE_CURVE;

    switch (elementType)
    {
    case PartitionOfUnity::L2:
        order_ = 1;
        numberOfQuadraturePoints_ = 1;
        numberOfNodes_ = 2;
        numberOfFaces_ = 2;
        vtkConnectivity_ = {0, 1};
        faceNodes_ = {{0}, {1}};
        faceVertices_ = {{0}, {1}};
        break;
    case PartitionOfUnity::L3:
        order_ = 2;
        numberOfQuadraturePoints_ = 2;
        numberOfNodes_ = 3;
        numberOfFaces_ = 2;
        vtkConnectivity_ = {0, 1, 2};
        faceNodes_ = {{0}, {1}};
        faceVertices_ = {{0}, {1}};
        break;
    case PartitionOfUnity::L4:
        order_ = 3;
        numberOfQuadraturePoints_ = 3;
        numberOfNodes_ = 4;
        numberOfFaces_ = 2;
        vtkConnectivity_ = {0, 1, 2, 3};
        faceNodes_ = {{0}, {1}};
        faceVertices_ = {{0}, {1}};
        break;
    default:
        std::cout << "Parametric line element not implemented!\n";
        exit(EXIT_FAILURE);
        break;
    }
    double **xsi, *weight;
    quadratures::lineQuadrature(numberOfQuadraturePoints_, xsi, weight);
    quadraturePoints_.reserve(numberOfQuadraturePoints_);
    for (int i = 0; i < numberOfQuadraturePoints_; i++)
    {
        double *phi, **dphi_dxsi;
        getShapeFunctions(xsi[i], phi);
        getShapeFunctionsDerivatives(xsi[i], dphi_dxsi);
        
        quadraturePoints_.emplace_back(new QuadraturePoint(xsi[i], phi, dphi_dxsi, weight[i], 1, numberOfNodes_));
    }
    setNodalParametricCoordinates();
}

ParametricLineElement::~ParametricLineElement()
{
    for (QuadraturePoint* qp : quadraturePoints_)
    {
        delete qp;
    }
    delete[] nodalParametricCoordinates_[0];
    delete[] nodalParametricCoordinates_;
}

void ParametricLineElement::setNodalParametricCoordinates()
{
    nodalParametricCoordinates_ = new double*[1];
    nodalParametricCoordinates_[0] = new double[numberOfNodes_];

    switch (elementType_)
    {
    case PartitionOfUnity::L2:
        nodalParametricCoordinates_[0][0] = -1.0;
        nodalParametricCoordinates_[0][1] =  1.0;
        break;
    case PartitionOfUnity::L3:
        nodalParametricCoordinates_[0][0] = -1.0;
        nodalParametricCoordinates_[0][1] =  1.0;
        nodalParametricCoordinates_[0][2] =  0.0;
        break;
    case PartitionOfUnity::L4:
        nodalParametricCoordinates_[0][0] = -1.0;
        nodalParametricCoordinates_[0][1] =  1.0;
        nodalParametricCoordinates_[0][2] = -1.0 / 3.0;
        nodalParametricCoordinates_[0][3] =  1.0 / 3.0;
        break;
    default:
        break;
    }
}

void ParametricLineElement::getShapeFunctions(double* xsi, double*& phi) const
{
    phi = new double[order_ + 1];
    double xsi1 = xsi[0];

    switch (elementType_)
    {
    case PartitionOfUnity::L2:
        phi[0] = 0.5 * (1.0 - xsi1);
	    phi[1] = 0.5 * (1.0 + xsi1);
        break;
    case PartitionOfUnity::L3:
        phi[0] = xsi1 * (xsi1 - 1.0) * 0.5;
        phi[1] = xsi1 * (xsi1 + 1.0) * 0.5;
        phi[2] = (1.0 + xsi1) * (1.0 - xsi1);
        break;
    case PartitionOfUnity::L4:
        phi[0] = -9.0 / 16.0 * (xsi1 + 1.0 / 3.0) * (xsi1 - 1.0 / 3.0) * (xsi1 - 1.0);
        phi[1] = 9.0 / 16.0 * (xsi1 + 1.0) * (xsi1 + 1.0 / 3.0) * (xsi1 - 1.0 / 3.0);
        phi[2] = 27.0 / 16.0 * (xsi1 + 1.0) * (xsi1 - 1.0 / 3.0) * (xsi1 - 1.0);
        phi[3] = -27.0 / 16.0 * (xsi1 + 1.0) * (xsi1 + 1.0 / 3.0) * (xsi1 - 1.0);
        break;
    default:
        break;
    }
}

void ParametricLineElement::getShapeFunctionsDerivatives(double* xsi, double**& dphi_dxsi) const
{
    dphi_dxsi = new double*;
    dphi_dxsi[0] = new double[order_ + 1];
    double xsi1 = xsi[0];

    switch (elementType_)
    {
    case PartitionOfUnity::L2:
        dphi_dxsi[0][0] = -0.5;
        dphi_dxsi[0][1] = 0.5;
        break;
    case PartitionOfUnity::L3:
        dphi_dxsi[0][0] = xsi1 - 0.5;
        dphi_dxsi[0][1] = xsi1 + 0.5;
        dphi_dxsi[0][2] = -2.0 * xsi1;
        break;
    case PartitionOfUnity::L4:
        dphi_dxsi[0][0] = -9.0/16.0  * ((xsi1-1.0/3.0) * (xsi1-1.0) + (xsi1+1.0/3.0) * (xsi1-1.0) + (xsi1+1.0/3.0) * (xsi1-1.0/3.0));
        dphi_dxsi[0][1] =  9.0/16.0  * ((xsi1+1.0/3.0) * (xsi1-1.0/3.0) + (xsi1+1.0) * (xsi1-1.0/3.0) + (xsi1+1.0) * (xsi1+1.0/3.0));
        dphi_dxsi[0][2] = 27.0/16.0  * ((xsi1-1.0/3.0) * (xsi1-1.0) + (xsi1+1.0) * (xsi1-1.0) + (xsi1+1.0) * (xsi1-1.0/3.0));
        dphi_dxsi[0][3] = -27.0/16.0 * ((xsi1+1.0/3.0) * (xsi1-1.0) + (xsi1+1.0) * (xsi1-1.0) + (xsi1+1.0) * (xsi1+1.0/3.0));
        break;
    default:
        break;
    }
}