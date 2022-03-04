#include "ParametricSurfaceElement.h"

ParametricSurfaceElement ParametricSurfaceElement::T3(PartitionOfUnity::T3);
ParametricSurfaceElement ParametricSurfaceElement::T6(PartitionOfUnity::T6);
ParametricSurfaceElement ParametricSurfaceElement::T10(PartitionOfUnity::T10);
ParametricSurfaceElement ParametricSurfaceElement::Q4(PartitionOfUnity::Q4);
ParametricSurfaceElement ParametricSurfaceElement::Q9(PartitionOfUnity::Q9);
ParametricSurfaceElement ParametricSurfaceElement::Q16(PartitionOfUnity::Q16);

ParametricSurfaceElement::ParametricSurfaceElement(const PartitionOfUnity elementType)
    :ParametricElement(elementType)
{
    double **xsi, *weight;
    switch (elementType)
    {
    case PartitionOfUnity::T3:
    {
        order_ = 1;
        numberOfQuadraturePoints_ = 1;
        numberOfNodes_ = 3;
        numberOfFaces_ = 3;
        numberOfEdges_ = 3;
        vtkCellType_ = VTK_LAGRANGE_TRIANGLE;
        vtkConnectivity_ = {0, 1, 2};
        faceNodes_ = {{1, 2}, {2, 0}, {0, 1}};
        faceVertices_ = {{1, 2}, {2, 0}, {0, 1}};
        edgeNodes_ = {{1, 2}, {2, 0}, {0, 1}};
        edgeVertices_ = {{1, 2}, {2, 0}, {0, 1}};
        quadratures::triangleQuadrature(numberOfQuadraturePoints_, xsi, weight);
        break;
    }
    case PartitionOfUnity::T6:
    {
        order_ = 2;
        numberOfQuadraturePoints_ = 6;
        numberOfNodes_ = 6;
        numberOfFaces_ = 3;
        numberOfEdges_ = 3;
        vtkCellType_ = VTK_LAGRANGE_TRIANGLE;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5};
        faceNodes_ = {{1, 4, 2}, {2, 5, 0}, {0, 3, 1}};
        faceVertices_ = {{1, 2}, {2, 0}, {0, 1}};
        edgeNodes_ = {{1, 4, 2}, {2, 5, 0}, {0, 3, 1}};
        edgeVertices_ = {{1, 2}, {2, 0}, {0, 1}};
        quadratures::triangleQuadrature(numberOfQuadraturePoints_, xsi, weight);
        break;
    }
    case PartitionOfUnity::T10:
    {
        order_ = 3;
        numberOfQuadraturePoints_ = 12;
        numberOfNodes_ = 10;
        numberOfFaces_ = 3;
        numberOfEdges_ = 3;
        vtkCellType_ = VTK_LAGRANGE_TRIANGLE;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
        faceNodes_ = {{1, 5, 6, 2}, {2, 7, 8, 0}, {0, 3, 4, 1}};
        faceVertices_ = {{1, 2}, {2, 0}, {0, 1}};
        edgeNodes_ = {{1, 5, 6, 2}, {2, 7, 8, 0}, {0, 3, 4, 1}};
        edgeVertices_ = {{1, 2}, {2, 0}, {0, 1}};
        quadratures::triangleQuadrature(numberOfQuadraturePoints_, xsi, weight);
        break;
    }
    case PartitionOfUnity::Q4:
    {
        order_ = 1;
        const int numberOfQuadraturePoints_aux = 2;
        numberOfQuadraturePoints_ = pow(numberOfQuadraturePoints_aux, 2);
        numberOfNodes_ = 4;
        numberOfFaces_ = 4;
        numberOfEdges_ = 4;
        vtkCellType_ = VTK_LAGRANGE_QUADRILATERAL;
        vtkConnectivity_ = {0, 1, 2, 3};
        faceNodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        faceVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        edgeNodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        edgeVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        double **xsi_aux, *weight_aux;
        quadratures::lineQuadrature(numberOfQuadraturePoints_aux, xsi_aux, weight_aux);
        xsi = new double*[numberOfQuadraturePoints_];
        for (int i = 0; i < numberOfQuadraturePoints_; i++)
        {
            xsi[i] = new double[2];
        }
        weight = new double[numberOfQuadraturePoints_];
        int sum = -1;
        for (int j = 0; j < numberOfQuadraturePoints_aux; j++)
        {
            for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
            {
                ++sum;
                xsi[sum][0] = xsi_aux[i][0];
                xsi[sum][1] = xsi_aux[j][0];
                weight[sum] = weight_aux[i] * weight_aux[j];
            }
        }
        for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
        {
            delete[] xsi_aux[i];
        }
        delete[] xsi_aux;
        delete[] weight_aux;
        break;
    }
    case PartitionOfUnity::Q9:
    {
        order_ = 2;
        const int numberOfQuadraturePoints_aux = 3;
        numberOfQuadraturePoints_ = pow(numberOfQuadraturePoints_aux, 2);
        numberOfNodes_ = 9;
        numberOfFaces_ = 4;
        numberOfEdges_ = 4;
        vtkCellType_ = VTK_LAGRANGE_QUADRILATERAL;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 8};
        faceNodes_ = {{0, 4, 1}, {1, 5, 2}, {2, 6, 3}, {3, 7, 0}};
        faceVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        edgeNodes_ = {{0, 4, 1}, {1, 5, 2}, {2, 6, 3}, {3, 7, 0}};
        edgeVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        double **xsi_aux, *weight_aux;
        quadratures::lineQuadrature(numberOfQuadraturePoints_aux, xsi_aux, weight_aux);
        xsi = new double*[numberOfQuadraturePoints_];
        for (int i = 0; i < numberOfQuadraturePoints_; i++)
        {
            xsi[i] = new double[2];
        }
        weight = new double[numberOfQuadraturePoints_];
        int sum = -1;
        for (int j = 0; j < numberOfQuadraturePoints_aux; j++)
        {
            for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
            {
                ++sum;
                xsi[sum][0] = xsi_aux[i][0];
                xsi[sum][1] = xsi_aux[j][0];
                weight[sum] = weight_aux[i] * weight_aux[j];
            }
        }
        for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
        {
            delete[] xsi_aux[i];
        }
        delete[] xsi_aux;
        delete[] weight_aux;
        break;
    }
    case PartitionOfUnity::Q16:
    {
        order_ = 3;
        const int numberOfQuadraturePoints_aux = 4;
        numberOfQuadraturePoints_ = pow(numberOfQuadraturePoints_aux, 2);
        numberOfNodes_ = 16;
        numberOfFaces_ = 4;
        numberOfEdges_ = 4;
        vtkCellType_ = VTK_LAGRANGE_QUADRILATERAL;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8, 11, 10, 12, 13, 15, 14};
        faceNodes_ = {{0, 4, 5, 1}, {1, 6, 7, 2}, {2, 8, 9, 3}, {3, 10, 11, 0}};
        faceVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        edgeNodes_ = {{0, 4, 5, 1}, {1, 6, 7, 2}, {2, 8, 9, 3}, {3, 10, 11, 0}};
        edgeVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
        double **xsi_aux, *weight_aux;
        quadratures::lineQuadrature(numberOfQuadraturePoints_aux, xsi_aux, weight_aux);
        xsi = new double*[numberOfQuadraturePoints_];
        for (int i = 0; i < numberOfQuadraturePoints_; i++)
        {
            xsi[i] = new double[2];
        }
        weight = new double[numberOfQuadraturePoints_];
        int sum = -1;
        for (int j = 0; j < numberOfQuadraturePoints_aux; j++)
        {
            for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
            {
                ++sum;
                xsi[sum][0] = xsi_aux[i][0];
                xsi[sum][1] = xsi_aux[j][0];
                weight[sum] = weight_aux[i] * weight_aux[j];
            }
        }
        for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
        {
            delete[] xsi_aux[i];
        }
        delete[] xsi_aux;
        delete[] weight_aux;
        break;
    }
    default:
        std::cout << "Parametric surface element not implemented!\n";
        exit(EXIT_FAILURE);
        break;
    }
    quadraturePoints_.reserve(numberOfQuadraturePoints_);
    for (int i = 0; i < numberOfQuadraturePoints_; i++)
    {
        double *phi, **dphi_dxsi;
        getShapeFunctions(xsi[i], phi);
        getShapeFunctionsDerivatives(xsi[i], dphi_dxsi);
        
        quadraturePoints_.emplace_back(new QuadraturePoint(xsi[i], phi, dphi_dxsi, weight[i], 2, numberOfNodes_));
    }
    setNodalParametricCoordinates();
}

ParametricSurfaceElement::~ParametricSurfaceElement()
{
    for (QuadraturePoint* qp : quadraturePoints_)
    {
        delete qp;
    }
    for (int i = 0; i < 2; i++)
        delete[] nodalParametricCoordinates_[i];
    delete[] nodalParametricCoordinates_;
}

void ParametricSurfaceElement::setNodalParametricCoordinates()
{
    nodalParametricCoordinates_ = new double*[2];
    for (int i = 0; i < 2; i++)
        nodalParametricCoordinates_[i] = new double[numberOfNodes_];
    
    switch (elementType_)
    {
    case PartitionOfUnity::T3:
        nodalParametricCoordinates_[0][0] = 0.0; 
        nodalParametricCoordinates_[0][1] = 1.0; 
        nodalParametricCoordinates_[0][2] = 0.0;

        nodalParametricCoordinates_[1][0] = 0.0;
        nodalParametricCoordinates_[1][1] = 0.0;
        nodalParametricCoordinates_[1][2] = 1.0;
        break;
    case PartitionOfUnity::T6:
        nodalParametricCoordinates_[0][0] = 0.0; 
        nodalParametricCoordinates_[0][1] = 1.0; 
        nodalParametricCoordinates_[0][2] = 0.0; 
        nodalParametricCoordinates_[0][3] = 0.5; 
        nodalParametricCoordinates_[0][4] = 0.5; 
        nodalParametricCoordinates_[0][5] = 0.0; 

        nodalParametricCoordinates_[1][0] = 0.0;
        nodalParametricCoordinates_[1][1] = 0.0;
        nodalParametricCoordinates_[1][2] = 1.0;
        nodalParametricCoordinates_[1][3] = 0.0;
        nodalParametricCoordinates_[1][4] = 0.5;
        nodalParametricCoordinates_[1][5] = 0.5;
        break;
    case PartitionOfUnity::T10:
        nodalParametricCoordinates_[0][0] = 0.0;       
        nodalParametricCoordinates_[0][1] = 1.0;       
        nodalParametricCoordinates_[0][2] = 0.0;       
        nodalParametricCoordinates_[0][3] = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][4] = 2.0 / 3.0; 
        nodalParametricCoordinates_[0][5] = 2.0 / 3.0; 
        nodalParametricCoordinates_[0][6] = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][7] = 0.0;       
        nodalParametricCoordinates_[0][8] = 0.0;       
        nodalParametricCoordinates_[0][9] = 1.0 / 3.0; 

        nodalParametricCoordinates_[1][0] = 0.0;
        nodalParametricCoordinates_[1][1] = 0.0;
        nodalParametricCoordinates_[1][2] = 1.0;
        nodalParametricCoordinates_[1][3] = 0.0;
        nodalParametricCoordinates_[1][4] = 0.0;
        nodalParametricCoordinates_[1][5] = 1.0 / 3.0;
        nodalParametricCoordinates_[1][6] = 2.0 / 3.0;
        nodalParametricCoordinates_[1][7] = 2.0 / 3.0;
        nodalParametricCoordinates_[1][8] = 1.0 / 3.0;
        nodalParametricCoordinates_[1][9] = 1.0 / 3.0;
        break;
    case PartitionOfUnity::Q4:
        nodalParametricCoordinates_[0][0] = -1.0; 
        nodalParametricCoordinates_[0][1] =  1.0; 
        nodalParametricCoordinates_[0][2] =  1.0; 
        nodalParametricCoordinates_[0][3] = -1.0; 

        nodalParametricCoordinates_[1][0] = -1.0;
        nodalParametricCoordinates_[1][1] = -1.0;
        nodalParametricCoordinates_[1][2] =  1.0;
        nodalParametricCoordinates_[1][3] =  1.0;
        break;
    case PartitionOfUnity::Q9:
        nodalParametricCoordinates_[0][0] = -1.0; 
        nodalParametricCoordinates_[0][1] =  1.0; 
        nodalParametricCoordinates_[0][2] =  1.0; 
        nodalParametricCoordinates_[0][3] = -1.0; 
        nodalParametricCoordinates_[0][4] =  0.0; 
        nodalParametricCoordinates_[0][5] =  1.0; 
        nodalParametricCoordinates_[0][6] =  0.0; 
        nodalParametricCoordinates_[0][7] = -1.0; 
        nodalParametricCoordinates_[0][8] =  0.0; 

        nodalParametricCoordinates_[1][0] = -1.0;
        nodalParametricCoordinates_[1][1] = -1.0;
        nodalParametricCoordinates_[1][2] =  1.0;
        nodalParametricCoordinates_[1][3] =  1.0;
        nodalParametricCoordinates_[1][4] = -1.0;
        nodalParametricCoordinates_[1][5] =  0.0;
        nodalParametricCoordinates_[1][6] =  1.0;
        nodalParametricCoordinates_[1][7] =  0.0;
        nodalParametricCoordinates_[1][8] =  0.0;
        break;
    case PartitionOfUnity::Q16:
        nodalParametricCoordinates_[0][0]  = -1.0;       
        nodalParametricCoordinates_[0][1]  =  1.0;       
        nodalParametricCoordinates_[0][2]  =  1.0;       
        nodalParametricCoordinates_[0][3]  = -1.0;       
        nodalParametricCoordinates_[0][4]  = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][5]  =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][6]  =  1.0;       
        nodalParametricCoordinates_[0][7]  =  1.0;       
        nodalParametricCoordinates_[0][8]  =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][9]  = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][10] = -1.0;       
        nodalParametricCoordinates_[0][11] = -1.0;       
        nodalParametricCoordinates_[0][12] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][13] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][14] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][15] = -1.0 / 3.0; 

        nodalParametricCoordinates_[1][0]  = -1.0;
        nodalParametricCoordinates_[1][1]  = -1.0;
        nodalParametricCoordinates_[1][2]  =  1.0;
        nodalParametricCoordinates_[1][3]  =  1.0;
        nodalParametricCoordinates_[1][4]  = -1.0;
        nodalParametricCoordinates_[1][5]  = -1.0;
        nodalParametricCoordinates_[1][6]  = -1.0 / 3.0;
        nodalParametricCoordinates_[1][7]  =  1.0 / 3.0;
        nodalParametricCoordinates_[1][8]  =  1.0;
        nodalParametricCoordinates_[1][9]  =  1.0;
        nodalParametricCoordinates_[1][10] =  1.0 / 3.0;
        nodalParametricCoordinates_[1][11] = -1.0 / 3.0;
        nodalParametricCoordinates_[1][12] = -1.0 / 3.0;
        nodalParametricCoordinates_[1][13] = -1.0 / 3.0;
        nodalParametricCoordinates_[1][14] =  1.0 / 3.0;
        nodalParametricCoordinates_[1][15] =  1.0 / 3.0;
        break;
    default:
        break;
    }
}

void ParametricSurfaceElement::getShapeFunctions(double* xsi, double*& phi) const
{
    phi = new double[numberOfNodes_];
    double xsi1 = xsi[0];
    double xsi2 = xsi[1];

    switch (elementType_)
    {
    case PartitionOfUnity::T3:
        phi[0] = xsi1;
        phi[1] = xsi2;
	    phi[2] = 1.0 - xsi1 - xsi2;
        break;
    case PartitionOfUnity::T6:
        phi[0] = xsi1 * (2.0 * xsi1 - 1.0);
        phi[1] = xsi2 * (2.0 * xsi2 - 1.0);
        phi[2] = (xsi2 + xsi1 - 1.0) * (2.0 * xsi2 + 2.0 * xsi1 - 1.0);
        phi[3] = 4.0 * xsi1 * xsi2;
        phi[4] = -4.0 * xsi2 * (xsi2 + xsi1 - 1.0);
        phi[5] = -4.0 * xsi1 * (xsi2 + xsi1 - 1.0);
        break;
    case PartitionOfUnity::T10:
        phi[0] = (xsi1 * (3.0 * xsi1 - 2.0) * (3.0 * xsi1 - 1.0)) / 2.0;
        phi[1] = (xsi2 * (3.0 * xsi2 - 2.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi[2] = -((xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0) * (3.0 * xsi2 + 3.0 * xsi1 - 1.0)) / 2.0;
        phi[3] = (9.0 * xsi1 * xsi2 * (3.0 * xsi1 - 1.0)) / 2.0;
        phi[4] = (9.0 * xsi1 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        phi[5] = -(9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 - 1.0)) / 2.0;
        phi[6] = (9.0 * xsi2 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi[7] = (9.0 * xsi1 * (xsi2 + xsi1 - 1.0) * (3.0 * xsi2 + 3.0 * xsi1 - 2.0)) / 2.0;
        phi[8] = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0) * (xsi2 + xsi1 - 1.0)) / 2.0;
        phi[9] = -27.0 * xsi1 * xsi2 * (xsi2 + xsi1 - 1.0);
        break;
    }
}

void ParametricSurfaceElement::getShapeFunctionsDerivatives(double* xsi, double**& dphi_dxsi) const
{
    dphi_dxsi = new double*[2];
    for (int i = 0; i < 2; i++)
    {
        dphi_dxsi[i] = new double[numberOfNodes_];
    }
    double xsi1 = xsi[0];
    double xsi2 = xsi[1];

    switch (elementType_)
    {
    case PartitionOfUnity::T3:
        dphi_dxsi[0][0] = 1.0;    
	    dphi_dxsi[0][1] = 0.0;    
        dphi_dxsi[0][2] = -1.0;   
        
        dphi_dxsi[1][0] = 0.0;
        dphi_dxsi[1][1] = 1.0;
        dphi_dxsi[1][2] = -1.0;
        break;
    case PartitionOfUnity::T6:
        dphi_dxsi[0][0] = 4.0 * xsi1 - 1.0;
        dphi_dxsi[0][1] = 0.0;
        dphi_dxsi[0][2] = 4.0 * xsi2 + 4.0 * xsi1 - 3.0;
        dphi_dxsi[0][3] = 4.0 * xsi2;
        dphi_dxsi[0][4] = -4.0 * xsi2;
        dphi_dxsi[0][5] = -4.0 * (xsi2 + 2.0 * xsi1 - 1.0);

        dphi_dxsi[1][0] = 0.0;
        dphi_dxsi[1][1] = 4.0 * xsi2 - 1.0;
        dphi_dxsi[1][2] = 4.0 * xsi2 + 4.0 * xsi1 - 3.0;
        dphi_dxsi[1][3] = 4.0 * xsi1;
        dphi_dxsi[1][4] = -4.0 * (2.0 * xsi2 + xsi1 - 1.0);
        dphi_dxsi[1][5] = -4.0 * xsi1;
        break;
    case PartitionOfUnity::T10:
        dphi_dxsi[0][0] = (27.0 * xsi1 * xsi1 - 18.0 * xsi1 + 2.0) / 2.0;
        dphi_dxsi[0][1] = 0.0;
        dphi_dxsi[0][2] = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi_dxsi[0][3] = (9.0 * xsi2 * (6.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi[0][4] = (9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi[0][5] = -(9.0 * xsi2 * (3.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi[0][6] = (9.0 * xsi2 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi_dxsi[0][7] = (9.0 * (3.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 5.0 * xsi2 + 9.0 * xsi1 * xsi1 - 10.0 * xsi1 + 2.0)) / 2.0;
        dphi_dxsi[0][8] = -(9.0 * (6.0 * xsi1 * xsi2 - xsi2 + 9.0 * xsi1 * xsi1 - 8.0 * xsi1 + 1.0)) / 2.0;
        dphi_dxsi[0][9] = -27.0 * xsi2 * (xsi2 + 2.0 * xsi1 - 1.0);

        dphi_dxsi[1][0] = 0.0;
        dphi_dxsi[1][1] = (27.0 * xsi2 * xsi2 - 18.0 * xsi2 + 2) / 2.0;
        dphi_dxsi[1][2] = -(27.0 * xsi2 * xsi2 + 54.0 * xsi1 * xsi2 - 36.0 * xsi2 + 27.0 * xsi1 * xsi1 - 36.0 * xsi1 + 11.0) / 2.0;
        dphi_dxsi[1][3] = (9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi[1][4] = (9.0 * xsi1 * (6.0 * xsi2 - 1.0)) / 2.0;
        dphi_dxsi[1][5] = -(9.0 * (9.0 * xsi2 * xsi2 + 6.0 * xsi1 * xsi2 - 8.0 * xsi2 - xsi1 + 1.0)) / 2.0;
        dphi_dxsi[1][6] = (9.0 * (9.0 * xsi2 * xsi2 + 12.0 * xsi1 * xsi2 - 10.0 * xsi2 + 3.0 * xsi1 * xsi1 - 5.0 * xsi1 + 2.0)) / 2.0;
        dphi_dxsi[1][7] = (9.0 * xsi1 * (6.0 * xsi2 + 6.0 * xsi1 - 5.0)) / 2.0;
        dphi_dxsi[1][8] = -(9.0 * xsi1 * (3.0 * xsi1 - 1.0)) / 2.0;
        dphi_dxsi[1][9] = -27.0 * xsi1 * (2.0 * xsi2 + xsi1 - 1.0);
        break;
    default:
        break;
    }
}