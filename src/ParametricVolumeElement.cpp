#include "ParametricVolumeElement.h"

ParametricVolumeElement ParametricVolumeElement::TET4(PartitionOfUnity::TET4);
ParametricVolumeElement ParametricVolumeElement::TET10(PartitionOfUnity::TET10);
ParametricVolumeElement ParametricVolumeElement::TET20(PartitionOfUnity::TET20);
ParametricVolumeElement ParametricVolumeElement::HEX8(PartitionOfUnity::HEX8);
ParametricVolumeElement ParametricVolumeElement::HEX27(PartitionOfUnity::HEX27);
ParametricVolumeElement ParametricVolumeElement::HEX64(PartitionOfUnity::HEX64);

ParametricVolumeElement::ParametricVolumeElement(const PartitionOfUnity elementType)
    : ParametricElement(elementType)
{
    double **xsi, *weight;
    switch (elementType)
    {
    case PartitionOfUnity::TET4:
    {
        order_ = 1;
        numberOfQuadraturePoints_ = 1;
        numberOfNodes_ = 4;
        numberOfFaces_ = 4;
        numberOfEdges_ = 6;
        vtkCellType_ = VTK_LAGRANGE_TETRAHEDRON;
        vtkConnectivity_ = {0, 1, 2, 3};
        faceNodes_ = {{1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {0, 2, 1}};
        faceVertices_ = {{1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {0, 2, 1}};
        edgeNodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 0}};
        edgeVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 0}};
        quadratures::tetrahedronQuadrature(numberOfQuadraturePoints_, xsi, weight);
        break;
    }
    case PartitionOfUnity::TET10:
    {
        order_ = 2;
        numberOfQuadraturePoints_ = 14;
        numberOfNodes_ = 10;
        numberOfFaces_ = 4;
        numberOfEdges_ = 6;
        vtkCellType_ = VTK_LAGRANGE_TETRAHEDRON;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8};
        faceNodes_ = {{1, 5, 2, 8, 3, 9}, {2, 6, 0, 7, 3, 8}, {0, 4, 1, 9, 3, 7}, {0, 6, 2, 5, 1, 4}};
        faceVertices_ = {{1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {0, 2, 1}};
        edgeNodes_ = {{0, 4, 1}, {1, 5, 2}, {2, 8, 3}, {3, 7, 0}, {1, 9, 3}, {2, 6, 0}};
        edgeVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 0}};
        quadratures::tetrahedronQuadrature(numberOfQuadraturePoints_, xsi, weight);
        break;
    }
    case PartitionOfUnity::TET20:
    {
        order_ = 3;
        numberOfQuadraturePoints_ = 24;
        numberOfNodes_ = 20;
        numberOfFaces_ = 4;
        numberOfEdges_ = 6;
        vtkCellType_ = VTK_LAGRANGE_TETRAHEDRON;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 10, 15, 14, 13, 12, 17, 19, 18, 16};
        faceNodes_ = {{}, {}, {}, {}};
        faceVertices_ = {{1, 2, 3}, {2, 0, 3}, {0, 1, 3}, {0, 2, 1}};
        edgeNodes_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 0}};
        edgeVertices_ = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {1, 3}, {2, 0}};
        quadratures::tetrahedronQuadrature(numberOfQuadraturePoints_, xsi, weight);
        break;
    }
    case PartitionOfUnity::HEX8:
    {
        order_ = 1;
        const int numberOfQuadraturePoints_aux = 2;
        numberOfQuadraturePoints_ = pow(numberOfQuadraturePoints_aux, 3);
        numberOfNodes_ = 8;
        numberOfFaces_ = 6;
        numberOfEdges_ = 12;
        vtkCellType_ = VTK_LAGRANGE_HEXAHEDRON;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7};
        faceNodes_ = {{}, {}, {}, {}, {}, {}, {}};
        faceVertices_ = {{}, {}, {}, {}, {}, {}, {}};
        double **xsi_aux, *weight_aux;
        quadratures::lineQuadrature(numberOfQuadraturePoints_aux, xsi_aux, weight_aux);
        xsi = new double*[numberOfQuadraturePoints_];
        for (int i = 0; i < numberOfQuadraturePoints_; i++)
        {
            xsi[i] = new double[3];
        }
        weight = new double[numberOfQuadraturePoints_];
        int sum = -1;
        for (int k = 0; k < numberOfQuadraturePoints_aux; k++)
        {
            for (int j = 0; j < numberOfQuadraturePoints_aux; j++)
            {
                for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
                {
                    ++sum;
                    xsi[sum][0] = xsi_aux[i][0];
                    xsi[sum][1] = xsi_aux[j][0];
                    xsi[sum][2] = xsi_aux[k][0];
                    weight[sum] = weight_aux[i] * weight_aux[j];
                }
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
    case PartitionOfUnity::HEX27:
    {
        order_ = 2;
        const int numberOfQuadraturePoints_aux = 3;
        numberOfQuadraturePoints_ = pow(numberOfQuadraturePoints_aux, 3);
        numberOfNodes_ = 27;
        numberOfFaces_ = 6;
        numberOfEdges_ = 12;
        vtkCellType_ = VTK_LAGRANGE_HEXAHEDRON;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 13, 9, 16, 18, 19, 17 ,10, 12, 15, 14,
                            22, 23, 21, 24, 20, 25, 26};
        faceNodes_ = {{}, {}, {}, {}, {}, {}, {}};
        faceVertices_ = {{}, {}, {}, {}, {}, {}, {}};
        double **xsi_aux, *weight_aux;
        quadratures::lineQuadrature(numberOfQuadraturePoints_aux, xsi_aux, weight_aux);
        xsi = new double*[numberOfQuadraturePoints_];
        for (int i = 0; i < numberOfQuadraturePoints_; i++)
        {
            xsi[i] = new double[3];
        }
        weight = new double[numberOfQuadraturePoints_];
        int sum = -1;
        for (int k = 0; k < numberOfQuadraturePoints_aux; k++)
        {
            for (int j = 0; j < numberOfQuadraturePoints_aux; j++)
            {
                for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
                {
                    ++sum;
                    xsi[sum][0] = xsi_aux[i][0];
                    xsi[sum][1] = xsi_aux[j][0];
                    xsi[sum][2] = xsi_aux[k][0];
                    weight[sum] = weight_aux[i] * weight_aux[j];
                }
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
    case PartitionOfUnity::HEX64:
    {
        order_ = 3;
        const int numberOfQuadraturePoints_aux = 4;
        numberOfQuadraturePoints_ = pow(numberOfQuadraturePoints_aux, 3);
        numberOfNodes_ = 64;
        numberOfFaces_ = 6;
        numberOfEdges_ = 12;
        vtkCellType_ = VTK_LAGRANGE_HEXAHEDRON;
        vtkConnectivity_ = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15, 19, 18, 10, 11 ,24, 25, 28, 29,
                            31, 30, 26, 27, 12, 13, 16, 17, 22, 23, 20, 21, 40, 43, 41, 42, 44,
                            45, 47, 46, 36, 37, 39, 38, 49, 48, 50, 51, 32, 35, 33, 34, 52, 53,
                            55, 54, 56, 57, 59, 58, 60, 61, 63, 62};
        faceNodes_ = {{}, {}, {}, {}, {}, {}, {}};
        faceVertices_ = {{}, {}, {}, {}, {}, {}, {}};
        double **xsi_aux, *weight_aux;
        quadratures::lineQuadrature(numberOfQuadraturePoints_aux, xsi_aux, weight_aux);
        xsi = new double*[numberOfQuadraturePoints_];
        for (int i = 0; i < numberOfQuadraturePoints_; i++)
        {
            xsi[i] = new double[3];
        }
        weight = new double[numberOfQuadraturePoints_];
        int sum = -1;
        for (int k = 0; k < numberOfQuadraturePoints_aux; k++)
        {
            for (int j = 0; j < numberOfQuadraturePoints_aux; j++)
            {
                for (int i = 0; i < numberOfQuadraturePoints_aux; i++)
                {
                    ++sum;
                    xsi[sum][0] = xsi_aux[i][0];
                    xsi[sum][1] = xsi_aux[j][0];
                    xsi[sum][2] = xsi_aux[k][0];
                    weight[sum] = weight_aux[i] * weight_aux[j];
                }
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
        std::cout << "Parametric volume element not implemented!\n";
        exit(EXIT_FAILURE);
    }
    quadraturePoints_.reserve(numberOfQuadraturePoints_);
    for (int i = 0; i < numberOfQuadraturePoints_; i++)
    {
        double *phi, **dphi_dxsi;
        getShapeFunctions(xsi[i], phi);
        getShapeFunctionsDerivatives(xsi[i], dphi_dxsi);
        
        quadraturePoints_.emplace_back(new QuadraturePoint(xsi[i], phi, dphi_dxsi, weight[i], 3, numberOfNodes_));
    }
    setNodalParametricCoordinates();
}

ParametricVolumeElement::~ParametricVolumeElement()
{
    for (QuadraturePoint* qp : quadraturePoints_)
    {
        delete qp;
    }
    for (int i = 0; i < 3; i++)
        delete[] nodalParametricCoordinates_[i];
    delete[] nodalParametricCoordinates_;
}

void ParametricVolumeElement::setNodalParametricCoordinates()
{
    nodalParametricCoordinates_ = new double*[3];
    for (int i = 0; i < 3; i++)
        nodalParametricCoordinates_[i] = new double[numberOfNodes_];

    switch (elementType_)
    {
    case PartitionOfUnity::TET4:
        nodalParametricCoordinates_[0][0] = 0.0; 
        nodalParametricCoordinates_[0][1] = 1.0; 
        nodalParametricCoordinates_[0][2] = 0.0; 
        nodalParametricCoordinates_[0][3] = 0.0; 

        nodalParametricCoordinates_[1][0] = 0.0; 
        nodalParametricCoordinates_[1][1] = 0.0; 
        nodalParametricCoordinates_[1][2] = 1.0; 
        nodalParametricCoordinates_[1][3] = 0.0; 

        nodalParametricCoordinates_[2][0] = 0.0;
        nodalParametricCoordinates_[2][1] = 0.0;
        nodalParametricCoordinates_[2][2] = 0.0;
        nodalParametricCoordinates_[2][3] = 1.0;
        break;
    case PartitionOfUnity::TET10:
        nodalParametricCoordinates_[0][0] = 0.0; 
        nodalParametricCoordinates_[0][1] = 1.0; 
        nodalParametricCoordinates_[0][2] = 0.0; 
        nodalParametricCoordinates_[0][3] = 0.0; 
        nodalParametricCoordinates_[0][4] = 0.5; 
        nodalParametricCoordinates_[0][5] = 0.5; 
        nodalParametricCoordinates_[0][6] = 0.0; 
        nodalParametricCoordinates_[0][7] = 0.0; 
        nodalParametricCoordinates_[0][8] = 0.0; 
        nodalParametricCoordinates_[0][9] = 0.5; 

        nodalParametricCoordinates_[1][0] = 0.0; 
        nodalParametricCoordinates_[1][1] = 0.0; 
        nodalParametricCoordinates_[1][2] = 1.0; 
        nodalParametricCoordinates_[1][3] = 0.0; 
        nodalParametricCoordinates_[1][4] = 0.0; 
        nodalParametricCoordinates_[1][5] = 0.5; 
        nodalParametricCoordinates_[1][6] = 0.5; 
        nodalParametricCoordinates_[1][7] = 0.0; 
        nodalParametricCoordinates_[1][8] = 0.5; 
        nodalParametricCoordinates_[1][9] = 0.0; 

        nodalParametricCoordinates_[2][0] = 0.0;
        nodalParametricCoordinates_[2][1] = 0.0;
        nodalParametricCoordinates_[2][2] = 0.0;
        nodalParametricCoordinates_[2][3] = 1.0;
        nodalParametricCoordinates_[2][4] = 0.0;
        nodalParametricCoordinates_[2][5] = 0.0;
        nodalParametricCoordinates_[2][6] = 0.0;
        nodalParametricCoordinates_[2][7] = 0.5;
        nodalParametricCoordinates_[2][8] = 0.5;
        nodalParametricCoordinates_[2][9] = 0.5;
        break;
    case PartitionOfUnity::TET20:
        nodalParametricCoordinates_[0][0]  = 0.0; 
        nodalParametricCoordinates_[0][1]  = 1.0; 
        nodalParametricCoordinates_[0][2]  = 0.0; 
        nodalParametricCoordinates_[0][3]  = 0.0; 
        nodalParametricCoordinates_[0][4]  = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][5]  = 2.0 / 3.0; 
        nodalParametricCoordinates_[0][6]  = 2.0 / 3.0; 
        nodalParametricCoordinates_[0][7]  = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][8]  = 0.0; 
        nodalParametricCoordinates_[0][9]  = 0.0; 
        nodalParametricCoordinates_[0][10] = 0.0; 
        nodalParametricCoordinates_[0][11] = 0.0; 
        nodalParametricCoordinates_[0][12] = 0.0; 
        nodalParametricCoordinates_[0][13] = 0.0; 
        nodalParametricCoordinates_[0][14] = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][15] = 2.0 / 3.0; 
        nodalParametricCoordinates_[0][16] = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][17] = 1.0 / 3.0; 
        nodalParametricCoordinates_[0][18] = 0.0; 
        nodalParametricCoordinates_[0][19] = 1.0 / 3.0; 

        nodalParametricCoordinates_[1][0]  = 0.0; 
        nodalParametricCoordinates_[1][1]  = 0.0; 
        nodalParametricCoordinates_[1][2]  = 1.0; 
        nodalParametricCoordinates_[1][3]  = 0.0; 
        nodalParametricCoordinates_[1][4]  = 0.0; 
        nodalParametricCoordinates_[1][5]  = 0.0; 
        nodalParametricCoordinates_[1][6]  = 1.0 / 3.0; 
        nodalParametricCoordinates_[1][7]  = 2.0 / 3.0; 
        nodalParametricCoordinates_[1][8]  = 2.0 / 3.0; 
        nodalParametricCoordinates_[1][9]  = 1.0 / 3.0; 
        nodalParametricCoordinates_[1][10] = 0.0; 
        nodalParametricCoordinates_[1][11] = 0.0; 
        nodalParametricCoordinates_[1][12] = 1.0 / 3.0; 
        nodalParametricCoordinates_[1][13] = 2.0 / 3.0; 
        nodalParametricCoordinates_[1][14] = 0.0; 
        nodalParametricCoordinates_[1][15] = 0.0; 
        nodalParametricCoordinates_[1][16] = 1.0 / 3.0; 
        nodalParametricCoordinates_[1][17] = 0.0; 
        nodalParametricCoordinates_[1][18] = 1.0 / 3.0; 
        nodalParametricCoordinates_[1][19] = 1.0 / 3.0; 

        nodalParametricCoordinates_[2][1]  = 0.0;
        nodalParametricCoordinates_[2][2]  = 0.0;
        nodalParametricCoordinates_[2][3]  = 1.0;
        nodalParametricCoordinates_[2][4]  = 0.0;
        nodalParametricCoordinates_[2][5]  = 0.0;
        nodalParametricCoordinates_[2][6]  = 0.0;
        nodalParametricCoordinates_[2][7]  = 0.0;
        nodalParametricCoordinates_[2][8]  = 0.0;
        nodalParametricCoordinates_[2][9]  = 0.0;
        nodalParametricCoordinates_[2][10] = 2.0 / 3.0;
        nodalParametricCoordinates_[2][11] = 1.0 / 3.0;
        nodalParametricCoordinates_[2][12] = 2.0 / 3.0;
        nodalParametricCoordinates_[2][13] = 1.0 / 3.0;
        nodalParametricCoordinates_[2][14] = 2.0 / 3.0;
        nodalParametricCoordinates_[2][15] = 1.0 / 3.0;
        nodalParametricCoordinates_[2][16] = 0.0;
        nodalParametricCoordinates_[2][17] = 1.0 / 3.0;
        nodalParametricCoordinates_[2][18] = 1.0 / 3.0;
        nodalParametricCoordinates_[2][19] = 1.0 / 3.0;
        break;
    case PartitionOfUnity::HEX8:
        nodalParametricCoordinates_[0][0] = -1.0; 
        nodalParametricCoordinates_[0][1] =  1.0; 
        nodalParametricCoordinates_[0][2] =  1.0; 
        nodalParametricCoordinates_[0][3] = -1.0; 
        nodalParametricCoordinates_[0][4] = -1.0; 
        nodalParametricCoordinates_[0][5] =  1.0; 
        nodalParametricCoordinates_[0][6] =  1.0; 
        nodalParametricCoordinates_[0][7] = -1.0; 

        nodalParametricCoordinates_[1][0] = -1.0; 
        nodalParametricCoordinates_[1][1] = -1.0; 
        nodalParametricCoordinates_[1][2] =  1.0; 
        nodalParametricCoordinates_[1][3] =  1.0; 
        nodalParametricCoordinates_[1][4] = -1.0; 
        nodalParametricCoordinates_[1][5] = -1.0; 
        nodalParametricCoordinates_[1][6] =  1.0; 
        nodalParametricCoordinates_[1][7] =  1.0; 

        nodalParametricCoordinates_[2][0] = -1.0;
        nodalParametricCoordinates_[2][1] = -1.0;
        nodalParametricCoordinates_[2][2] = -1.0;
        nodalParametricCoordinates_[2][3] = -1.0;
        nodalParametricCoordinates_[2][4] =  1.0;
        nodalParametricCoordinates_[2][5] =  1.0;
        nodalParametricCoordinates_[2][6] =  1.0;
        nodalParametricCoordinates_[2][7] =  1.0;
        break;
    case PartitionOfUnity::HEX27:
        nodalParametricCoordinates_[0][0]  = -1.0; 
        nodalParametricCoordinates_[0][1]  =  1.0; 
        nodalParametricCoordinates_[0][2]  =  1.0; 
        nodalParametricCoordinates_[0][3]  = -1.0; 
        nodalParametricCoordinates_[0][4]  = -1.0; 
        nodalParametricCoordinates_[0][5]  =  1.0; 
        nodalParametricCoordinates_[0][6]  =  1.0; 
        nodalParametricCoordinates_[0][7]  = -1.0; 
        nodalParametricCoordinates_[0][8]  =  0.0; 
        nodalParametricCoordinates_[0][9]  = -1.0; 
        nodalParametricCoordinates_[0][10] = -1.0; 
        nodalParametricCoordinates_[0][11] =  1.0; 
        nodalParametricCoordinates_[0][12] =  1.0; 
        nodalParametricCoordinates_[0][13] =  0.0; 
        nodalParametricCoordinates_[0][14] =  1.0; 
        nodalParametricCoordinates_[0][15] = -1.0; 
        nodalParametricCoordinates_[0][16] =  0.0; 
        nodalParametricCoordinates_[0][17] = -1.0; 
        nodalParametricCoordinates_[0][18] =  1.0; 
        nodalParametricCoordinates_[0][19] =  0.0; 
        nodalParametricCoordinates_[0][20] =  0.0; 
        nodalParametricCoordinates_[0][21] =  0.0; 
        nodalParametricCoordinates_[0][22] = -1.0; 
        nodalParametricCoordinates_[0][23] =  1.0; 
        nodalParametricCoordinates_[0][24] =  0.0; 
        nodalParametricCoordinates_[0][25] =  0.0; 
        nodalParametricCoordinates_[0][26] =  0.0; 

        nodalParametricCoordinates_[1][0]  = -1.0; 
        nodalParametricCoordinates_[1][1]  = -1.0; 
        nodalParametricCoordinates_[1][2]  =  1.0; 
        nodalParametricCoordinates_[1][3]  =  1.0; 
        nodalParametricCoordinates_[1][4]  = -1.0; 
        nodalParametricCoordinates_[1][5]  = -1.0; 
        nodalParametricCoordinates_[1][6]  =  1.0; 
        nodalParametricCoordinates_[1][7]  =  1.0; 
        nodalParametricCoordinates_[1][8]  = -1.0; 
        nodalParametricCoordinates_[1][9]  =  0.0; 
        nodalParametricCoordinates_[1][10] = -1.0; 
        nodalParametricCoordinates_[1][11] =  0.0; 
        nodalParametricCoordinates_[1][12] = -1.0; 
        nodalParametricCoordinates_[1][13] =  1.0; 
        nodalParametricCoordinates_[1][14] =  1.0; 
        nodalParametricCoordinates_[1][15] =  1.0; 
        nodalParametricCoordinates_[1][16] = -1.0; 
        nodalParametricCoordinates_[1][17] =  0.0; 
        nodalParametricCoordinates_[1][18] =  0.0; 
        nodalParametricCoordinates_[1][19] =  1.0; 
        nodalParametricCoordinates_[1][20] =  0.0; 
        nodalParametricCoordinates_[1][21] = -1.0; 
        nodalParametricCoordinates_[1][22] =  0.0; 
        nodalParametricCoordinates_[1][23] =  0.0; 
        nodalParametricCoordinates_[1][24] =  1.0; 
        nodalParametricCoordinates_[1][25] =  0.0; 
        nodalParametricCoordinates_[1][26] =  0.0; 

        nodalParametricCoordinates_[2][0]  = -1.0;
        nodalParametricCoordinates_[2][1]  = -1.0;
        nodalParametricCoordinates_[2][2]  = -1.0;
        nodalParametricCoordinates_[2][3]  = -1.0;
        nodalParametricCoordinates_[2][4]  =  1.0;
        nodalParametricCoordinates_[2][5]  =  1.0;
        nodalParametricCoordinates_[2][6]  =  1.0;
        nodalParametricCoordinates_[2][7]  =  1.0;
        nodalParametricCoordinates_[2][8]  = -1.0;
        nodalParametricCoordinates_[2][9]  = -1.0;
        nodalParametricCoordinates_[2][10] =  0.0;
        nodalParametricCoordinates_[2][11] = -1.0;
        nodalParametricCoordinates_[2][12] =  0.0;
        nodalParametricCoordinates_[2][13] = -1.0;
        nodalParametricCoordinates_[2][14] =  0.0;
        nodalParametricCoordinates_[2][15] =  0.0;
        nodalParametricCoordinates_[2][16] =  1.0;
        nodalParametricCoordinates_[2][17] =  1.0;
        nodalParametricCoordinates_[2][18] =  1.0;
        nodalParametricCoordinates_[2][19] =  1.0;
        nodalParametricCoordinates_[2][20] = -1.0;
        nodalParametricCoordinates_[2][21] =  0.0;
        nodalParametricCoordinates_[2][22] =  0.0;
        nodalParametricCoordinates_[2][23] =  0.0;
        nodalParametricCoordinates_[2][24] =  0.0;
        nodalParametricCoordinates_[2][25] =  1.0;
        nodalParametricCoordinates_[2][26] =  0.0;
        break;
    case PartitionOfUnity::HEX64:
        nodalParametricCoordinates_[0][0]  = -1.0; 
        nodalParametricCoordinates_[0][1]  = 1.0; 
        nodalParametricCoordinates_[0][2]  = 1.0; 
        nodalParametricCoordinates_[0][3]  = -1.0; 
        nodalParametricCoordinates_[0][4]  = -1.0; 
        nodalParametricCoordinates_[0][5]  = 1.0; 
        nodalParametricCoordinates_[0][6]  = 1.0; 
        nodalParametricCoordinates_[0][7]  = -1.0; 
        nodalParametricCoordinates_[0][8]  = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][9]  =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][10] = -1.0; 
        nodalParametricCoordinates_[0][11] = -1.0; 
        nodalParametricCoordinates_[0][12] = -1.0; 
        nodalParametricCoordinates_[0][13] = -1.0; 
        nodalParametricCoordinates_[0][14] = 1.0; 
        nodalParametricCoordinates_[0][15] = 1.0; 
        nodalParametricCoordinates_[0][16] = 1.0; 
        nodalParametricCoordinates_[0][17] = 1.0; 
        nodalParametricCoordinates_[0][18] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][19] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][20] = 1.0; 
        nodalParametricCoordinates_[0][21] = 1.0; 
        nodalParametricCoordinates_[0][22] = -1.0; 
        nodalParametricCoordinates_[0][23] = -1.0; 
        nodalParametricCoordinates_[0][24] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][25] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][26] = -1.0; 
        nodalParametricCoordinates_[0][27] = -1.0; 
        nodalParametricCoordinates_[0][28] = 1.0; 
        nodalParametricCoordinates_[0][29] = 1.0; 
        nodalParametricCoordinates_[0][30] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][31] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][32] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][33] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][34] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][35] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][36] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][37] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][38] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][39] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][40] = -1.0; 
        nodalParametricCoordinates_[0][41] = -1.0; 
        nodalParametricCoordinates_[0][42] = -1.0; 
        nodalParametricCoordinates_[0][43] = -1.0; 
        nodalParametricCoordinates_[0][44] = 1.0; 
        nodalParametricCoordinates_[0][45] = 1.0; 
        nodalParametricCoordinates_[0][46] = 1.0; 
        nodalParametricCoordinates_[0][47] = 1.0; 
        nodalParametricCoordinates_[0][48] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][49] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][50] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][51] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][52] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][53] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][54] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][55] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][56] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][57] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][58] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][59] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][60] = -1.0 / 3.0; 
        nodalParametricCoordinates_[0][61] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][62] =  1.0 / 3.0; 
        nodalParametricCoordinates_[0][63] = -1.0 / 3.0;

        nodalParametricCoordinates_[1][0]  = -1.0; 
        nodalParametricCoordinates_[1][1]  = -1.0; 
        nodalParametricCoordinates_[1][2]  = 1.0; 
        nodalParametricCoordinates_[1][3]  = 1.0; 
        nodalParametricCoordinates_[1][4]  = -1.0; 
        nodalParametricCoordinates_[1][5]  = -1.0; 
        nodalParametricCoordinates_[1][6]  = 1.0; 
        nodalParametricCoordinates_[1][7]  = 1.0; 
        nodalParametricCoordinates_[1][8]  = -1.0; 
        nodalParametricCoordinates_[1][9]  = -1.0; 
        nodalParametricCoordinates_[1][10] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][11] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][12] = -1.0; 
        nodalParametricCoordinates_[1][13] = -1.0; 
        nodalParametricCoordinates_[1][14] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][15] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][16] = -1.0; 
        nodalParametricCoordinates_[1][17] = -1.0; 
        nodalParametricCoordinates_[1][18] = 1.0; 
        nodalParametricCoordinates_[1][19] = 1.0; 
        nodalParametricCoordinates_[1][20] = 1.0; 
        nodalParametricCoordinates_[1][21] = 1.0; 
        nodalParametricCoordinates_[1][22] = 1.0; 
        nodalParametricCoordinates_[1][23] = 1.0; 
        nodalParametricCoordinates_[1][24] = -1.0; 
        nodalParametricCoordinates_[1][25] = -1.0; 
        nodalParametricCoordinates_[1][26] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][27] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][28] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][29] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][30] = 1.0; 
        nodalParametricCoordinates_[1][31] = 1.0; 
        nodalParametricCoordinates_[1][32] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][33] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][34] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][35] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][36] = -1.0; 
        nodalParametricCoordinates_[1][37] = -1.0; 
        nodalParametricCoordinates_[1][38] = -1.0; 
        nodalParametricCoordinates_[1][39] = -1.0; 
        nodalParametricCoordinates_[1][40] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][41] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][42] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][43] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][44] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][45] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][46] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][47] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][48] = 1.0; 
        nodalParametricCoordinates_[1][49] = 1.0; 
        nodalParametricCoordinates_[1][50] = 1.0; 
        nodalParametricCoordinates_[1][51] = 1.0; 
        nodalParametricCoordinates_[1][52] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][53] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][54] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][55] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][56] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][57] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][58] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][59] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][60] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][61] = -1.0 / 3.0; 
        nodalParametricCoordinates_[1][62] =  1.0 / 3.0; 
        nodalParametricCoordinates_[1][63] =  1.0 / 3.0;

        nodalParametricCoordinates_[2][0]  = -1.0;
        nodalParametricCoordinates_[2][1]  = -1.0;
        nodalParametricCoordinates_[2][2]  = -1.0;
        nodalParametricCoordinates_[2][3]  = -1.0;
        nodalParametricCoordinates_[2][4]  = 1.0;
        nodalParametricCoordinates_[2][5]  = 1.0;
        nodalParametricCoordinates_[2][6]  = 1.0;
        nodalParametricCoordinates_[2][7]  = 1.0;
        nodalParametricCoordinates_[2][8]  = -1.0;
        nodalParametricCoordinates_[2][9]  = -1.0;
        nodalParametricCoordinates_[2][10] = -1.0;
        nodalParametricCoordinates_[2][11] = -1.0;
        nodalParametricCoordinates_[2][12] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][13] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][14] = -1.0;
        nodalParametricCoordinates_[2][15] = -1.0;
        nodalParametricCoordinates_[2][16] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][17] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][18] = -1.0;
        nodalParametricCoordinates_[2][19] = -1.0;
        nodalParametricCoordinates_[2][20] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][21] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][22] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][23] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][24] = 1.0;
        nodalParametricCoordinates_[2][25] = 1.0;
        nodalParametricCoordinates_[2][26] = 1.0;
        nodalParametricCoordinates_[2][27] = 1.0;
        nodalParametricCoordinates_[2][28] = 1.0;
        nodalParametricCoordinates_[2][29] = 1.0;
        nodalParametricCoordinates_[2][30] = 1.0;
        nodalParametricCoordinates_[2][31] = 1.0;
        nodalParametricCoordinates_[2][32] = -1.0;
        nodalParametricCoordinates_[2][33] = -1.0;
        nodalParametricCoordinates_[2][34] = -1.0;
        nodalParametricCoordinates_[2][35] = -1.0;
        nodalParametricCoordinates_[2][36] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][37] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][38] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][39] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][40] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][41] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][42] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][43] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][44] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][45] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][46] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][47] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][48] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][49] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][50] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][51] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][52] = 1.0;
        nodalParametricCoordinates_[2][53] = 1.0;
        nodalParametricCoordinates_[2][54] = 1.0;
        nodalParametricCoordinates_[2][55] = 1.0;
        nodalParametricCoordinates_[2][56] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][57] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][58] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][59] = -1.0 / 3.0;
        nodalParametricCoordinates_[2][60] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][61] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][62] =  1.0 / 3.0;
        nodalParametricCoordinates_[2][63] =  1.0 / 3.0; 
        break;
    default:
        break;
    }
}
    

void ParametricVolumeElement::getShapeFunctions(double* xsi, double*& phi) const
{
    phi = new double[numberOfNodes_];
    double xsi1 = xsi[0];
    double xsi2 = xsi[1];
    double xsi3 = xsi[2];

    switch (elementType_)
    {
    case PartitionOfUnity::TET4:
        phi[0] = 1.0 - xsi1 - xsi2 - xsi3;
        phi[1] = xsi1;
        phi[2] = xsi2;
        phi[3] = xsi3;
        break;
    case PartitionOfUnity::TET10:
        phi[0] = 1.0 - 3.0 * xsi1 + 2.0 * xsi1 * xsi1 - 3.0 * xsi2 + 4.0 * xsi1 * xsi2 + 2.0 * xsi2 * xsi2 - 3.0 * xsi3 + 4.0 * xsi1 * xsi3 + 4.0 * xsi2 * xsi3 + 2.0 * xsi3 * xsi3;
        phi[1] = (2.0 * xsi1 - 1.0) * xsi1;
        phi[2] = (2.0 * xsi2 - 1.0) * xsi2;
        phi[3] = (2.0 * xsi3 - 1.0) * xsi3;
        phi[4] = 4.0 * xsi1 - 4.0 * xsi1 * xsi1 - 4.0 * xsi1 * xsi2 - 4.0 * xsi1 * xsi3;
        phi[5] = 4.0 * xsi1 * xsi2;
        phi[6] = 4.0 * xsi2 - 4.0 * xsi1 * xsi2 - 4.0 * xsi2 * xsi2 - 4.0 * xsi2 * xsi3;
        phi[7] = 4.0 * xsi3 - 4.0 * xsi1 * xsi3 - 4.0 * xsi2 * xsi3 - 4.0 * xsi3 * xsi3;
        phi[8] = 4.0 * xsi2 * xsi3;
        phi[9] = 4.0 * xsi1 * xsi3;
        break;
    case PartitionOfUnity::TET20:
        double xsi4 = 1.0 - xsi1 - xsi2 - xsi3;

        phi[0] = (1.0 / 2.0) * (3.0 * xsi4 - 1.0) * (3.0 * xsi4 - 2.0) * xsi4;
        phi[1] = (1.0 / 2.0) * (3.0 * xsi1 - 1.0) * (3.0 * xsi1 - 2.0) * xsi1;
        phi[2] = (1.0 / 2.0) * (3.0 * xsi2 - 1.0) * (3.0 * xsi2 - 2.0) * xsi2;
        phi[3] = (1.0 / 2.0) * (3.0 * xsi3 - 1.0) * (3.0 * xsi3 - 2.0) * xsi3;
        phi[4] = (9.0 / 2.0) * (3.0 * xsi4 - 1.0) * xsi1 * xsi4;
        phi[5] = (9.0 / 2.0) * (3.0 * xsi1 - 1.0) * xsi1 * xsi4;
        phi[6] = (9.0 / 2.0) * (3.0 * xsi1 - 1.0) * xsi1 * xsi2;
        phi[7] = (9.0 / 2.0) * (3.0 * xsi2 - 1.0) * xsi1 * xsi2;
        phi[8] = (9.0 / 2.0) * (3.0 * xsi2 - 1.0) * xsi2 * xsi4;
        phi[9] = (9.0 / 2.0) * (3.0 * xsi4 - 1.0) * xsi2 * xsi4;
        phi[10] = (9.0 / 2.0) * (3.0 * xsi3 - 1.0) * xsi3 * xsi4;
        phi[11] = (9.0 / 2.0) * (3.0 * xsi4 - 1.0) * xsi3 * xsi4;
        phi[12] = (9.0 / 2.0) * (3.0 * xsi3 - 1.0) * xsi2 * xsi3;
        phi[13] = (9.0 / 2.0) * (3.0 * xsi2 - 1.0) * xsi2 * xsi3;
        phi[14] = (9.0 / 2.0) * (3.0 * xsi3 - 1.0) * xsi1 * xsi3;
        phi[15] = (9.0 / 2.0) * (3.0 * xsi1 - 1.0) * xsi1 * xsi3;
        phi[16] = 27.0 * xsi1 * xsi2 * xsi4;
        phi[17] = 27.0 * xsi1 * xsi3 * xsi4;
        phi[18] = 27.0 * xsi2 * xsi3 * xsi4;
        phi[19] = 27.0 * xsi1 * xsi2 * xsi3;
        break;
    }
}

void ParametricVolumeElement::getShapeFunctionsDerivatives(double* xsi, double**& dphi_dxsi) const
{
    dphi_dxsi = new double*[3];
    for (int i = 0; i < 3; i++)
    {
        dphi_dxsi[i] = new double[numberOfNodes_];
    }
    double xsi1 = xsi[0];
    double xsi2 = xsi[1];
    double xsi3 = xsi[2];

    switch (elementType_)
    {
    case PartitionOfUnity::TET4:
        dphi_dxsi[0][0] = -1.0;
        dphi_dxsi[0][1] = 1.0;
        dphi_dxsi[0][2] = 0.0;
        dphi_dxsi[0][3] = 0.0;

        dphi_dxsi[1][0] = -1.0;
        dphi_dxsi[1][1] = 0.0;
        dphi_dxsi[1][2] = 1.0;
        dphi_dxsi[1][3] = 0.0;

        dphi_dxsi[2][0] = -1.0;
        dphi_dxsi[2][1] = 0.0;
        dphi_dxsi[2][2] = 0.0;
        dphi_dxsi[2][3] = 1.0;
        break;
    case PartitionOfUnity::TET10:
        dphi_dxsi[0][0] = -3.0 + 4.0 * xsi1 + 4.0 * xsi2 + 4.0 * xsi3;
        dphi_dxsi[0][1] = -1.0 + 4.0 * xsi1;
        dphi_dxsi[0][2] = 0.0;
        dphi_dxsi[0][3] = 0.0;
        dphi_dxsi[0][4] = 4.0 - 8.0 * xsi1 - 4.0 * xsi2 - 4.0 * xsi3;
        dphi_dxsi[0][5] = 4.0 * xsi2;
        dphi_dxsi[0][6] = -4.0 * xsi2;
        dphi_dxsi[0][7] = -4.0 * xsi3;
        dphi_dxsi[0][8] = 0.0;
        dphi_dxsi[0][9] = 4.0 * xsi3;

        dphi_dxsi[1][0] = -3.0 + 4.0 * xsi1 + 4.0 * xsi2 + 4.0 * xsi3;
        dphi_dxsi[1][1] = 0.0;
        dphi_dxsi[1][2] = -1.0 + 4.0 * xsi2;
        dphi_dxsi[1][3] = 0.0;
        dphi_dxsi[1][4] = -4.0 * xsi1;
        dphi_dxsi[1][5] = 4.0 * xsi1;
        dphi_dxsi[1][6] = 4.0 - 4.0 * xsi1 - 8.0 * xsi2 - 4.0 * xsi3;
        dphi_dxsi[1][7] = -4.0 * xsi3;
        dphi_dxsi[1][8] = 4.0 * xsi3;
        dphi_dxsi[1][9] = 0.0;

        dphi_dxsi[2][0] = -3.0 + 4.0 * xsi1 + 4.0 * xsi2 + 4.0 * xsi3;
        dphi_dxsi[2][1] = 0.0;
        dphi_dxsi[2][2] = 0.0;
        dphi_dxsi[2][3] = -1.0 + 4.0 * xsi3;
        dphi_dxsi[2][4] = -4.0 * xsi1;
        dphi_dxsi[2][5] = 0.0;
        dphi_dxsi[2][6] = -4.0 * xsi2;
        dphi_dxsi[2][7] = 4.0 - 4.0 * xsi1 - 4.0 * xsi2 - 8.0 * xsi3;
        dphi_dxsi[2][8] = 4.0 * xsi2;
        dphi_dxsi[2][9] = 4.0 * xsi1;
        break;
    case PartitionOfUnity::TET20:
        dphi_dxsi[0][0] = 0.50 * (-11.0 + 36.0 * xsi1 - 27.0 * xsi1 * xsi1 + 36.0 * xsi2 - 54.0 * xsi1 * xsi2 - 27.0 * xsi2 * xsi2 + 36.0 * xsi3 - 54.0 * xsi1 * xsi3 - 54.0 * xsi2 * xsi3 - 27.0 * xsi3 * xsi3);
        dphi_dxsi[0][1] = 0.50 * (2.0 - 18.0 * xsi1 + 27.0 * xsi1 * xsi1);
        dphi_dxsi[0][2] = 0.0;
        dphi_dxsi[0][3] = 0.0;
        dphi_dxsi[0][4] = 9.0 / 2.0 * (2.0 - 10.0 * xsi1 + 9.0 * xsi1 * xsi1 - 5.0 * xsi2 + 12.0 * xsi1 * xsi2 + 3.0 * xsi2 * xsi2 - 5.0 * xsi3 + 12.0 * xsi1 * xsi3 + 6.0 * xsi2 * xsi3 + 3.0 * xsi3 * xsi3);
        dphi_dxsi[0][5] = -9.0 / 2.0 * (1.0 - 8.0 * xsi1 + 9.0 * xsi1 * xsi1 - xsi2 + 6.0 * xsi1 * xsi2 - xsi3 + 6.0 * xsi1 * xsi3);
        dphi_dxsi[0][6] = 9.0 / 2.0 * (-1.0 + 6.0 * xsi1) * xsi2;
        dphi_dxsi[0][7] = 9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
        dphi_dxsi[0][8] = -9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
        dphi_dxsi[0][9] = 9.0 / 2.0 * xsi2 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
        dphi_dxsi[0][10] = -9.0 / 2.0 * xsi3 * (-1.0 + 3.0 * xsi3);
        dphi_dxsi[0][11] = 9.0 / 2.0 * xsi3 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
        dphi_dxsi[0][12] = 0.0;
        dphi_dxsi[0][13] = 0.0;
        dphi_dxsi[0][14] = 9.0 / 2.0 * xsi3 * (-1.0 + 3.0 * xsi3);
        dphi_dxsi[0][15] = 9.0 / 2.0 * (-1.0 + 6.0 * xsi1) * xsi3;
        dphi_dxsi[0][16] = -27.0 * xsi2 * (-1.0 + 2.0 * xsi1 + xsi2 + xsi3);
        dphi_dxsi[0][17] = -27.0 * xsi3 * (-1.0 + 2.0 * xsi1 + xsi2 + xsi3);
        dphi_dxsi[0][18] = -27.0 * xsi2 * xsi3;
        dphi_dxsi[0][19] = 27.0 * xsi2 * xsi3;

        dphi_dxsi[1][0] = 0.50 * (-11.0 + 36.0 * xsi1 - 27.0 * xsi1 * xsi1 + 36.0 * xsi2 - 54.0 * xsi1 * xsi2 - 27.0 * xsi2 * xsi2 + 36.0 * xsi3 - 54.0 * xsi1 * xsi3 - 54.0 * xsi2 * xsi3 - 27.0 * xsi3 * xsi3);
        dphi_dxsi[1][1] = 0.0;
        dphi_dxsi[1][2] = 0.50 * (2.0 - 18.0 * xsi2 + 27.0 * xsi2 * xsi2);
        dphi_dxsi[1][3] = 0.0;
        dphi_dxsi[1][4] = 9.0 / 2.0 * xsi1 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
        dphi_dxsi[1][5] = -9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
        dphi_dxsi[1][6] = 9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
        dphi_dxsi[1][7] = 9.0 / 2.0 * xsi1 * (-1.0 + 6.0 * xsi2);
        dphi_dxsi[1][8] = -9.0 / 2.0 * (1.0 - xsi1 - 8.0 * xsi2 - xsi3 + 6.0 * xsi1 * xsi2 + 6.0 * xsi2 * xsi3 + 9.0 * xsi2 * xsi2);
        dphi_dxsi[1][9] = 9.0 / 2.0 * (2.0 - 5.0 * xsi1 + 3.0 * xsi1 * xsi1 - 10.0 * xsi2 + 12.0 * xsi1 * xsi2 + 9.0 * xsi2 * xsi2 - 5.0 * xsi3 + 6.0 * xsi1 * xsi3 + 12.0 * xsi2 * xsi3 + 3.0 * xsi3 * xsi3);
        dphi_dxsi[1][10] = -9.0 / 2.0 * xsi3 * (-1.0 + 3.0 * xsi3);
        dphi_dxsi[1][11] = 9.0 / 2.0 * xsi3 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
        dphi_dxsi[1][12] = 9.0 / 2.0 * (-1.0 + 3.0 * xsi3) * xsi3;
        dphi_dxsi[1][13] = 9.0 / 2.0 * (-1.0 + 6.0 * xsi2) * xsi3;
        dphi_dxsi[1][14] = 0.0;
        dphi_dxsi[1][15] = 0.0;
        dphi_dxsi[1][16] = -27.0 * xsi1 * (-1.0 + xsi1 + 2.0 * xsi2 + xsi3);
        dphi_dxsi[1][17] = -27.0 * xsi1 * xsi3;
        dphi_dxsi[1][18] = -27.0 * xsi3 * (-1.0 + xsi1 + 2.0 * xsi2 + xsi3);
        dphi_dxsi[1][19] = 27.0 * xsi1 * xsi3;

        dphi_dxsi[2][0] = 0.50 * (-11.0 + 36.0 * xsi1 - 27.0 * xsi1 * xsi1 + 36.0 * xsi2 - 54.0 * xsi1 * xsi2 - 27.0 * xsi2 * xsi2 + 36.0 * xsi3 - 54.0 * xsi1 * xsi3 - 54.0 * xsi2 * xsi3 - 27.0 * xsi3 * xsi3);
        dphi_dxsi[2][1] = 0.0;
        dphi_dxsi[2][2] = 0.0;
        dphi_dxsi[2][3] = 0.50 * (2.0 - 18.0 * xsi3 + 27.0 * xsi3 * xsi3);
        dphi_dxsi[2][4] = 9.0 / 2.0 * xsi1 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
        dphi_dxsi[2][5] = -9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
        dphi_dxsi[2][6] = 0.0;
        dphi_dxsi[2][7] = 0.0;
        dphi_dxsi[2][8] = -9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
        dphi_dxsi[2][9] = 9.0 / 2.0 * xsi2 * (-5.0 + 6.0 * xsi1 + 6.0 * xsi2 + 6.0 * xsi3);
        dphi_dxsi[2][10] = -9.0 / 2.0 * (1.0 - xsi1 - xsi2 - 8.0 * xsi3 + 6.0 * xsi1 * xsi3 + 6.0 * xsi2 * xsi3 + 9.0 * xsi3 * xsi3);
        dphi_dxsi[2][11] = 9.0 / 2.0 * (2.0 - 5.0 * xsi1 + 3.0 * xsi1 * xsi1 - 5.0 * xsi2 + 6.0 * xsi1 * xsi2 + 3.0 * xsi2 * xsi2 - 10.0 * xsi3 + 12.0 * xsi1 * xsi3 + 12.0 * xsi2 * xsi3 + 9.0 * xsi3 * xsi3);
        dphi_dxsi[2][12] = 9.0 / 2.0 * xsi2 * (-1.0 + 6.0 * xsi3);
        dphi_dxsi[2][13] = 9.0 / 2.0 * xsi2 * (-1.0 + 3.0 * xsi2);
        dphi_dxsi[2][14] = 9.0 / 2.0 * xsi1 * (-1.0 + 6.0 * xsi3);
        dphi_dxsi[2][15] = 9.0 / 2.0 * xsi1 * (-1.0 + 3.0 * xsi1);
        dphi_dxsi[2][16] = -27.0 * xsi1 * xsi2;
        dphi_dxsi[2][17] = -27.0 * xsi1 * (-1.0 + xsi1 + 2.0 * xsi3 + xsi2);
        dphi_dxsi[2][18] = -27.0 * xsi2 * (-1.0 + xsi1 + 2.0 * xsi3 + xsi2);
        dphi_dxsi[2][19] = 27.0 * xsi1 * xsi2;
        break;
    }
}