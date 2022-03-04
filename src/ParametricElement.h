#pragma once
#include "../include/vtkCellType.h"
#include "Quadratures.h"
#include "QuadraturePoint.h"
#include <vector>
#include <cmath>
#include <iostream>

enum PartitionOfUnity
{
    //point
    VERTEX,
    // line partition of unity
    L2, L3, L4,
    // triangular partition of unity
    T3, T6, T10,
    // quadrilateral partition of unity
    Q4, Q9, Q16,
    // tetrahedron partition of unity
    TET4, TET10, TET20,
    // hexahedron partition of unity
    HEX8, HEX27, HEX64
};

class ParametricElement
{
    public:
    ParametricElement(const PartitionOfUnity elementType);

    virtual ~ParametricElement() = default;

    void setOrder(const int order);

    void setNumberOfQuadraturePoints(const int numberOfQuadraturePoints);

    void setVTKCellType(const VTKCellType vtkCellType);

    void setQuadraturePoints(const std::vector<QuadraturePoint*>& quadraturePoints);

    void setVTKConnectivity(const std::vector<int>& vtkConnectivity);

    virtual void setNodalParametricCoordinates() = 0;

    PartitionOfUnity getElementType() const;

    int getOrder() const;

    int getNumberOfQuadraturePoints() const;

    int getNumberOfNodes() const;

    int getNumberOfFaces() const;

    int getNumberOfEdges() const;

    double** getNodalParametricCoordinates() const;

    VTKCellType getVTKCellType() const;

    const std::vector<QuadraturePoint*>& getQuadraturePoints() const;

    const std::vector<int>& getVTKConnectivity() const;

    const std::vector<int>& getFaceNodes(const int faceNumber) const;

    const std::vector<int>& getFaceVertices(const int faceNumber) const;

    const std::vector<int>& getEdgeNodes(const int edgeNumber) const;

    const std::vector<int>& getEdgeVertices(const int edgeNumber) const;

    virtual void getShapeFunctions(double* xsi, double*& phi) const = 0;

    virtual void getShapeFunctionsDerivatives(double* xsi, double**& dphi_dxsi) const = 0;

    protected:
    const PartitionOfUnity elementType_;
    int order_;
    int numberOfQuadraturePoints_;
    int numberOfNodes_;
    int numberOfFaces_;
    int numberOfEdges_;
    VTKCellType vtkCellType_;
    std::vector<int> vtkConnectivity_;
    double** nodalParametricCoordinates_;
    std::vector<QuadraturePoint*> quadraturePoints_;
    std::vector<std::vector<int>> faceNodes_;
    std::vector<std::vector<int>> faceVertices_;
    std::vector<std::vector<int>> edgeNodes_;
    std::vector<std::vector<int>> edgeVertices_;
};