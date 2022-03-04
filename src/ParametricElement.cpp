#include "ParametricElement.h"

ParametricElement::ParametricElement(const PartitionOfUnity elementType)
                   : elementType_(elementType) {}

void ParametricElement::setOrder(const int order)
{
    order_ = order;
}

void ParametricElement::setNumberOfQuadraturePoints(const int numberOfQuadraturePoints)
{
    numberOfQuadraturePoints_ = numberOfQuadraturePoints;
}

void ParametricElement::setVTKCellType(const VTKCellType vtkCellType)
{
    vtkCellType_ = vtkCellType;
}

void ParametricElement::setQuadraturePoints(const std::vector<QuadraturePoint*>& quadraturePoints)
{
    quadraturePoints_ = quadraturePoints;
}

void ParametricElement::setVTKConnectivity(const std::vector<int>& vtkConnectivity)
{
    vtkConnectivity_ = vtkConnectivity;
}

PartitionOfUnity ParametricElement::getElementType() const
{
    return elementType_;
}

int ParametricElement::getOrder() const
{
    return order_;
}

int ParametricElement::getNumberOfQuadraturePoints() const
{
    return numberOfQuadraturePoints_;
}

int ParametricElement::getNumberOfNodes() const
{
    return numberOfNodes_;
}

int ParametricElement::getNumberOfFaces() const
{
    return numberOfFaces_;
}

int ParametricElement::getNumberOfEdges() const
{
    return numberOfEdges_;
}

double** ParametricElement::getNodalParametricCoordinates() const
{
    return nodalParametricCoordinates_;
}

VTKCellType ParametricElement::getVTKCellType() const
{
    return vtkCellType_;
}

const std::vector<QuadraturePoint*>& ParametricElement::getQuadraturePoints() const
{
    return quadraturePoints_;
}

const std::vector<int>& ParametricElement::getVTKConnectivity() const
{
    return vtkConnectivity_;
}

const std::vector<int>& ParametricElement::getFaceNodes(const int faceNumber) const
{
    return faceNodes_[faceNumber];
}

const std::vector<int>& ParametricElement::getFaceVertices(const int faceNumber) const
{
    return faceVertices_[faceNumber];
}

const std::vector<int>& ParametricElement::getEdgeNodes(const int edgeNumber) const
{
    return edgeNodes_[edgeNumber];
}

const std::vector<int>& ParametricElement::getEdgeVertices(const int edgeNumber) const
{
    return edgeVertices_[edgeNumber];
}