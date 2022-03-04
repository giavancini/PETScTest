#pragma once

#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include "Node.h"
#include "PlaneElement.h"
#include "VolumeElement.h"
#include "AnalysisParameters.h"

class Mesher
{
public:
    Mesher();

    virtual ~Mesher();
    
    virtual void executePreMeshingProcesses(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    virtual void executePostMeshingProcesses(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    virtual void execute(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param) = 0;

    void setMeshInfo(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    void removeMeshNodes(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    void generateNewNodes(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    void selectMeshElements(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    void generateNewElements(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    void buildFreeSurface(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    void deactivateSlivers(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param);

    virtual bool alphaShape(std::vector<Node*>& nodes, const double& alpha, const double& meanMeshLength) = 0;

    struct MeshContainer
    {
    protected:
        double *pointList_;
        int *elementList_;

        double *elementSizeList_;
        int *elementNeighbourList_;

        int numberOfPoints_;
        int numberOfElements_;

    public:
        //flags to set when the pointers are created (true) or deleted (false)
        bool pointListFlag_;
        bool elementListFlag_;
        bool elementSizeListFlag_;
        bool elementNeighbourListFlag_;

        void setPointList(double*& pointList) { pointList_ = pointList; }
        void setElementList(int*& elementList) { elementList_ = elementList; };
        void setElementSizeList(double*& elementSizeList) { elementSizeList_ = elementSizeList; };
        void setElementNeighbourList(int*& elementNeighbourList) { elementNeighbourList_ = elementNeighbourList; };
        void setNumberOfPoints(int& numberOfPoints) { numberOfPoints_ = numberOfPoints; };
        void setNumberOfElements(int& numberOfElements) { numberOfElements_ = numberOfElements; };

        double* getPointList() { return pointList_; };
        int* getElementList() { return elementList_; };
        double* getElementSizeList() { return elementSizeList_; };
        int* getElementNeighbourList() { return elementNeighbourList_; };

        int& getNumberOfPoints() { return numberOfPoints_; };
        int& getNumberOfElements() { return numberOfElements_; };

        void createPointList(const unsigned int numberOfPoints,
                             const unsigned int dimension)
        {
            if (pointList_)
            {
                delete[] pointList_;
            }
            numberOfPoints_ = numberOfPoints;
            pointList_ = new double[numberOfPoints * dimension];
            pointListFlag_ = true;
        }

        void createElementList(const unsigned int numberOfElements,
                               const unsigned int numberOfVertices)
        {
            if (elementList_)
            {
                delete[] elementList_;
            }
            numberOfElements_ = numberOfElements;
            elementList_ = new int[numberOfElements * numberOfVertices];
            elementListFlag_ = true;
        }

        void createElementSizeList(const unsigned int numberOfElements)
        {
            if (elementSizeList_)
            {
                delete[] elementSizeList_;
            }
            elementSizeList_ = new double[numberOfElements];
            elementSizeListFlag_ = true;
        }

        void createElementNeighbourList(const unsigned int numberOfElements,
                                        const unsigned int numberOfFaces)
        {
            if (elementNeighbourList_)
            {
                delete[] elementNeighbourList_;
            }
            elementNeighbourList_ = new int[numberOfElements * numberOfFaces];
            elementNeighbourListFlag_ = true;
        }

        void initialize()
        {
            pointList_ = (double *)NULL;
            elementList_ = (int *)NULL;
            elementSizeList_ = (double *)NULL;
            elementNeighbourList_ = (int *)NULL;
            numberOfPoints_ = 0;
            numberOfElements_ = 0;

            pointListFlag_ = false;
            elementListFlag_ = false;
            elementSizeListFlag_ = false;
            elementNeighbourListFlag_ = false;
        }

        void finalize()
        {
            if (pointList_ && pointListFlag_)
            {
                delete[] pointList_;
            }

            if (elementList_ && elementListFlag_)
            {
                delete[] elementList_;
            }

            if (elementSizeList_ && elementSizeListFlag_)
            {
                delete[] elementSizeList_;
            }

            if (elementNeighbourList_ && elementNeighbourListFlag_)
            {
                delete[] elementNeighbourList_;
            }

            initialize();
        }
    };

    struct MeshingInfo
    {
        public:
        
        unsigned int numberOfElements_;
        unsigned int numberOfNodes_;
        unsigned int numberOfInitialNodes_;

        unsigned int numberOfNewElements_;
        unsigned int numberOfNewNodes_;

        unsigned int insertedNodes_;
        unsigned int removedNodes_;

        double initialMeshVolume_;

        void initialize()
        {
            numberOfElements_ = 0;
            numberOfNodes_ = 0;

            numberOfNewElements_ = 0;
            numberOfNewNodes_ = 0;

            insertedNodes_ = 0;
            removedNodes_ = 0;
        }

        void setNumberOfNodes(unsigned int numberOfNodes)
        {
            numberOfNodes_ = numberOfNodes;
        }

        void setNumberOfElements(unsigned int numberOfElements)
        {
            numberOfElements_ = numberOfElements;
        }

        unsigned int getNumberOfNodes()
        {
            return numberOfNodes_;
        }

        unsigned int getNumberOfElements()
        {
            return numberOfElements_;
        }

        void setNumberOfNewNodes(unsigned int numberOfNodes)
        {
            numberOfNewNodes_ = numberOfNodes;
        }

        void setNumberOfNewElements(unsigned int numberOfElements)
        {
            numberOfNewElements_ = numberOfElements;
        }

        void setInitialMeshVolume(double volume)
        {
            initialMeshVolume_ = volume;
        }

        double getInitialMeshVolume()
        {
            return initialMeshVolume_;
        }
    };

protected:

    MeshContainer inMesh_;
    MeshContainer outMesh_;

    MeshingInfo info_;

    std::vector<int> preservedElements_;
    int numberOfPreservedElements_;
};