#pragma once
#include "DegreeOfFreedom.h"
#include "Material.h"
#include <vector>
#include <limits>

class Element;

class Node
{

    public:
        Node(const unsigned int& index, 
             const std::vector<DegreeOfFreedom*>& degreesOfFreedom);

        ~Node();

        Node(const Node& node) = delete;

        Node operator=(const Node& node) = delete;

        bool operator==(const Node& node) const;

        bool operator<(const Node& node) const;

        void setIndex(const unsigned int& index);

        void setPermutedIndex(const unsigned int& index);

        void setRank(const unsigned int& rank);

        void setMaterial(Material* material);

        void setDegreesOfFreedom(const std::vector<DegreeOfFreedom*>& degreesOfFreedom);

        void setCauchyStress(const unsigned int& size, 
                             double*& cauchyStress);

        void clearCauchyStress(const unsigned int& size);

        void incrementCauchyStress(const unsigned int& size, 
                                   double*& cauchyStress);

        void setCloudArea(const double& cloudArea);

        void setMeshSize(const double& meshSize);

        void setBoundary(const bool& isBoundary);

        void setFreeSurface(const bool& isFreeSurface);

        void setConstrain(const bool& isConstrained);

        void setBlocked(const bool& isBlocked);

        void setIsolated(const bool& isIsolated);

        void setInterface(const bool& isInterface);

        void setPreviouslyFreeSurface(const bool& wasFreeSurface);

        void setPreviouslyIsolated(const bool& wasIsolated);

        void setToRemove(const bool& isToRemove);

        void setNewEntity(const bool& isNewEntity);

        void setInterfaceNode(Node* interfaceNode);

        unsigned int getIndex() const;

        unsigned int getPermutedIndex() const;

        unsigned int getRank() const;

        Material* getMaterial() const;

        const std::vector<DegreeOfFreedom*>& getDegreesOfFreedom() const;

        double* getCauchyStress() const;

        unsigned int getNumberOfDegreesOfFreedom() const;

        DegreeOfFreedom* getDegreeOfFreedom(const unsigned int& index) const;

        double getCloudArea() const;

        double getMeshSize() const;

        bool isBoundary() const;

        bool isFreeSurface() const;

        bool isConstrained() const;

        bool isBlocked() const;

        bool isIsolated() const;

        bool isInterface() const;

        bool wasIsolated() const;

        bool wasFreeSurface() const;

        bool isToRemove() const;

        bool isNewEntity() const;

        const std::vector<Node*>& getNeighborNodes() const;

        const std::vector<Element*>& getNeighborElements() const;

        Node* getInterfaceNode() const;

        void addDegreeOfFreedom(DegreeOfFreedom* dof);

        void addNeighborNode(Node* node);

        void addNeighborElement(Element* el);

        void clearNeighborNodes();

        void clearNeighborElements();

        double distanceToNode(const Node& node, 
                              const unsigned int& dimension);

        double squareDistanceToNode(const Node& node, 
                                    const unsigned int& dimension);

        void searchNodesInRadius(const double& radius, const int& dimension, int& numberOfFoundNodes, std::vector<Node*>& foundNodes, std::vector<double>& distances);

    private:
        unsigned int index_;
        unsigned int permutedIndex_;
        unsigned int rank_;
        double cloudArea_;
        double meshSize_;
        bool isBoundary_;
        bool isFreeSurface_;
        bool isConstrained_;
        bool isBlocked_;
        bool isIsolated_;
        bool isInterface_;
        bool wasFreeSurface_;
        bool wasIsolated_;
        bool isToRemove_;
        bool isNewEntity_;
        Material* material_;
        double* cauchyStress_;
        Node* interfaceNode_;
        std::vector<DegreeOfFreedom*> degreesOfFreedom_;
        std::vector<Node*> neighborNodes_;
        std::vector<Element*> neighborElements_;
};
