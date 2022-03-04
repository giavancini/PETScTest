#pragma once
#include "BaseLineElement.h"
#include "BaseSurfaceElement.h"
#include "BaseVolumeElement.h"

enum class ReferenceConfiguration
{
    INITIAL,
    PAST,
    CURRENT
};

class Element
{
    public:
    Element(const int index,
            const std::vector<DegreeOfFreedom*>& degreesOfFreedom);
        
    virtual ~Element() = default;

    void setIndex(const int& index);

    void setRank(const int& rank);

    void setBoundary(const bool& isBoundary);

    void setActive(const bool& isActive);

    void setReferenceConfiguration(const ReferenceConfiguration reference);

    void setDegreesOfFreedom(const std::vector<DegreeOfFreedom*>& degreesOfFreedom);

    int getIndex() const;

    int getRank() const;
    
    bool isBoundary() const;

    bool isActive() const;

    ReferenceConfiguration getReferenceConfiguration() const;

    const std::vector<DegreeOfFreedom*>& getDegreesOfFreedom() const;

    const unsigned int getNumberOfDOFs() const;

    virtual const std::vector<Node*>& getNodes() const = 0;

    virtual BaseElement* getBaseElement() const = 0;

    virtual ParametricElement* getParametricElement() const = 0;

    virtual void getDOFIndexes(unsigned int& ndof,
                               int*& indexes) const = 0;

    virtual void getCauchyStress(double**& nodalCauchyStress) const = 0;

    virtual void elementContributions(int& ndofs1,
                                        int& ndofs2,
                                        int*& indexes,
                                        double*& rhsValues,
                                        double*& hessianValues) const = 0;

    virtual void clearNeighborElements() = 0;

    void addNeighborElement(Element* el);

    const std::vector<Element*>& getNeighborElements() const;

    protected:
    int index_;
    int rank_;
    bool isBoundary_;
    bool isActive_;
    ReferenceConfiguration referenceConfiguration_;
    std::vector<DegreeOfFreedom*> degreesOfFreedom_;
    std::vector<Element*> neighborElements_;
};