#pragma once
#include "ParametricLineElement.h"
#include "ParametricSurfaceElement.h"
#include "ParametricVolumeElement.h"
#include "Node.h"

class BaseElement
{
    public:
    BaseElement(const int index,
                ParametricElement& parametricElement,     
                const std::vector<Node*>& nodes);

    virtual ~BaseElement() = default;

    int getIndex() const;

    int getNumberOfNodes() const;

    ParametricElement* getParametricElement() const;

    const std::vector<Node*>& getNodes() const;

    Node* getNode(const int index) const;

    bool getPlot() const;

    void setPlot(const bool plot);

    virtual double getRadius() const = 0;

    virtual double getJacobianIntegration() const = 0;
    
    void setIndex(const int index);
    
    protected:
    int index_;
    ParametricElement* parametricElement_;
    std::vector<Node*> nodes_;
    bool plot_;
};