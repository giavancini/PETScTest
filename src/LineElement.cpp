#include "LineElement.h"

LineElement::LineElement(const int& index,
                         const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                         Material* material,
                         BaseLineElement* base,
                         AnalysisParameters* params)
        : Element(index, degreesOfFreedons),
          material_(material),
          base_(base),
          parameters_(params) {}

LineElement::~LineElement()
{
    delete base_;
}

void LineElement::setMaterial(Material* material)
{
    material_ = material;
}

void LineElement::setBaseElement(BaseLineElement* base)
{
    base_ = base;
}

void LineElement::setAnalysisParameters(AnalysisParameters* params)
{
    parameters_ = params;
}
        
Material* LineElement::getMaterial() const
{
    return material_;
}

BaseLineElement* LineElement::getBaseElement() const
{
    return base_;
}

ParametricElement* LineElement::getParametricElement() const
{
    return base_->getParametricElement();
}

const std::vector<Node*>& LineElement::getNodes() const
{
    return base_->getNodes();
}

void LineElement::getDOFIndexes(unsigned int& ndof,
                                 int*& indexes) const
{
    ndof = degreesOfFreedom_.size();
    indexes = new int[ndof];

    for (int i = 0; i < ndof; i++)
    {   
        indexes[i] = degreesOfFreedom_[i]->getIndex();
    }
}

void LineElement::elementContributions(int& ndofs1,
                                        int& ndofs2,
                                        int*& indexes,
                                        double*& rhsValues,
                                        double*& hessianValues) const {}

void LineElement::clearNeighborElements() {}

void LineElement::getCauchyStress(double**& nodalCauchyStress) const {}