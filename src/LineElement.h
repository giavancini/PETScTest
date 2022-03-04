#pragma once

#include "Element.h"
#include "Material.h"
#include "AnalysisParameters.h"

class LineElement : public Element
{
    public:
    LineElement(const int& index,
                const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                Material* material,
                BaseLineElement* base,
                AnalysisParameters* params);

    ~LineElement() override;

    void setMaterial(Material* material);

    void setBaseElement(BaseLineElement* base);

    void setAnalysisParameters(AnalysisParameters* params);

    Material* getMaterial() const;

    BaseLineElement* getBaseElement() const override;

    ParametricElement* getParametricElement() const override;

    const std::vector<Node*>& getNodes() const override;

    void getDOFIndexes(unsigned int& ndof,
                       int*& indexes) const override;

    void getCauchyStress(double**& nodalCauchyStress) const override;

    void elementContributions(int& ndofs1,
                                        int& ndofs2,
                                        int*& indexes,
                                        double*& rhsValues,
                                        double*& hessianValues) const override;
    
    void clearNeighborElements() override;

    private:
    Material* material_;
    BaseLineElement* base_;
    AnalysisParameters* parameters_;
};