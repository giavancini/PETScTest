#pragma once
#include "ParametricElement.h"

class ParametricSurfaceElement : public ParametricElement
{
    public:
    ParametricSurfaceElement(const PartitionOfUnity elementType);

    ~ParametricSurfaceElement();

    void setNodalParametricCoordinates() override;

    void getShapeFunctions(double* xsi, double*& phi) const override;

    void getShapeFunctionsDerivatives(double* xsi, double**& dphi_dxsi) const override;

    static ParametricSurfaceElement T3;

    static ParametricSurfaceElement T6;

    static ParametricSurfaceElement T10;

    static ParametricSurfaceElement Q4;

    static ParametricSurfaceElement Q9;

    static ParametricSurfaceElement Q16;
};