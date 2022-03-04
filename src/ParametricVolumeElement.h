#pragma once
#include "ParametricElement.h"

class ParametricVolumeElement : public ParametricElement
{
    public:
    ParametricVolumeElement(const PartitionOfUnity elementType);

    ~ParametricVolumeElement();

    void setNodalParametricCoordinates() override;

    void getShapeFunctions(double* xsi, double*& phi) const override;

    void getShapeFunctionsDerivatives(double* xsi, double**& dphi_dxsi) const override;

    static ParametricVolumeElement TET4;

    static ParametricVolumeElement TET10;

    static ParametricVolumeElement TET20;

    static ParametricVolumeElement HEX8;

    static ParametricVolumeElement HEX27;

    static ParametricVolumeElement HEX64;
};