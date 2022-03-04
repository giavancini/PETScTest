#include "GeometricDirichlet.h"

GeometricDirichlet::GeometricDirichlet(const int index, const ConstrainedDOF dof, const double value)
    : index_(index), dof_(dof), value_(value) {}

GeometricDirichlet::~GeometricDirichlet() {}

int GeometricDirichlet::getIndex() const
{
    return index_;
}

ConstrainedDOF GeometricDirichlet::getDegreeOfFreedom() const
{
    return dof_;
}

double GeometricDirichlet::getValue() const
{
    return value_;
}

GeometricPointDirichlet::GeometricPointDirichlet(const int index, Point* const point, const ConstrainedDOF dof, const double value)
    : GeometricDirichlet(index, dof, value), point_(point) {}

GeometricPointDirichlet::~GeometricPointDirichlet() {}

Point* GeometricPointDirichlet::getPoint() const
{
    return point_;
}

Line* GeometricPointDirichlet::getLine() const
{
    return nullptr;
}

Surface* GeometricPointDirichlet::getSurface() const
{
    return nullptr;
}

Volume* GeometricPointDirichlet::getVolume() const
{
    return nullptr;
}

GeometricLineDirichlet::GeometricLineDirichlet(const int index, Line* const line, const ConstrainedDOF dof, const double value)
    : GeometricDirichlet(index, dof, value), line_(line) {}

GeometricLineDirichlet::~GeometricLineDirichlet() {}

Point* GeometricLineDirichlet::getPoint() const
{
    return nullptr;
}

Line* GeometricLineDirichlet::getLine() const
{
    return line_;
}

Surface* GeometricLineDirichlet::getSurface() const
{
    return nullptr;
}

Volume* GeometricLineDirichlet::getVolume() const
{
    return nullptr;
}

GeometricSurfaceDirichlet::GeometricSurfaceDirichlet(const int index, Surface* const surface, const ConstrainedDOF dof, const double value)
    : GeometricDirichlet(index, dof, value), surface_(surface) {}

GeometricSurfaceDirichlet::~GeometricSurfaceDirichlet() {}

Point* GeometricSurfaceDirichlet::getPoint() const
{
    return nullptr;
}

Line* GeometricSurfaceDirichlet::getLine() const
{
    return nullptr;
}

Surface* GeometricSurfaceDirichlet::getSurface() const
{
    return surface_;
}

Volume* GeometricSurfaceDirichlet::getVolume() const
{
    return nullptr;
}

GeometricVolumeDirichlet::GeometricVolumeDirichlet(const int index, Volume* const volume, const ConstrainedDOF dof, const double value)
    : GeometricDirichlet(index, dof, value), volume_(volume) {}

GeometricVolumeDirichlet::~GeometricVolumeDirichlet() {}

Point* GeometricVolumeDirichlet::getPoint() const
{
    return nullptr;
}

Line* GeometricVolumeDirichlet::getLine() const
{
    return nullptr;
}

Surface* GeometricVolumeDirichlet::getSurface() const
{
    return nullptr;
}

Volume* GeometricVolumeDirichlet::getVolume() const
{
    return volume_;
}