#include "GeometricNeumann.h"

GeometricNeumann::GeometricNeumann(const int index, const double valueX, const double valueY, const double valueZ)
    : index_(index), valueX_(valueX), valueY_(valueY), valueZ_(valueZ) {}

GeometricNeumann::~GeometricNeumann() {}

int GeometricNeumann::getIndex() const
{
    return index_;
}

double GeometricNeumann::getValueX() const
{
    return valueX_;
}
        
double GeometricNeumann::getValueY() const
{
    return valueY_;
}

double GeometricNeumann::getValueZ() const
{
    return valueZ_;
}

GeometricPointNeumann::GeometricPointNeumann(const int index, Point* const point, const double valueX, const double valueY, const double valueZ)
    : GeometricNeumann(index, valueX, valueY, valueZ), point_(point) {}

GeometricPointNeumann::~GeometricPointNeumann() {}

Point* GeometricPointNeumann::getPoint() const
{
    return point_;
}

Line* GeometricPointNeumann::getLine() const
{
    return nullptr;
}

Surface* GeometricPointNeumann::getSurface() const
{
    return nullptr;
}

GeometricLineNeumann::GeometricLineNeumann(const int index, Line* const line, const double valueX, const double valueY, const double valueZ)
    : GeometricNeumann(index, valueX, valueY, valueZ), line_(line) {}

GeometricLineNeumann::~GeometricLineNeumann() {}

Point* GeometricLineNeumann::getPoint() const
{
    return nullptr;
}

Line* GeometricLineNeumann::getLine() const
{
    return line_;
}

Surface* GeometricLineNeumann::getSurface() const
{
    return nullptr;
}

GeometricSurfaceNeumann::GeometricSurfaceNeumann(const int index, Surface* const surface, const double valueX, const double valueY, const double valueZ)
    : GeometricNeumann(index, valueX, valueY, valueZ), surface_(surface) {}

GeometricSurfaceNeumann::~GeometricSurfaceNeumann() {}

Point* GeometricSurfaceNeumann::getPoint() const
{
    return nullptr;
}

Line* GeometricSurfaceNeumann::getLine() const
{
    return nullptr;
}

Surface* GeometricSurfaceNeumann::getSurface() const
{
    return surface_;
}