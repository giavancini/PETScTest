#include "GeometricInterface.h"

GeometricInterface::GeometricInterface(const int index)
    : index_(index) {}

GeometricInterface::~GeometricInterface() {}

int GeometricInterface::getIndex() const
{
    return index_;
}

GeometricPointInterface::GeometricPointInterface(const int index, Point* const point)
    : GeometricInterface(index), point_(point) {}

GeometricPointInterface::~GeometricPointInterface() {}

Point* GeometricPointInterface::getPoint() const
{
    return point_;
}

Line* GeometricPointInterface::getLine() const
{
    return nullptr;
}

Surface* GeometricPointInterface::getSurface() const
{
    return nullptr;
}

Volume* GeometricPointInterface::getVolume() const
{
    return nullptr;
}

GeometricLineInterface::GeometricLineInterface(const int index, Line* const line)
    : GeometricInterface(index), line_(line) {}

GeometricLineInterface::~GeometricLineInterface() {}

Point* GeometricLineInterface::getPoint() const
{
    return nullptr;
}

Line* GeometricLineInterface::getLine() const
{
    return line_;
}

Surface* GeometricLineInterface::getSurface() const
{
    return nullptr;
}

Volume* GeometricLineInterface::getVolume() const
{
    return nullptr;
}

GeometricSurfaceInterface::GeometricSurfaceInterface(const int index, Surface* const surface)
    : GeometricInterface(index), surface_(surface) {}

GeometricSurfaceInterface::~GeometricSurfaceInterface() {}

Point* GeometricSurfaceInterface::getPoint() const
{
    return nullptr;
}

Line* GeometricSurfaceInterface::getLine() const
{
    return nullptr;
}

Surface* GeometricSurfaceInterface::getSurface() const
{
    return surface_;
}

Volume* GeometricSurfaceInterface::getVolume() const
{
    return nullptr;
}

GeometricVolumeInterface::GeometricVolumeInterface(const int index, Volume* const volume)
    : GeometricInterface(index), volume_(volume) {}

GeometricVolumeInterface::~GeometricVolumeInterface() {}

Point* GeometricVolumeInterface::getPoint() const
{
    return nullptr;
}

Line* GeometricVolumeInterface::getLine() const
{
    return nullptr;
}

Surface* GeometricVolumeInterface::getSurface() const
{
    return nullptr;
}

Volume* GeometricVolumeInterface::getVolume() const
{
    return volume_;
}