#pragma once
#include "Point.h"
#include "Line.h"
#include "Surface.h"
#include "Volume.h"

class GeometricInterface
{
    protected:
        int index_;

    public:

        GeometricInterface(const int index);

        virtual ~GeometricInterface() = 0;

        int getIndex() const;

        virtual Point* getPoint() const = 0;

        virtual Line* getLine() const = 0;

        virtual Surface* getSurface() const = 0;

        virtual Volume* getVolume() const = 0;
};

class GeometricPointInterface : public GeometricInterface
{
    private:
        Point* point_;

    public:

        GeometricPointInterface(const int index, Point* const point);

        ~GeometricPointInterface() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};

class GeometricLineInterface : public GeometricInterface
{
    private:
        Line* line_;

    public:

        GeometricLineInterface(const int index, Line* const line);

        ~GeometricLineInterface() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};

class GeometricSurfaceInterface : public GeometricInterface
{
    private:
        Surface* surface_;

    public:

        GeometricSurfaceInterface(const int index, Surface* const surface);

        ~GeometricSurfaceInterface() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};

class GeometricVolumeInterface : public GeometricInterface
{
    private:
        Volume* volume_;

    public:

        GeometricVolumeInterface(const int index, Volume* const volume);

        ~GeometricVolumeInterface() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};