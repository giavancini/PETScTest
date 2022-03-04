#pragma once
#include "Point.h"
#include "Line.h"
#include "Surface.h"
#include "Volume.h"

enum ConstrainedDOF
{
    X,
    Y,
    Z,
    ALL
};

class GeometricDirichlet
{
    protected:
        int index_;
        ConstrainedDOF dof_;
        double value_;

    public:

        GeometricDirichlet(const int index, const ConstrainedDOF dof, const double value);

        virtual ~GeometricDirichlet() = 0;

        int getIndex() const;

        ConstrainedDOF getDegreeOfFreedom() const;

        double getValue() const;

        virtual Point* getPoint() const = 0;

        virtual Line* getLine() const = 0;

        virtual Surface* getSurface() const = 0;

        virtual Volume* getVolume() const = 0;
};

class GeometricPointDirichlet : public GeometricDirichlet
{
    private:
        Point* point_;

    public:

        GeometricPointDirichlet(const int index, Point* const point, const ConstrainedDOF dof, const double value);

        ~GeometricPointDirichlet() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};

class GeometricLineDirichlet : public GeometricDirichlet
{
    private:
        Line* line_;

    public:

        GeometricLineDirichlet(const int index, Line* const line, const ConstrainedDOF dof, const double value);

        ~GeometricLineDirichlet() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};

class GeometricSurfaceDirichlet : public GeometricDirichlet
{
    private:
        Surface* surface_;

    public:

        GeometricSurfaceDirichlet(const int index, Surface* const surface, const ConstrainedDOF dof, const double value);

        ~GeometricSurfaceDirichlet() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};

class GeometricVolumeDirichlet : public GeometricDirichlet
{
    private:
        Volume* volume_;

    public:

        GeometricVolumeDirichlet(const int index, Volume* const volume, const ConstrainedDOF dof, const double value);

        ~GeometricVolumeDirichlet() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;

        Volume* getVolume() const override;
};