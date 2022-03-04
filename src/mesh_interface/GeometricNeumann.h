#pragma once
#include "Point.h"
#include "Line.h"
#include "Surface.h"

class GeometricNeumann
{
    protected:
        int index_;
        double valueX_;
        double valueY_;
        double valueZ_;

    public:

        GeometricNeumann(const int index, const double valueX, const double valueY, const double valueZ);

        virtual ~GeometricNeumann() = 0;

        int getIndex() const;

        double getValueX() const;

        double getValueY() const;

        double getValueZ() const;

        virtual Point* getPoint() const = 0;

        virtual Line* getLine() const = 0;

        virtual Surface* getSurface() const = 0;
};

class GeometricPointNeumann : public GeometricNeumann
{
    private:
        Point* point_;

    public:

        GeometricPointNeumann(const int index, Point* const point, const double valueX, const double valueY, const double valueZ);

        ~GeometricPointNeumann() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;
};

class GeometricLineNeumann : public GeometricNeumann
{
    private:
        Line* line_;

    public:

        GeometricLineNeumann(const int index, Line* const line, const double valueX, const double valueY, const double valueZ);

        ~GeometricLineNeumann() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;
};

class GeometricSurfaceNeumann : public GeometricNeumann
{
    private:
        Surface* surface_;
        
    public:

        GeometricSurfaceNeumann(const int index, Surface* const surface, const double valueX, const double valueY, const double valueZ);

        ~GeometricSurfaceNeumann() override;

        Point* getPoint() const override;

        Line* getLine() const override;

        Surface* getSurface() const override;
};