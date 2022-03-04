#pragma once

#include "Volume.h"
#include "GeometricDirichlet.h"
#include "GeometricNeumann.h"
#include "GeometricInterface.h"
#include <unordered_map>
#include <functional>

class Geometry
{
public:
	Geometry();

	Geometry(const int& index);

	~Geometry();

	int getIndex();

	int getNumberOfPoints();

	int getNumberOfLines();

	int getNumberOfLineLoops();

	int getNumberOfSurfaces();

	int getNumberOfSurfaceLoops();

	int getNumberOfVolumes();

	int getNumberOfBoundaryConditions(const std::string& type);

	Point* getPoint(const std::string& name);

	Line* getLine(const std::string& name);

	LineLoop* getLineLoop(const std::string& name);

	Surface* getSurface(const std::string& name);

	SurfaceLoop* getSurfaceLoop(const std::string& name);

	Volume* getVolume(const std::string& name);

	const std::unordered_map<std::string, Line*>& getLines();

	const std::unordered_map<std::string, Surface*>& getSurfaces();

	const std::unordered_map<std::string, Volume*>& getVolumes();

	const std::vector<GeometricDirichlet*>& getDirichletBoundaryConditions();

	const std::vector<GeometricNeumann*>& getNeumannBoundaryConditions();

	const std::vector<GeometricInterface*>& getInterfaceBoundaryConditions();

	std::string getGmshCode();

	Point* addPoint(std::vector<double> coordinates, const double& lcar = 1.0, const bool& discretization = true);

	Line* addLine(std::vector<Point*> points, const bool& discretization = true);

	Circle* addCircle(std::vector<Point*> points, const bool& discretization = true);

	Spline* addSpline(std::vector<Point*> points, std::function<double(double)> function, const int& ndiv, const bool& discretization = true);

	LineLoop* addLineLoop(std::vector<Line*> lines);

	Surface* addSurface(LineLoop* lineLoop);

	Surface* addSurface(std::vector<Line*> lines);

	PlaneSurface* addPlaneSurface(LineLoop* lineLoop);

	PlaneSurface* addPlaneSurface(std::vector<Line*> lines);

	SurfaceLoop* addSurfaceLoop(std::vector<Surface*> surfaces);

	Volume* addVolume(SurfaceLoop* surfaceLoop);

	Volume* addVolume(std::vector<Surface*> surfaces);

	void appendGmshCode(std::string text);

	void transfiniteLine(std::vector<Line*> lines, const int& divisions, const double& progression = 1);

	void transfiniteSurface(std::vector<Surface*> surfaces, std::string oientation = "Left", std::vector<Point*> points = std::vector<Point*>());

	void transfiniteVolume(std::vector<Volume*> volumes);

	void addDirichletBoundaryCondition(const std::vector<Point*>& points, const ConstrainedDOF dof, const double value);
	
	void addDirichletBoundaryCondition(const std::vector<Line*>& lines, const ConstrainedDOF dof, const double value);

	void addDirichletBoundaryCondition(const std::vector<Surface*>& surfaces, const ConstrainedDOF dof, const double value);

	void addDirichletBoundaryCondition(const std::vector<Volume*>& volumes, const ConstrainedDOF dof, const double value);

	void addNeumannBoundaryCondition(const std::vector<Point*>& points, const double valueX, const double valueY, const double valueZ);
	
	void addNeumannBoundaryCondition(const std::vector<Line*>& lines, const double valueX, const double valueY, const double valueZ);

	void addNeumannBoundaryCondition(const std::vector<Surface*>& surfaces, const double valueX, const double valueY, const double valueZ);

	void addNeumannBoundaryCondition(const std::vector<Line*>& lines, const std::vector<double>& valuesX, const std::vector<double>& valuesY, const std::vector<double>& valuesZ);

	void addNeumannBoundaryCondition(const std::vector<Surface*>& surfaces, const std::vector<double>& valuesX, const std::vector<double>& valuesY, const std::vector<double>& valuesZ);

	void addInterfaceBoundaryCondition(const std::vector<Point*>& points);
	
	void addInterfaceBoundaryCondition(const std::vector<Line*>& lines);

	void addInterfaceBoundaryCondition(const std::vector<Surface*>& surfaces);

private:
	int index_;
	std::unordered_map<std::string, Point*> points_;
	std::unordered_map<std::string, Line*> lines_;
	std::unordered_map<std::string, LineLoop*> lineLoops_;
	std::unordered_map<std::string, Surface*> surfaces_;
	std::unordered_map<std::string, SurfaceLoop*> surfaceLoops_;
	std::unordered_map<std::string, Volume*> volumes_;
	std::vector<GeometricDirichlet*> dirichlet_;
	std::vector<GeometricNeumann*> neumann_;
	std::vector<GeometricInterface*> interface_;
	std::string gmshCode_;

public:
	friend class FluidDomain;
	friend class SolidDomain;
};

