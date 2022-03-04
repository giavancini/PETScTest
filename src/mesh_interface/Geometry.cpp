#include "Geometry.h"

Geometry::Geometry() {}

Geometry::Geometry(const int& index)
{
	index_ = index;
}

Geometry::~Geometry() {}

int Geometry::getIndex()
{
	return index_;
}

int Geometry::getNumberOfPoints()
{
	return points_.size();
}

int Geometry::getNumberOfLines()
{
	return lines_.size();
}

int Geometry::getNumberOfLineLoops()
{
	return lineLoops_.size();
}

int Geometry::getNumberOfSurfaces()
{
	return surfaces_.size();
}

int Geometry::getNumberOfSurfaceLoops()
{
	return surfaceLoops_.size();
}

int Geometry::getNumberOfVolumes()
{
	return volumes_.size();
}

Point* Geometry::getPoint(const std::string& name)
{
	return points_[name];
}

Line* Geometry::getLine(const std::string& name)
{
	return lines_[name];
}

LineLoop* Geometry::getLineLoop(const std::string& name)
{
	return lineLoops_[name];
}

Surface* Geometry::getSurface(const std::string& name)
{
	return surfaces_[name];
}

SurfaceLoop* Geometry::getSurfaceLoop(const std::string& name)
{
	return surfaceLoops_[name];
}

Volume* Geometry::getVolume(const std::string& name)
{
	return volumes_[name];
}

const std::unordered_map<std::string, Line*>& Geometry::getLines()
{
	return lines_;
}

const std::unordered_map<std::string, Surface*>& Geometry::getSurfaces()
{
	return surfaces_;
}

const std::unordered_map<std::string, Volume*>& Geometry::getVolumes()
{
	return volumes_;
}

const std::vector<GeometricDirichlet*>& Geometry::getDirichletBoundaryConditions()
{
	return dirichlet_;
}

const std::vector<GeometricNeumann*>& Geometry::getNeumannBoundaryConditions()
{
	return neumann_;
}

const std::vector<GeometricInterface*>& Geometry::getInterfaceBoundaryConditions()
{
	return interface_;
}

std::string Geometry::getGmshCode()
{
	return gmshCode_;
}

Point* Geometry::addPoint(std::vector<double> coordinates, const double& lcar, const bool& discretization)
{
	int index = getNumberOfPoints();
	std::stringstream name;
	name << "p" << index;
	Point* p = new Point(index, name.str(), coordinates, lcar, discretization);
	points_[p->getName()] = p;
	gmshCode_ += p->getGmshCode();
	return p;
}

Line* Geometry::addLine(std::vector<Point*> points, const bool& discretization)
{
	int index = getNumberOfLines();
	std::stringstream name;
	name << "l" << index;
	Line* l = new Line(index, name.str(), points, discretization);
	lines_[l->getName()] = l;
	gmshCode_ += l->getGmshCode();
	return l;
}

Circle* Geometry::addCircle(std::vector<Point*> points, const bool& discretization)
{
	int index = getNumberOfLines();
	std::stringstream name;
	name << "l" << index;
	Circle* l = new Circle(index, name.str(), points, discretization);
	lines_[l->getName()] = l;
	gmshCode_ += l->getGmshCode();
	return l;
}

Spline* Geometry::addSpline(std::vector<Point*> edgePoints, std::function<double(double)> function, const int& ndiv, const bool& discretization)
{
	int index = getNumberOfLines();
	std::stringstream name;
	name << "l" << index;
	double dx = (edgePoints[1]->getX() - edgePoints[0]->getX()) / (ndiv - 1);
	double dy = (edgePoints[1]->getY() - edgePoints[0]->getY()) / (ndiv - 1);
	double x = edgePoints[0]->getX();
	double z = edgePoints[0]->getZ();
	std::vector<Point*> points; points.reserve(ndiv);
	points.push_back(edgePoints[0]);
	for (int i = 0; i < ndiv - 2; i++)
	{
		x += dx;
		double y = function(x);
		int index = getNumberOfPoints();
		std::stringstream name;
		name << "p" << index;
		Point* p = new Point(index, name.str(), {x, y, z}, 1.0, false);
		points_[p->getName()] = p;
		gmshCode_ += p->getGmshCode();
		points.push_back(p);
	}
	points.push_back(edgePoints[1]);

	Spline* l = new Spline(index, name.str(), points, discretization);
	lines_[l->getName()] = l;
	gmshCode_ += l->getGmshCode();
	return l;
}

LineLoop* Geometry::addLineLoop(std::vector<Line*> lines)
{
	int index = getNumberOfLineLoops();
	std::stringstream name;
	name << "ll" << index;
	LineLoop* ll = new LineLoop(index, name.str(), lines);
	ll->verification();
	lineLoops_[ll->getName()] = ll;
	gmshCode_ += ll->getGmshCode();
	return ll;
}

Surface* Geometry::addSurface(LineLoop* lineLoop)
{
	int index = getNumberOfSurfaces();
	std::stringstream name;
	name << "s" << index;
	Surface* s = new Surface(index, name.str(), lineLoop);
	surfaces_[s->getName()] = s;
	gmshCode_ += s->getGmshCode();
	return s;
}

Surface* Geometry::addSurface(std::vector<Line*> lines)
{
	int index = getNumberOfSurfaces();
	std::stringstream name;
	name << "s" << index;
	LineLoop* ll = addLineLoop(lines);
	Surface* s = new Surface(index, name.str(), ll);
	surfaces_[s->getName()] = s;
	gmshCode_ += s->getGmshCode();
	return s;
}

PlaneSurface* Geometry::addPlaneSurface(LineLoop* lineLoop)
{
	int index = getNumberOfSurfaces();
	std::stringstream name;
	name << "s" << index;
	PlaneSurface* s = new PlaneSurface(index, name.str(), lineLoop);
	surfaces_[s->getName()] = s;
	gmshCode_ += s->getGmshCode();
	return s;
}

PlaneSurface* Geometry::addPlaneSurface(std::vector<Line*> lines)
{
	int index = getNumberOfSurfaces();
	std::stringstream name;
	name << "s" << index;
	LineLoop* ll = addLineLoop(lines);
	PlaneSurface* s = new PlaneSurface(index, name.str(), ll);
	surfaces_[s->getName()] = s;
	gmshCode_ += s->getGmshCode();
	return s;
}

SurfaceLoop* Geometry::addSurfaceLoop(std::vector<Surface*> surfaces)
{
	int index = getNumberOfSurfaceLoops();
	std::stringstream name;
	name << "sl" << index;
	SurfaceLoop* sl = new SurfaceLoop(index, name.str(), surfaces);
	sl->verification();
	surfaceLoops_[sl->getName()] = sl;
	gmshCode_ += sl->getGmshCode();
	return sl;
}

Volume* Geometry::addVolume(SurfaceLoop* surfaceLoop)
{
	int index = getNumberOfVolumes();
	std::stringstream name;
	name << "v" << index;
	Volume* v = new Volume(index, name.str(), surfaceLoop);
	volumes_[v->getName()] = v;
	gmshCode_ += v->getGmshCode();
	return v;
}

Volume* Geometry::addVolume(std::vector<Surface*> surfaces)
{
	int index = getNumberOfVolumes();
	std::stringstream name;
	name << "v" << index;
	SurfaceLoop* sl = addSurfaceLoop(surfaces);
	Volume* v = new Volume(index, name.str(), sl);
	volumes_[v->getName()] = v;
	gmshCode_ += v->getGmshCode();
	return v;
}

void Geometry::appendGmshCode(std::string text)
{
	gmshCode_ += text;
}

void Geometry::transfiniteLine(std::vector<Line*> lines, const int& divisions, const double& progression)
{
	std::stringstream text;
	text << "Transfinite Line {";
	for (size_t i = 0; i < lines.size(); i++)
	{
		text << lines[i]->getName();
		if (i != (lines.size() - 1))
			text << ", ";
	}
	text << "} = " << divisions << " Using Progression " << progression << ";\n//\n";
	gmshCode_ += text.str();
}

void Geometry::transfiniteSurface(std::vector<Surface*> Surfaces, std::string orientation, std::vector<Point*> points)
{
	std::stringstream text;
	text << "Transfinite Surface {";
	for (size_t i = 0; i < Surfaces.size(); i++)
	{
		text << Surfaces[i]->getName();
		if (i != (Surfaces.size() - 1))
			text << ", ";
	}
	text << "} ";
	if (points.size() != 0)
	{
		text << "= {";
		for (size_t i = 0; i < points.size(); i++)
		{
			text << points[i]->getName();
			if (i != (points.size() - 1))
				text << ", ";
		}
		text << "} " << orientation << ";\n//\n";
	}
	else
	{
		text << orientation << ";\n//\n";
	}
	gmshCode_ += text.str();
}

void Geometry::transfiniteVolume(std::vector<Volume*> volumes)
{
	std::stringstream text;
	text << "Transfinite Volume {";
	for (size_t i = 0; i < volumes.size(); i++)
	{
		text << volumes[i]->getName();
		if (i != (volumes.size() - 1))
			text << ", ";
	}
	text << "};\n//\n";
	
	gmshCode_ += text.str();
}

void Geometry::addDirichletBoundaryCondition(const std::vector<Point*>& points, const ConstrainedDOF dof, const double value)
{
	int index = dirichlet_.size() - 1;
    for (Point* point : points)
        dirichlet_.emplace_back(new GeometricPointDirichlet(++index, point, dof, value));
}
	
void Geometry::addDirichletBoundaryCondition(const std::vector<Line*>& lines, const ConstrainedDOF dof, const double value)
{
	int index = dirichlet_.size() - 1;
    for (Line* line : lines)
        dirichlet_.emplace_back(new GeometricLineDirichlet(++index, line, dof, value));
}

void Geometry::addDirichletBoundaryCondition(const std::vector<Surface*>& surfaces, const ConstrainedDOF dof, const double value)
{
	int index = dirichlet_.size() - 1;
    for (Surface* surface : surfaces)
        dirichlet_.emplace_back(new GeometricSurfaceDirichlet(++index, surface, dof, value));
}

void Geometry::addDirichletBoundaryCondition(const std::vector<Volume*>& volumes, const ConstrainedDOF dof, const double value)
{
	int index = dirichlet_.size() - 1;
    for (Volume* volume : volumes)
        dirichlet_.emplace_back(new GeometricVolumeDirichlet(++index, volume, dof, value));
}

void Geometry::addNeumannBoundaryCondition(const std::vector<Point*>& points, const double valueX, const double valueY, const double valueZ)
{
	int index = neumann_.size() - 1;
    for (Point* point : points)
        neumann_.emplace_back(new GeometricPointNeumann(++index, point, valueX, valueY, valueZ));
}
	
void Geometry::addNeumannBoundaryCondition(const std::vector<Line*>& lines, const double valueX, const double valueY, const double valueZ)
{
	int index = neumann_.size() - 1;
    for (Line* line : lines)
        neumann_.emplace_back(new GeometricLineNeumann(++index, line, valueX, valueY, valueZ));
}

void Geometry::addNeumannBoundaryCondition(const std::vector<Surface*>& surfaces, const double valueX, const double valueY, const double valueZ)
{
	int index = neumann_.size() - 1;
    for (Surface* surface : surfaces)
        neumann_.emplace_back(new GeometricSurfaceNeumann(++index, surface, valueX, valueY, valueZ));
}

void Geometry::addInterfaceBoundaryCondition(const std::vector<Point*>& points)
{
	int index = interface_.size() - 1;
    for (Point* point : points)
        interface_.emplace_back(new GeometricPointInterface(++index, point));
}
	
void Geometry::addInterfaceBoundaryCondition(const std::vector<Line*>& lines)
{
	int index = interface_.size() - 1;
    for (Line* line : lines)
        interface_.emplace_back(new GeometricLineInterface(++index, line));
}

void Geometry::addInterfaceBoundaryCondition(const std::vector<Surface*>& surfaces)
{
	int index = interface_.size() - 1;
    for (Surface* surface : surfaces)
        interface_.emplace_back(new GeometricSurfaceInterface(++index, surface));
}
