#include "Point.h"

Point::Point() {}

Point::Point(const int index,
			 const std::string name,
			 const std::vector<double> coordinates,
			 const double lcar,
			 const bool discretization)
	: index_(index),
	  name_(name),
	  coordinates_(coordinates),
	  lcar_(lcar),
	  discretization_(discretization),
	  node_(nullptr) {}

Point::~Point() {}

int Point::getIndex()
{
	return index_;
}
std::string Point::getName()
{
	return name_;
}
double Point::getX()
{
	return coordinates_[0];
}
double Point::getY()
{
	return coordinates_[1];
}
double Point::getZ()
{
	return coordinates_[2];
}
double Point::getLcar()
{
	return lcar_;
}
bool Point::getDiscretization()
{
	return discretization_;
}
Node* Point::getNode()
{
	return node_;
}
std::string Point::getGmshCode()
{
	std::stringstream text;
	text << std::fixed;
	if (discretization_) {
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", " << coordinates_[2] << ", " << lcar_
			<< "}; Physical Point('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newp; Point(" << name_ << ") = {" << coordinates_[0] << ", " << coordinates_[1] << ", " << coordinates_[2] << ", " << lcar_ << "};\n//\n";
		return text.str();
	}
}
void Point::setIndex(const int& index)
{
	index_ = index;
}
void Point::setName(const std::string& name)
{
	name_ = name;
}
void Point::setX(const double& x)
{
	coordinates_[0] = x;
}
void Point::setY(const double& y)
{
	coordinates_[1] = y;
}
void Point::setZ(const double& z)
{
	coordinates_[2] = z;
}
void Point::setLcar(const double& lcar)
{
	lcar_ = lcar;
}
void Point::setDiscretization(const bool& discretization)
{
	discretization_ = discretization;
}
void Point::addNode(Node* node)
{
	node_ = node;
}