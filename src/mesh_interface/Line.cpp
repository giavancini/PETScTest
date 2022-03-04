#include "Line.h"

Line::Line() {}

Line::Line(const int index, const std::string name, std::vector<Point*> points, const bool discretization)
	: index_(index),
	  numberOfElements_(0),
	  numberOfNodes_(0),
	  name_(name),
	  points_(points),
	  discretization_(discretization),
      material_(nullptr) {}

Line::~Line() {}

Line* Line::operator-()
{
	Line* copy = new Line(index_, name_, {points_});
	copy->setName("-" + name_);
	return copy;
}

int Line::getIndex()
{
	return index_;
}
std::string Line::getName()
{
	return name_;
}
Point* Line::getInitialPoint()
{
	return points_[0];
}
Point* Line::getEndPoint()
{
	int last = points_.size() -1;
	return points_[last];
}
bool Line::getDiscretization()
{
	return discretization_;
}
const std::vector<Node*>& Line::getNodes()
{
	return nodes_;
}

const std::vector<BaseLineElement*>& Line::getBaseElements()
{
	return baseElements_;
}

const std::vector<Element*>& Line::getElements()
{
	return elements_;
}

int Line::getNumberOfNodes() const
{
	return nodes_.size();
}

int Line::getNumberOfBaseElements() const
{
	return baseElements_.size();
}

int Line::getNumberOfElements() const
{
	return elements_.size();
}

Material* Line::getMaterial()
{
	return material_;
}

std::string Line::getGmshCode()
{
	std::stringstream text;
	if (discretization_) {
		text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
			<< "}; Physical Line('" << name_ << "') = {" << name_ << "};\n//\n";
		return text.str();
	}
	else {
		text << name_ << " = newl; Line(" << name_ << ") = {" << points_[0]->getName() << ", " << points_[1]->getName()
			<< "};\n//\n";
		return text.str();
	}
}
void Line::setIndex(const int& index)
{
	index_ = index;
}
void Line::setName(const std::string& name)
{
	name_ = name;
}
void Line::setInitialPoint(Point& point)
{
	//points_[0] = point;
}
void Line::setEndPoint(Point& point)
{
	//points_[1] = point;
}
void Line::setDiscretization(const bool& discretization)
{
	discretization_ = discretization;
}
void Line::setMaterial(Material* const material)
{
	material_ = material;
}
void Line::addNodes(const std::vector<Node*>& nodes)
{
	for (auto& node1 : nodes)
	{
		bool notDuplicate = true;
		for (auto& node2 : nodes_)
		{
			if (node1->getIndex() == node2->getIndex())
			{
				notDuplicate = false;
				break;
			}
		}
		if (notDuplicate)
			nodes_.push_back(node1);
	}
}
void Line::addBaseElement(BaseLineElement* element)
{
	baseElements_.push_back(element);
}

void Line::addElement(Element* element)
{
	elements_.push_back(element);
}

void Line::incrementNumberOfElements()
{
	numberOfElements_++;
}

void Line::setNumberOfElements(const unsigned int numberOfElements)
{
	numberOfElements_ = numberOfElements;
}