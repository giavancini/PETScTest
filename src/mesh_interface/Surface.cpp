#include "Surface.h"

Surface::Surface() {}

Surface::Surface(const int index, const std::string name, LineLoop* lineLoop)
	: index_(index),
	  numberOfElements_(0),
	  numberOfNodes_(0),
	  name_(name),
	  lineLoop_(lineLoop),
	  material_(nullptr) {}

Surface::~Surface() {}

int Surface::getIndex()
{
	return index_;
}

std::string Surface::getName()
{
	return name_;
}

LineLoop* Surface::getLineLoop()
{
	return lineLoop_;
}

Material* Surface::getMaterial()
{
	return material_;
}

const std::vector<BaseSurfaceElement*>& Surface::getBaseElements()
{
	return baseElements_;
}

const std::vector<Element*>& Surface::getElements()
{
	return elements_;
}

const std::vector<Node*>& Surface::getNodes()
{
	return nodes_;
}

int Surface::getNumberOfNodes() const
{
	return nodes_.size();
}

int Surface::getNumberOfBaseElements() const
{
	return baseElements_.size();
}

int Surface::getNumberOfElements() const
{
	return elements_.size();
}

void Surface::setMaterial(Material* material)
{
	material_ = material;
}

void Surface::addElement(Element* element)
{
	elements_.push_back(element);
}

void Surface::addBaseElement(BaseSurfaceElement* element)
{
	baseElements_.push_back(element);
}

std::string Surface::getGmshCode()
{
	std::stringstream text;
	text << name_ << " = news; Surface(" << name_ << ") = {" << lineLoop_->getName() << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	return text.str();
}

void Surface::addNodes(const std::vector<Node*>& nodes)
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

void Surface::clearAllElements()
{
	elements_.clear();
	baseElements_.clear();
	std::vector<Element*>().swap(elements_);
	std::vector<BaseSurfaceElement*>().swap(baseElements_);
}

void Surface::clearAllNodes()
{
	nodes_.clear();
	std::vector<Node*>().swap(nodes_);
}

void Surface::removeElement(const int& index)
{
	elements_.erase(elements_.begin() + index);
	baseElements_.erase(baseElements_.begin() + index);
}

void Surface::incrementNumberOfElements()
{
	numberOfElements_++;
}

void Surface::setNumberOfElements(const unsigned int numberOfElements)
{
	numberOfElements_ = numberOfElements;
}