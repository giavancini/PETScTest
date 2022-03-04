#include "Volume.h"

Volume::Volume() {}

Volume::Volume(const int index, const std::string name, SurfaceLoop* surfaceLoop)
	: index_(index),
	  numberOfElements_(0),
	  numberOfNodes_(0),
	  name_(name),
	  surfaceLoop_(surfaceLoop),
	  material_(nullptr) {}

Volume::~Volume() {}

int Volume::getIndex()
{
	return index_;
}

std::string Volume::getName()
{
	return name_;
}

SurfaceLoop* Volume::getSurfaceLoop()
{
	return surfaceLoop_;
}

Material* Volume::getMaterial()
{
	return material_;
}

const std::vector<BaseVolumeElement*>& Volume::getBaseElements()
{
	return baseElements_;
}

const std::vector<Element*>& Volume::getElements()
{
	return elements_;
}

const std::vector<Node*>& Volume::getNodes()
{
	return nodes_;
}

int Volume::getNumberOfNodes() const
{
	return nodes_.size();
}

int Volume::getNumberOfBaseElements() const
{
	return baseElements_.size();
}

int Volume::getNumberOfElements() const
{
	return elements_.size();
}

void Volume::setMaterial(Material* material)
{
	material_ = material;
}

void Volume::addElement(Element* element)
{
	elements_.push_back(element);
}

void Volume::addBaseElement(BaseVolumeElement* element)
{
	baseElements_.push_back(element);
}

void Volume::addNodes(const std::vector<Node*>& nodes)
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

std::string Volume::getGmshCode()
{
	std::stringstream text;
	text << name_ << " = newv; Volume(" << name_ << ") = {" << surfaceLoop_->getName() << "}; Physical Volume('" << name_ << "') = {" << name_ << "};\n//\n";
	return text.str();
}

void Volume::incrementNumberOfElements()
{
	numberOfElements_++;
}

void Volume::setNumberOfElements(const unsigned int numberOfElements)
{
	numberOfElements_ = numberOfElements;
}