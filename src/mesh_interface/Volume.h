#pragma once

#include "SurfaceLoop.h"
#include "../Element.h"
#include "../BaseVolumeElement.h"
#include "../Material.h"

class Volume
{
public:
	Volume();

	Volume(const int index, const std::string name, SurfaceLoop* surfaceLoop);

	~Volume();

	int getIndex();

	std::string getName();

	SurfaceLoop* getSurfaceLoop();

	Material* getMaterial();

	const std::vector<BaseVolumeElement*>& getBaseElements();

	const std::vector<Element*>& getElements();

	const std::vector<Node*>& getNodes();

	int getNumberOfNodes() const;

	int getNumberOfBaseElements() const;

	int getNumberOfElements() const;

	void setMaterial(Material* material);

	void addElement(Element* element);

	void addBaseElement(BaseVolumeElement* element);

	void addNodes(const std::vector<Node*>& nodes);

	std::string getGmshCode();

	void incrementNumberOfElements();

	void setNumberOfElements(const unsigned int numberOfElements);

private:
	int index_;
	unsigned int numberOfElements_;
	unsigned int numberOfNodes_;
	std::string name_;
	SurfaceLoop* surfaceLoop_;
	Material* material_;
	std::vector<Element*> elements_;
	std::vector<Node*> nodes_;
	std::vector<BaseVolumeElement*> baseElements_;

public:
	friend class FluidDomain;
	friend class SolidDomain;
};

