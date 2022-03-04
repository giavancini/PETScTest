#pragma once
#include "LineLoop.h"
#include "../Element.h"
#include "../BaseSurfaceElement.h"
#include "../Material.h"
#include <algorithm>

class Surface
{
public:
	Surface();

	Surface(const int index, const std::string name, LineLoop* lineLoop);

	~Surface();

	int getIndex();

	std::string getName();

	LineLoop* getLineLoop();

	Material* getMaterial();

	const std::vector<BaseSurfaceElement*>& getBaseElements();

	const std::vector<Element*>& getElements();

	const std::vector<Node*>& getNodes();

	int getNumberOfNodes() const;

	int getNumberOfBaseElements() const;

	int getNumberOfElements() const;

	void setMaterial(Material* material);

	void addElement(Element* element);

	void addBaseElement(BaseSurfaceElement* element);

	virtual std::string getGmshCode();

	void addNodes(const std::vector<Node*>& nodes);

	void clearAllElements();

	void clearAllNodes();

	void removeElement(const int& index);

	void incrementNumberOfElements();

	void setNumberOfElements(const unsigned int numberOfElements);

protected:
	int index_;
	unsigned int numberOfElements_;
	unsigned int numberOfNodes_;
	std::string name_;
	LineLoop* lineLoop_;
	Material* material_;
	std::vector<BaseSurfaceElement*> baseElements_;
	std::vector<Element*> elements_;
	std::vector<Node*> nodes_;

public:
	friend class FluidDomain;
	friend class SolidDomain;
};

