#pragma once

#include "Point.h"
#include "../Element.h"
#include "../BaseLineElement.h"
#include "../Material.h"

class Line
{
public:
	Line();

	Line(const int index, const std::string name, std::vector<Point*> points, const bool discretization = true);

	~Line();

	Line* operator-();

	int getIndex();

	std::string getName();

	Point* getInitialPoint();

	Point* getEndPoint();

	bool getDiscretization();

	const std::vector<Node*>& getNodes();

	const std::vector<BaseLineElement*>& getBaseElements();

	const std::vector<Element*>& getElements();

	int getNumberOfNodes() const;

	int getNumberOfBaseElements() const;

	int getNumberOfElements() const;

	Material* getMaterial();

	std::string virtual getGmshCode();

	void setIndex(const int& index);

	void setName(const std::string& name);

	void setInitialPoint(Point& point);

	void setEndPoint(Point& point);

	void setDiscretization(const bool& discretization);

	void addNodes(const std::vector<Node*>& nodes);

	void addBaseElement(BaseLineElement* element);

	void addElement(Element* element);

	void setMaterial(Material* material);

	void incrementNumberOfElements();

	void setNumberOfElements(const unsigned int numberOfElements);

protected:
	int index_;
	unsigned int numberOfElements_;
	unsigned int numberOfNodes_;
	std::string name_;
	std::vector<Point*> points_;
	bool discretization_;
	Material* material_;
	std::vector<Node*> nodes_;
	std::vector<BaseLineElement*> baseElements_;
	std::vector<Element*> elements_;

public:
	friend class FluidDomain;
	friend class SolidDomain;
};

