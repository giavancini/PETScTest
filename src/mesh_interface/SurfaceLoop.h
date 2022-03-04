#pragma once

#include "PlaneSurface.h"

class SurfaceLoop
{
public:
	SurfaceLoop();

	SurfaceLoop(const int& index, const std::string& name, std::vector<Surface*> surfaces);

	~SurfaceLoop();

	int getIndex();

	std::string getName();

	Surface* getSurface(const int& index);

	std::vector<Surface*> getSurfaces();

	std::string getGmshCode();

	void setIndex(const int& index);

	void setName(const std::string& name);

	void verification();

private:
	int index_;
	std::string name_;
	std::vector<Surface*> surfaces_;
};

