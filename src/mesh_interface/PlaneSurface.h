#pragma once

#include "Surface.h"

class PlaneSurface : public Surface
{
public:
	PlaneSurface();

	PlaneSurface(const int& index, const std::string& name, LineLoop* lineLoop);

	~PlaneSurface();

	std::string getGmshCode() override;
};

