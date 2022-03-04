#include "PlaneSurface.h"

PlaneSurface::PlaneSurface() {}

PlaneSurface::PlaneSurface(const int& index, const std::string& name, LineLoop* lineLoop)
	: Surface(index, name, lineLoop) {}

PlaneSurface::~PlaneSurface() {}

std::string PlaneSurface::getGmshCode()
{
	std::stringstream text;
	text << name_ << " = news; Plane Surface(" << name_ << ") = {" << lineLoop_->getName() << "}; Physical Surface('" << name_ << "') = {" << name_ << "};\n//\n";
	return text.str();
}