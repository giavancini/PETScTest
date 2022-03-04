#include "SurfaceLoop.h"

SurfaceLoop::SurfaceLoop() {}

SurfaceLoop::SurfaceLoop(const int& index, const std::string& name, std::vector<Surface*> surfaces)
{
	index_ = index;
	name_ = name;
	surfaces_.reserve(surfaces.size());
	for (Surface* surface : surfaces)
	{
		surfaces_.push_back(surface);
	}
		
}

SurfaceLoop::~SurfaceLoop() {}

int SurfaceLoop::getIndex()
{
	return index_;
}

std::string SurfaceLoop::getName()
{
	return name_;
}

Surface* SurfaceLoop::getSurface(const int& index)
{
	return surfaces_[index];
}

std::vector<Surface*> SurfaceLoop::getSurfaces()
{
	return surfaces_;
}

std::string SurfaceLoop::getGmshCode()
{
	std::stringstream text;
	text << name_ << " = newsl; Surface Loop(" << name_ << ") = {";
	for (size_t i = 0; i<surfaces_.size(); i++)
	{
		text << surfaces_[i]->getName();
		if (i != (surfaces_.size() - 1))
			text << ", ";
	}
	text << "};\n//\n";
	return text.str();
}

void SurfaceLoop::setIndex(const int& index)
{
	index_ = index;
}

void SurfaceLoop::setName(const std::string& name)
{
	name_ = name;
}

void SurfaceLoop::verification()
{
    if (surfaces_.size() < 4)
    {
        std::stringstream text;
        text << "The surfaces ";
        for (size_t j = 0; j < surfaces_.size(); j++)
        {
            text << surfaces_[j]->getName();
            if (j != (surfaces_.size() - 1))
                text << ", ";
        }
        text << " do not form a closed volume." << std::endl;
        std::cout << text.str();
    }
}