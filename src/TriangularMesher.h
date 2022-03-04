#pragma once

// If SINGLE is defined when triangle.o is compiled, it should also be defined here
// If SINGLE is NOT defined in compilation, it should not be defined here.
// #define SINGLE

#ifdef   SINGLE
#define  REAL float
#else    // not SINGLE
#define  REAL double
#endif   // not SINGLE

// External includes
#ifndef TRILIBRARY
#define TRILIBRARY
#endif

#include "../external_libraries/triangle/triangle.h"
#include "Mesher.h"
#include <string>

extern "C"
{
	void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
	void trifree(void *);
}

class TriangularMesher : public Mesher
{
protected:

    enum TriangleErrors {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

public:

    TriangularMesher();

    virtual ~TriangularMesher();

    void execute(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param) override;

    bool alphaShape(std::vector<Node*>& nodes, const double& alpha, const double& meanMeshLength) override;

    void clearTrianglesList(struct triangulateio& tr);

    void buildInput(std::vector<Node*>& nodes, AnalysisParameters* param, struct triangulateio &in);

    void getFromContainer(struct triangulateio& tr);

    void setToContainer(struct triangulateio& tr);

    int generateTesselation(struct triangulateio& in, struct triangulateio& out);

    void deleteInContainer(struct triangulateio& tr);

    void deleteOutContainer(struct triangulateio& tr);
};

