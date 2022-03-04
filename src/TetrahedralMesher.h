#pragma once

// External includes
#ifndef TETLIBRARY
#define TETLIBRARY
#endif

#include "../external_libraries/tetgen/tetgen.h"
#include "Mesher.h"
#include <string>

class TetrahedralMesher : public Mesher
{
protected:

    enum TetgenErrors {INPUT_MEMORY_ERROR=1, INTERNAL_ERROR=2, INVALID_GEOMETRY_ERROR=3};

public:

    TetrahedralMesher();

    ~TetrahedralMesher();

    void execute(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param) override;

    bool alphaShape(std::vector<Node*>& nodes, const double& alpha, const double& meanMeshLength) override;

    void buildInput(std::vector<Node*>& nodes, AnalysisParameters* param, tetgenio& in);

    void getFromContainer(tetgenio& tr);

    void setToContainer(tetgenio& tr);

    int generateTesselation(tetgenio& in, tetgenio& out);

    void deleteInContainer(tetgenio& tr);

    void deleteOutContainer(tetgenio& tr);

    void clearTetgenIO(tetgenio& tr);
};
