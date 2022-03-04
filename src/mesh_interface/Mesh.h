#pragma once

#include <fstream>
#include <stdio.h>
#include "Geometry.h"
#include <unistd.h>

enum MeshAlgorithm
{
    AUTO,
    FRONT,
    DELAUNAY,
    ADAPT,
    PACK,
    QUAD,
    BAMG
};

std::pair<std::string, bool> createMesh(Geometry* geometry, const PartitionOfUnity& elementType, const MeshAlgorithm& algorithm = AUTO, std::string geofile = std::string(), const std::string& gmshPath = std::string(), const bool& plotMesh = true, const bool& showInfo = false);