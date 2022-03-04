#pragma once

namespace quadratures
{
    void lineQuadrature(const int& numberOfPoints, double**& xsi, double*& weight);

    void triangleQuadrature(const int& numberOfPoints, double**& xsi, double*& weight);

    void tetrahedronQuadrature(const int& numberOfPoints, double**& xsi, double*& weight);

    void pyramidQuadrature(const int& numberOfPoints, double**& xsi, double*& weight);
}