#include "TriangularMesher.h"

TriangularMesher::TriangularMesher()
    : Mesher() {}

TriangularMesher::~TriangularMesher() {}

void TriangularMesher::execute(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    executePreMeshingProcesses(nodes, elements, param);

    //Creating the containers for the input and output
    struct triangulateio in;
    struct triangulateio out;
    clearTrianglesList(out);

    buildInput(nodes, param, in);

    int triangle_error = generateTesselation(in, out);

    setToContainer(out);

    executePostMeshingProcesses(nodes, elements, param);

    deleteInContainer(in);
    deleteOutContainer(out);
}

bool TriangularMesher::alphaShape(std::vector<Node*>& nodes, const double& alpha, const double& meanMeshLength)
{
    BaseSurfaceElement* triangle = new BaseSurfaceElement(0, ParametricSurfaceElement::T3, nodes);

    double circumRadius = triangle->getRadius();

    delete triangle;

    double alphaRadius = alpha * meanMeshLength;

    if (circumRadius <= 0.0)
        return false; // degenerated element or a sliver
    else if (circumRadius < alphaRadius)
        return true; // passed the test
    else
        return false; // probably a big or distorted element
}

void TriangularMesher::clearTrianglesList(struct triangulateio& tr)
{
    tr.pointlist                  = (REAL*) NULL;
    tr.pointattributelist         = (REAL*) NULL;
    tr.pointmarkerlist            = (int*) NULL;
    tr.numberofpoints             = 0;
    tr.numberofpointattributes    = 0;

    tr.trianglelist               = (int*) NULL;
    tr.triangleattributelist      = (REAL*) NULL;
    tr.trianglearealist           = (REAL*) NULL;
    tr.neighborlist               = (int*) NULL;
    tr.numberoftriangles          = 0;
    tr.numberofcorners            = 3; //for three node triangles
    tr.numberoftriangleattributes = 0;

    tr.segmentlist                = (int*) NULL;
    tr.segmentmarkerlist          = (int*) NULL;
    tr.numberofsegments           = 0;

    tr.holelist                   = (REAL*) NULL;
    tr.numberofholes              = 0;

    tr.regionlist                 = (REAL*) NULL;
    tr.numberofregions            = 0;

    tr.edgelist                   = (int*) NULL;
    tr.edgemarkerlist             = (int*) NULL;
    tr.normlist                   = (REAL*) NULL;
    tr.numberofedges              = 0;
}

void TriangularMesher::buildInput(std::vector<Node*>& nodes, AnalysisParameters* param, struct triangulateio &in)
{
    int dimension = param->getDimension();

    inMesh_.createPointList(nodes.size(), dimension);

    double* pointList = inMesh_.getPointList();
    int& npoints = inMesh_.getNumberOfPoints();

    int base = 0;
    for (Node* const& node : nodes)
    {
        for (int j = 0; j < dimension; j++)
        {
            double coord = node->getDegreeOfFreedom(j)->getCurrentValue();
            pointList[base+j] = coord;
        }
        base += dimension;
    }

    clearTrianglesList(in);
    getFromContainer(in);
}

void TriangularMesher::getFromContainer(struct triangulateio& tr)
{
    //get pointers
    tr.pointlist        = inMesh_.getPointList();
    tr.trianglelist     = inMesh_.getElementList();
    tr.trianglearealist = inMesh_.getElementSizeList();
    tr.neighborlist     = inMesh_.getElementNeighbourList();

    if( inMesh_.getNumberOfPoints() != 0 )
      tr.numberofpoints = inMesh_.getNumberOfPoints();

    if( inMesh_.getNumberOfElements() != 0 )
      tr.numberoftriangles = inMesh_.getNumberOfElements();
}

void TriangularMesher::setToContainer(struct triangulateio& tr)
{
    //set pointers
    outMesh_.setPointList(tr.pointlist);
    outMesh_.setElementList(tr.trianglelist);
    outMesh_.setElementSizeList(tr.trianglearealist);
    outMesh_.setElementNeighbourList(tr.neighborlist);

    // copy the numbers
    if( tr.numberofpoints != 0 ){
      outMesh_.setNumberOfPoints(tr.numberofpoints);
    }

    if( tr.numberoftriangles != 0 ){
      outMesh_.setNumberOfElements(tr.numberoftriangles);
    }
}

int TriangularMesher::generateTesselation(struct triangulateio& in, struct triangulateio& out)
{
    std::string triflags = "znQP";
	char *triswitches = new char[triflags.size()+1];
	std::strcpy(triswitches, triflags.c_str());

    int triangle_error = 0;

    try
    {
      triangulate(triswitches, &in, &out, (struct triangulateio *)NULL);
    }

    catch(int error_code)
    {
        switch(TriangleErrors(error_code))
	    {
	        case INPUT_MEMORY_ERROR:
                triangle_error=1;
	            break;
	        case INTERNAL_ERROR:
                triangle_error=2;
	            break;
	        case INVALID_GEOMETRY_ERROR:
                triangle_error=3;
	            break;
	        default:
                triangle_error=0;
	            break;
	    }
    }

    delete [] triswitches;

	return triangle_error;
}

void TriangularMesher::deleteInContainer(struct triangulateio& tr)
{
    clearTrianglesList(tr);

    inMesh_.finalize();
}

void TriangularMesher::deleteOutContainer(struct triangulateio& tr)
{
    if (tr.numberoftriangles)
    {
        if (tr.trianglelist)
            trifree(tr.trianglelist);
        if (tr.triangleattributelist)
            trifree(tr.triangleattributelist);
        if (tr.trianglearealist)
            trifree(tr.trianglearealist);
        if (tr.neighborlist)
            trifree(tr.neighborlist);
    }

    if (tr.segmentlist)
        trifree(tr.segmentlist);
    if (tr.segmentmarkerlist)
        trifree(tr.segmentmarkerlist);
    
    if (tr.holelist)
    {
        delete [] tr.holelist;
        tr.numberofholes = 0;
    }

    if (tr.regionlist)
    {
        delete [] tr.regionlist;
        tr.numberofregions = 0;
    }
    
    if (tr.edgelist)
        trifree(tr.edgelist);
    if (tr.edgemarkerlist)
        trifree(tr.edgemarkerlist);
    if (tr.normlist)
        trifree(tr.normlist);

    if (tr.numberofpoints)
    {
      if (tr.pointlist)
        trifree(tr.pointlist);
      if (tr.pointattributelist)
        trifree(tr.pointattributelist);
      if (tr.pointmarkerlist)
        trifree(tr.pointmarkerlist);
    }

    clearTrianglesList(tr);

    outMesh_.finalize();
}