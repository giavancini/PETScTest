#include "TetrahedralMesher.h"

TetrahedralMesher::TetrahedralMesher()
    : Mesher() {}

TetrahedralMesher::~TetrahedralMesher() {}

void TetrahedralMesher::execute(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    executePreMeshingProcesses(nodes, elements, param);

    //Creating the containers for the input and output
    tetgenio in;
    tetgenio out;

    buildInput(nodes, param, in);

    int tetgen_error = generateTesselation(in, out);

    setToContainer(out);

    executePostMeshingProcesses(nodes, elements, param);

    deleteInContainer(in);
    deleteOutContainer(out);
}

bool TetrahedralMesher::alphaShape(std::vector<Node*>& nodes, const double& alpha, const double& meanMeshLength)
{
    BaseVolumeElement* tetrahedron = new BaseVolumeElement(0, ParametricVolumeElement::TET4, nodes);

    double circumRadius = tetrahedron->getRadius();

    delete tetrahedron;

    double alphaRadius = alpha * meanMeshLength;

    if (circumRadius <= 0.0)
        return false; // degenerated element or a sliver
    else if (circumRadius < alphaRadius)
        return true; // passed the test
    else
        return false; // probably a big or distorted element
}

void  TetrahedralMesher::buildInput(std::vector<Node*>& nodes, AnalysisParameters* param, tetgenio &in)
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

    in.firstnumber = 0;
    in.mesh_dim = 3;
    getFromContainer(in);
}

void TetrahedralMesher::getFromContainer(tetgenio& tr)
{
    //get pointers
    tr.pointlist             = inMesh_.getPointList();
    tr.tetrahedronlist       = inMesh_.getElementList();
    tr.tetrahedronvolumelist = inMesh_.getElementSizeList();
    tr.neighborlist          = inMesh_.getElementNeighbourList();

    if( inMesh_.getNumberOfPoints() != 0 )
      tr.numberofpoints = inMesh_.getNumberOfPoints();

    if( inMesh_.getNumberOfElements() != 0 )
      tr.numberoftetrahedra = inMesh_.getNumberOfElements();
}

void TetrahedralMesher::setToContainer(tetgenio& tr)
{
    //set pointers
    outMesh_.setPointList(tr.pointlist);
    outMesh_.setElementList(tr.tetrahedronlist);
    outMesh_.setElementSizeList(tr.tetrahedronvolumelist);
    outMesh_.setElementNeighbourList(tr.neighborlist);

    // copy the numbers
    if( tr.numberofpoints != 0 ){
      outMesh_.setNumberOfPoints(tr.numberofpoints);
    }

    if( tr.numberoftetrahedra != 0 ){
      outMesh_.setNumberOfElements(tr.numberoftetrahedra);
    }
}

int TetrahedralMesher::generateTesselation(tetgenio& in, tetgenio& out)
{
    std::string tetflags = "znJQF";
	char *tetswitches = new char[tetflags.size()+1];
	std::strcpy(tetswitches, tetflags.c_str());

    int tetgen_error = 0;

    try
    {
      tetrahedralize(tetswitches, &in, &out);
    }

    catch(int error_code)
    {
        switch(TetgenErrors(error_code))
	    {
	        case INPUT_MEMORY_ERROR:
                tetgen_error=1;
	            break;
	        case INTERNAL_ERROR:
                tetgen_error=2;
	            break;
	        case INVALID_GEOMETRY_ERROR:
                tetgen_error=3;
	            break;
	        default:
                tetgen_error=0;
	            break;
	    }
    }

    delete [] tetswitches;

	return tetgen_error;
}

void TetrahedralMesher::deleteInContainer(tetgenio& tr)
{
    inMesh_.finalize();
    clearTetgenIO(tr);
}

void TetrahedralMesher::deleteOutContainer(tetgenio& tr)
{
    tr.deinitialize();
    tr.initialize();

    outMesh_.finalize();
}

void TetrahedralMesher::clearTetgenIO(tetgenio& tr)
  {
    tr.pointlist                     = (REAL*) NULL;
    tr.numberofpoints                = 0;
    tr.numberofpointattributes       = 0;

    tr.tetrahedronlist               = (int*) NULL;
    tr.tetrahedronvolumelist         = (REAL*) NULL;
    tr.neighborlist                  = (int*) NULL;
    tr.numberoftetrahedra            = 0;
  }