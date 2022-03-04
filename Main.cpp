#include "src/CoupledDomain.h"

int main(int argc, char **args)
{
    PetscInitialize(&argc, &args, (char*)0, NULL);

    //Solid Examples
    //#include "examples/plate_with_hole.h"
    //#include "examples/cooks_membrane.h"
    //#include "examples/nonlinear_clamped_beam.h"
    //#include "examples/3dnonlinear_clamped_beam.h"
    //#include "examples/cantilever_beam.h"
    //#include "examples/plane_strain_cantilever_beam.h"
    //#include "examples/3dplane_strain_cantilever_beam.h"
    //#include "examples/tensioned_bar.h"
    //#include "examples/3dtest_static.h"
    //#include "examples/3dtest_transient.h"

    //Fluid Examples
    //#include "examples/3dtank.h"
    //#include "examples/dam.h"
    //#include "examples/3ddam.h"
    //#include "examples/slosh.h"
    //#include "examples/3dslosh.h"
    //#include "examples/wave.h"
    //#include "examples/3dwave.h"
    #include "examples/damPFEM.h"
    //#include "examples/sloshPFEM.h"
    //#include "examples/waterDropPFEM.h"
    //#include "examples/twoWaterColumnsPFEM.h"
    //#include "examples/testTriangularMesher.hpp"
    //#include "examples/testTetrahedralMesher.h"
    //#include "examples/sloshPFEM3D.h"

    //FSI Examples
    //#include "examples/flexibleTank.h"
    //#include "examples/flexibleTankPFEM.h"
    //#include "examples/3dflexibleTank.h"
    //#include "examples/3dflexibleTankPFEM.h"
    //#include "examples/elastic_gatePFEM.h"
    //#include "examples/damElasticPFEM.h"
    //#include "examples/3ddamElasticPFEM.h"
    //#include "examples/damElasticPFEMHard.h"
    //#include "examples/fallingCylinderPFEM.h"
    //#include "examples/fillingElasticContainerPFEM.h"

    PetscFinalize();
	return 0;
}
