//Fluid problem

double h = 4.0e-3;

Geometry *fluid_geo = new Geometry(0);
Point *p0 = fluid_geo->addPoint({0.0 , 0.0  ,  0.0});
Point *p1 = fluid_geo->addPoint({0.1 , 0.0  ,  0.0});
Point *p2 = fluid_geo->addPoint({0.1 , 0.08 ,  0.0});
Point *p3 = fluid_geo->addPoint({0.0 , 0.08 ,  0.0});
Point *p4 = fluid_geo->addPoint({0.0 , 0.0  , -0.05});
Point *p5 = fluid_geo->addPoint({0.1 , 0.0  , -0.05});
Point *p6 = fluid_geo->addPoint({0.1 , 0.08 , -0.05});
Point *p7 = fluid_geo->addPoint({0.0 , 0.08 , -0.05});
Line  *l0 = fluid_geo->addLine({p0, p1});
Line  *l1 = fluid_geo->addLine({p1, p2});
Line  *l2 = fluid_geo->addLine({p2, p3});
Line  *l3 = fluid_geo->addLine({p3, p0});
Line  *l4 = fluid_geo->addLine({p4, p5});
Line  *l5 = fluid_geo->addLine({p5, p6});
Line  *l6 = fluid_geo->addLine({p6, p7});
Line  *l7 = fluid_geo->addLine({p7, p4});
Line  *l8 = fluid_geo->addLine({p0, p4});
Line  *l9 = fluid_geo->addLine({p1, p5});
Line *l10 = fluid_geo->addLine({p2, p6});
Line *l11 = fluid_geo->addLine({p3, p7});
Surface *s0 = fluid_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = fluid_geo->addPlaneSurface({l4, l5, l6, l7});
Surface *s2 = fluid_geo->addPlaneSurface({l0, l9, -(*l4), -(*l8)});
Surface *s3 = fluid_geo->addPlaneSurface({-(*l2), l10, l6, -(*l11)});
Surface *s4 = fluid_geo->addPlaneSurface({l9, l5, -(*l10), -(*l1)});
Surface *s5 = fluid_geo->addPlaneSurface({l8, -(*l7), -(*l11), l3});
Volume  *v0 = fluid_geo->addVolume({s0, s1, s2, s3, s4, s5});
fluid_geo->transfiniteLine({l0}, 0.1 / h + 1);
fluid_geo->transfiniteLine({l1}, 0.08 / h + 1);
fluid_geo->transfiniteLine({l2}, 0.1 / h + 1);
fluid_geo->transfiniteLine({l3}, 0.08 / h + 1);
fluid_geo->transfiniteLine({l4}, 0.1 / h + 1);
fluid_geo->transfiniteLine({l5}, 0.08 / h + 1);
fluid_geo->transfiniteLine({l6}, 0.1 / h + 1);
fluid_geo->transfiniteLine({l7}, 0.08 / h + 1);
fluid_geo->transfiniteLine({l8}, 0.05 / h + 1);
fluid_geo->transfiniteLine({l9}, 0.05 / h + 1);
fluid_geo->transfiniteLine({l10}, 0.05 / h + 1);
fluid_geo->transfiniteLine({l11}, 0.05 / h + 1);
fluid_geo->transfiniteSurface({s4}, "Right");
fluid_geo->addDirichletBoundaryCondition({s2}, Y, 0.0);
fluid_geo->addDirichletBoundaryCondition({s5}, X, 0.0);
fluid_geo->addDirichletBoundaryCondition({s0}, Z, 0.0);
fluid_geo->addDirichletBoundaryCondition({s1}, Z, 0.0);
fluid_geo->addInterfaceBoundaryCondition({s4});
Material *fluid_mat = new NewtonianFluid(1.0, 1000.0);
FluidDomain *fluid_problem = new FluidDomain(fluid_geo);
fluid_problem->applyMaterial({v0}, fluid_mat);
fluid_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
fluid_problem->setGravity(0.0, -9.81, 0.0);
fluid_problem->setMeshLength(h);

//Solid problem
Geometry *solid_geo = new Geometry(1);
Point *p8 =  solid_geo->addPoint({0.1   , 0.0  ,  0.0});
Point *p9 =  solid_geo->addPoint({0.112 , 0.0  ,  0.0});
Point *p10 = solid_geo->addPoint({0.112 , 0.1  ,  0.0});
Point *p11 = solid_geo->addPoint({0.1   , 0.1  ,  0.0});
Point *p12 = solid_geo->addPoint({0.1   , 0.08 ,  0.0});
Point *p13 = solid_geo->addPoint({0.1   , 0.0  , -0.05});
Point *p14 = solid_geo->addPoint({0.112 , 0.0  , -0.05});
Point *p15 = solid_geo->addPoint({0.112 , 0.1  , -0.05});
Point *p16 = solid_geo->addPoint({0.1   , 0.1  , -0.05});
Point *p17 = solid_geo->addPoint({0.1   , 0.08 , -0.05});
Line  *l12 = solid_geo->addLine({p8, p9});
Line  *l13 = solid_geo->addLine({p9, p10});
Line  *l14 = solid_geo->addLine({p10, p11});
Line  *l15 = solid_geo->addLine({p11, p12});
Line  *l16 = solid_geo->addLine({p12, p8});
Line  *l17 = solid_geo->addLine({p13, p14});
Line  *l18 = solid_geo->addLine({p14, p15});
Line  *l19 = solid_geo->addLine({p15, p16});
Line  *l20 = solid_geo->addLine({p16, p17});
Line  *l21 = solid_geo->addLine({p17, p13});
Line  *l22 = solid_geo->addLine({p8, p13});
Line  *l23 = solid_geo->addLine({p9, p14});
Line  *l24 = solid_geo->addLine({p10, p15});
Line  *l25 = solid_geo->addLine({p11, p16});
Line  *l26 = solid_geo->addLine({p12, p17});
Surface *s6  = solid_geo->addPlaneSurface({l12, l13, l14, l15, l16});
Surface *s7  = solid_geo->addPlaneSurface({l17, l18, l19, l20, l21});
Surface *s8  = solid_geo->addPlaneSurface({l12, l23, -(*l17), -(*l22)});
Surface *s9  = solid_geo->addPlaneSurface({-(*l14), l24, l19, -(*l25)});
Surface *s10 = solid_geo->addPlaneSurface({l23, l18, -(*l24), -(*l13)});
Surface *s11 = solid_geo->addPlaneSurface({l26, -(*l20), -(*l25), l15});
Surface *s12 = solid_geo->addPlaneSurface({l22, -(*l21), -(*l26), l16});
Volume *v1   = solid_geo->addVolume({s6, s7, s8, s9, s10, s11, s12});
solid_geo->transfiniteLine({l12}, 0.012 / h + 1);
solid_geo->transfiniteLine({l13}, 0.1 / h + 1);
solid_geo->transfiniteLine({l14}, 0.012 / h + 1);
solid_geo->transfiniteLine({l15}, 0.02 / h + 1);
solid_geo->transfiniteLine({l16}, 0.08 / h + 1);
solid_geo->transfiniteLine({l17}, 0.012 / h + 1);
solid_geo->transfiniteLine({l18}, 0.1 / h + 1);
solid_geo->transfiniteLine({l19}, 0.012 / h + 1);
solid_geo->transfiniteLine({l20}, 0.02 / h + 1);
solid_geo->transfiniteLine({l21}, 0.08 / h + 1);
solid_geo->transfiniteLine({l22}, 0.05 / h + 1);
solid_geo->transfiniteLine({l23}, 0.05 / h + 1);
solid_geo->transfiniteLine({l24}, 0.05 / h + 1);
solid_geo->transfiniteLine({l25}, 0.05 / h + 1);
solid_geo->transfiniteLine({l26}, 0.05 / h + 1);
solid_geo->transfiniteSurface({s6},  "Left", {p8, p9, p10, p11, p12});
solid_geo->transfiniteSurface({s7},  "Left", {p13, p14, p15, p16, p17});
solid_geo->transfiniteSurface({s8},  "Right", {p8, p9, p13, p14});
solid_geo->transfiniteSurface({s9},  "Right", {p11, p10, p15, p16});
solid_geo->transfiniteSurface({s10}, "Right", {p9, p14, p15, p10});
solid_geo->transfiniteSurface({s11}, "Right", {p12, p17, p16, p11});
solid_geo->transfiniteSurface({s12}, "Right", {p8, p13, p17, p12});
solid_geo->addDirichletBoundaryCondition({s8}, ALL, 0.0);
// solid_geo->addDirichletBoundaryCondition({s6}, Z, 0.0);
// solid_geo->addDirichletBoundaryCondition({s7}, Z, 0.0);
solid_geo->addInterfaceBoundaryCondition({s12});
Material *solid_mat = new ElasticSolid(8.0e6, 0.0, 2500.0);
SolidDomain *solid_problem = new SolidDomain(solid_geo);
solid_problem->applyMaterial({v1}, solid_mat);
solid_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
solid_problem->setGravity(0.0, 0.0, 0.0);

//Coupled problem
CoupledDomain *coupled_problem = new CoupledDomain(fluid_problem, solid_problem);
coupled_problem->setNumberOfSteps(1500);
coupled_problem->setDeltat(0.001);
coupled_problem->setMaxNonlinearIterations(6);
coupled_problem->setNonlinearTolerance(1.0e-4);
coupled_problem->setSpectralRadius(0.5);
coupled_problem->solveMonolithicCoupledProblem();