double h = 0.0006;
double l = 0.02;
double a = 0.0025;
double pi = 3.14;

Geometry *fluid_geo = new Geometry(0);
Point *p0 = fluid_geo->addPoint({0.0   , 0.0 });
Point *p1 = fluid_geo->addPoint({2.0*l , 0.0 });
Point *p2 = fluid_geo->addPoint({2.0*l , 4.0*l});
Point *p3 = fluid_geo->addPoint({0.0   , 4.0*l});
Point *p4 = fluid_geo->addPoint({1.0*l , 3.0*l}, 1.0, false);
Point *p5 = fluid_geo->addPoint({1.0*l+a , 3.0*l });
// Point *p6 = fluid_geo->addPoint({1.0*l , 3.0*l+a });
Point *p7 = fluid_geo->addPoint({1.0*l-a , 3.0*l });
// Point *p8 = fluid_geo->addPoint({1.0*l , 3.0*l-a });
Line *l0 = fluid_geo->addLine({p0, p1});
Line *l1 = fluid_geo->addLine({p1, p2});
Line *l2 = fluid_geo->addLine({p2, p3});
Line *l3 = fluid_geo->addLine({p3, p0});
Line *l4 = fluid_geo->addCircle({p5, p4, p7});
Line *l5 = fluid_geo->addCircle({p7, p4, p5});
// Line *l6 = fluid_geo->addCircle({p7, p4, p8});
// Line *l7 = fluid_geo->addCircle({p8, p4, p5});
Surface *s0 = fluid_geo->addPlaneSurface({l0, l1, l2, l3, -(*l4), -(*l5)});
fluid_geo->transfiniteLine({l0}, 2.0*l / h + 1);
fluid_geo->transfiniteLine({l1}, 4.0*l / h + 1);
fluid_geo->transfiniteLine({l2}, 2.0*l / h + 1);
fluid_geo->transfiniteLine({l3}, 4.0*l / h + 1);
fluid_geo->transfiniteLine({l4}, (2.0*pi*a/2.0) / h + 1);
fluid_geo->transfiniteLine({l5}, (2.0*pi*a/2.0) / h + 1);
// fluid_geo->transfiniteLine({l6}, (2.0*pi*a/4.0) / h + 1);
// fluid_geo->transfiniteLine({l7}, (2.0*pi*a/4.0) / h + 1);
fluid_geo->addDirichletBoundaryCondition({l0}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({l1}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({l3}, ALL, 0.0);
fluid_geo->addInterfaceBoundaryCondition({l4});
fluid_geo->addInterfaceBoundaryCondition({l5});
// fluid_geo->addInterfaceBoundaryCondition({l6});
// fluid_geo->addInterfaceBoundaryCondition({l7});
Material *fluid_mat = new NewtonianFluid(0.1, 1000.0);
FluidDomain *fluid_problem = new FluidDomain(fluid_geo);
fluid_problem->applyMaterial({s0}, fluid_mat);
fluid_problem->generateMesh(T3, FRONT, "teste", "", true, false);
fluid_problem->setGravity(0.0, -9.81, 0.0);
fluid_problem->setMeshLength(h);

//Solid problem
Geometry *solid_geo = new Geometry(1);
Point *p9 = solid_geo->addPoint({1.0*l , 3.0*l}, 1.0, false);
Point *p10 = solid_geo->addPoint({1.0*l+a , 3.0*l });
// Point *p11 = solid_geo->addPoint({1.0*l , 3.0*l+a });
Point *p12 = solid_geo->addPoint({1.0*l-a , 3.0*l });
// Point *p13 = solid_geo->addPoint({1.0*l , 3.0*l-a });
Line *l8 = solid_geo->addCircle({p10, p9, p12});
Line *l9 = solid_geo->addCircle({p12, p9, p10});
// Line *l10 = solid_geo->addCircle({p12, p9, p13});
// Line *l11 = solid_geo->addCircle({p13, p9, p10});
Surface *s1 = solid_geo->addPlaneSurface({l8, l9});
solid_geo->transfiniteLine({l8}, (2.0*pi*a/2.0) / h + 1);
solid_geo->transfiniteLine({l9}, (2.0*pi*a/2.0) / h + 1);
// solid_geo->transfiniteLine({l10}, (2.0*pi*a/4.0) / h + 1);
// solid_geo->transfiniteLine({l11}, (2.0*pi*a/4.0) / h + 1);
solid_geo->addInterfaceBoundaryCondition({l8});
solid_geo->addInterfaceBoundaryCondition({l9});
// solid_geo->addInterfaceBoundaryCondition({l10});
// solid_geo->addInterfaceBoundaryCondition({l11});
Material *solid_mat = new ElasticSolid(1.0e16, 0.35, 1200.0);
SolidDomain *solid_problem = new SolidDomain(solid_geo);
solid_problem->applyMaterial({s1}, solid_mat);
solid_problem->generateMesh(T3, FRONT, "teste", "", true, false);
solid_problem->setGravity(0.0, -9.81, 0.0);

//Coupled problem
CoupledDomain *coupled_problem = new CoupledDomain(fluid_problem, solid_problem);
coupled_problem->setNumberOfSteps(1200);
coupled_problem->setDeltat(0.001);
coupled_problem->setMaxNonlinearIterations(6);
coupled_problem->setNonlinearTolerance(1.0e-6);
coupled_problem->setSpectralRadius(0.8);
coupled_problem->solveMonolithicCoupledPFEMProblem();