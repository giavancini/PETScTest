//Fluid problem

double h = 4.0e-3;

Geometry *fluid_geo = new Geometry(0);
Point *p0 = fluid_geo->addPoint({0.0 , 0.0  , 0.0 });
Point *p1 = fluid_geo->addPoint({0.1 , 0.0  , 0.0 });
Point *p2 = fluid_geo->addPoint({0.1 , 0.08 , 0.0});
Point *p3 = fluid_geo->addPoint({0.0 , 0.08 , 0.0});
Line *l0 = fluid_geo->addLine({p0, p1});
Line *l1 = fluid_geo->addLine({p1, p2});
Line *l2 = fluid_geo->addLine({p2, p3});
Line *l3 = fluid_geo->addLine({p3, p0});
Surface *s0 = fluid_geo->addPlaneSurface({l0, l1, l2, l3});
fluid_geo->transfiniteLine({l0}, 0.1 / h + 1);
fluid_geo->transfiniteLine({l1}, 0.08 / h + 1);
fluid_geo->transfiniteLine({l2}, 0.1 / h + 1);
fluid_geo->transfiniteLine({l3}, 0.08 / h + 1);
fluid_geo->addDirichletBoundaryCondition({l0}, Y, 0.0);
fluid_geo->addDirichletBoundaryCondition({l3}, X, 0.0);
fluid_geo->addInterfaceBoundaryCondition({l1});
Material *fluid_mat = new NewtonianFluid(1.0, 1000.0);
FluidDomain *fluid_problem = new FluidDomain(fluid_geo);
fluid_problem->applyMaterial({s0}, fluid_mat);
fluid_problem->generateMesh(T3, FRONT, "teste", "", true, false);
fluid_problem->setGravity(0.0, -9.81, 0.0);

//Solid problem
Geometry *solid_geo = new Geometry(1);
Point *p4 = solid_geo->addPoint({0.1   , 0.0  , 0.0});
Point *p5 = solid_geo->addPoint({0.112 , 0.0  , 0.0});
Point *p6 = solid_geo->addPoint({0.112 , 0.1  , 0.0});
Point *p7 = solid_geo->addPoint({0.1   , 0.1  , 0.0});
Point *p8 = solid_geo->addPoint({0.1   , 0.08 , 0.0});
Line *l4 = solid_geo->addLine({p4, p5});
Line *l5 = solid_geo->addLine({p5, p6});
Line *l6 = solid_geo->addLine({p6, p7});
Line *l7 = solid_geo->addLine({p7, p8});
Line *l8 = solid_geo->addLine({p8, p4});
Surface *s1 = solid_geo->addPlaneSurface({l4, l5, l6, l7, l8});
solid_geo->transfiniteLine({l4}, 0.012 / h + 1);
solid_geo->transfiniteLine({l5}, 0.1 / h + 1);
solid_geo->transfiniteLine({l6}, 0.012 / h + 1);
solid_geo->transfiniteLine({l7}, 0.02 / h + 1);
solid_geo->transfiniteLine({l8}, 0.08 / h + 1);
//solid_geo->transfiniteSurface({s1}, "Left", {p4, p5, p6, p7});
solid_geo->addDirichletBoundaryCondition({l4}, ALL, 0.0);
solid_geo->addInterfaceBoundaryCondition({l8});
Material *solid_mat = new ElasticSolid(8.0e6, 0.0, 2500.0);
SolidDomain *solid_problem = new SolidDomain(solid_geo);
solid_problem->applyMaterial({s1}, solid_mat);
solid_problem->generateMesh(T3, FRONT, "teste", "", true, false);
solid_problem->setGravity(0.0, 0.0, 0.0);

//Coupled problem
CoupledDomain *coupled_problem = new CoupledDomain(fluid_problem, solid_problem);
coupled_problem->setNumberOfSteps(1500);
coupled_problem->setDeltat(0.001);
coupled_problem->setMaxNonlinearIterations(6);
coupled_problem->setNonlinearTolerance(1.0e-4);
coupled_problem->setSpectralRadius(0.5);
coupled_problem->solveMonolithicCoupledProblem();