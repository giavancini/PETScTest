double h = 0.012;
double D = 0.5;
double H1 = 0.35;
double H2 = 0.15;

Geometry *slosh_geo = new Geometry(0);

Point *p0 = slosh_geo->addPoint({0.0 , 0.0});
Point *p1 = slosh_geo->addPoint({D   , 0.0});
Point *p2 = slosh_geo->addPoint({D   , H2 });
Point *p3 = slosh_geo->addPoint({0.0 , H1 });
Point *p4 = slosh_geo->addPoint({0.0 , D});
Point *p5 = slosh_geo->addPoint({D   , D});

Line *l0 = slosh_geo->addLine({p0, p1});
Line *l1 = slosh_geo->addLine({p1, p2});
Line *l2 = slosh_geo->addLine({p2, p3});
Line *l3 = slosh_geo->addLine({p3, p0});
Line *l4 = slosh_geo->addLine({p2, p5});
Line *l5 = slosh_geo->addLine({p3, p4});

Surface *s0 = slosh_geo->addPlaneSurface({l0, l1, l2, l3});

slosh_geo->transfiniteLine({l0}, D / h + 1);
slosh_geo->transfiniteLine({l1}, H2 / h + 1);
slosh_geo->transfiniteLine({l2}, sqrt(D*D + (H1-H2)*(H1-H2)) / h + 1);
slosh_geo->transfiniteLine({l3}, H1 / h + 1);
slosh_geo->transfiniteLine({l4}, (D-H2) / h + 1);
slosh_geo->transfiniteLine({l5}, (D-H1) / h + 1);
//slosh_geo->transfiniteSurface({s0});

slosh_geo->addDirichletBoundaryCondition({l0}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({l1}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({l3}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({l4}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({l5}, ALL, 0.0);

Material *mat = new NewtonianFluid(5.0, 1000.0);

FluidDomain *slosh_problem = new FluidDomain(slosh_geo);

slosh_problem->applyMaterial({s0}, mat);
slosh_problem->generateMesh(T3, FRONT, "teste", "", true, false);
slosh_problem->setNumberOfSteps(10000);
slosh_problem->setDeltat(0.001);
slosh_problem->setMaxNonlinearIterations(6);
slosh_problem->setNonlinearTolerance(1.0e-6);
slosh_problem->setGravity(0.0, -9.81, 0.0);
slosh_problem->setSpectralRadius(0.8);
slosh_problem->setMeshLength(h);
slosh_problem->solvePFEMProblem();