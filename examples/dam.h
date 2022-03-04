Geometry *dam_geo = new Geometry(0);

Point *p0 = dam_geo->addPoint({0.0, 0.0, 0.0}, 0.1);
Point *p1 = dam_geo->addPoint({0.35, 0.0, 0.0}, 0.001);
Point *p2 = dam_geo->addPoint({0.35, 0.7, 0.0}, 0.1);
Point *p3 = dam_geo->addPoint({0.0, 0.7, 0.0}, 0.1);

Line *l0 = dam_geo->addLine({p0, p1});
Line *l1 = dam_geo->addLine({p1, p2});
Line *l2 = dam_geo->addLine({p2, p3});
Line *l3 = dam_geo->addLine({p3, p0});

Surface *s0 = dam_geo->addPlaneSurface({l0, l1, l2, l3});

dam_geo->transfiniteLine({l0}, 11);
dam_geo->transfiniteLine({l1}, 21);
dam_geo->transfiniteLine({l2}, 11);
dam_geo->transfiniteLine({l3}, 21);
//dam_geo->transfiniteSurface({ s0 }, "Right", {p0,p1,p2,p3});

dam_geo->addDirichletBoundaryCondition({l0}, Y, 0.0);
dam_geo->addDirichletBoundaryCondition({l3}, X, 0.0);
//dam_geo->addDirichletBoundaryCondition({l1}, X, 0.0);

Material *mat = new NewtonianFluid(0.001, 1.0);

FluidDomain *dam_problem = new FluidDomain(dam_geo);

dam_problem->applyMaterial({s0}, mat);
dam_problem->generateMesh(T3, FRONT, "teste", "", true, false);
dam_problem->setNumberOfSteps(335);
dam_problem->setDeltat(0.005);
dam_problem->setMaxNonlinearIterations(6);
dam_problem->setNonlinearTolerance(1.0e-4);
dam_problem->setGravity(0.0, -1.0, 0.0);
//dam_problem->setLumpedMass(false);
dam_problem->setSpectralRadius(0.5);
dam_problem->solveTransientProblem();