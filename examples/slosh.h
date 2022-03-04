Geometry *slosh_geo = new Geometry(0);

Point *p0 = slosh_geo->addPoint({-0.5, 0.0, 0.0});
Point *p1 = slosh_geo->addPoint({0.5, 0.0, 0.0});
Point *p2 = slosh_geo->addPoint({0.5, 1.01, 0.0});
Point *p3 = slosh_geo->addPoint({-0.5, 0.99, 0.0});

Line *l0 = slosh_geo->addLine({p0, p1});
Line *l1 = slosh_geo->addLine({p1, p2});
Line *l2 = slosh_geo->addSpline({p2, p3}, [](double x) { return (0.01 * std::sin(std::acos(-1.0) * x) + 1.0); }, 1000);
Line *l3 = slosh_geo->addLine({p3, p0});

Surface *s0 = slosh_geo->addPlaneSurface({l0, l1, l2, l3});

slosh_geo->transfiniteLine({l0}, 21);
slosh_geo->transfiniteLine({l1}, 21);
slosh_geo->transfiniteLine({l2}, 21);
slosh_geo->transfiniteLine({l3}, 21);
slosh_geo->transfiniteSurface({s0}, "RIGHT", {p0, p1, p2, p3});

slosh_geo->addDirichletBoundaryCondition({l0}, Y, 0.0);
slosh_geo->addDirichletBoundaryCondition({l1}, X, 0.0);
slosh_geo->addDirichletBoundaryCondition({l3}, X, 0.0);

Material *mat = new NewtonianFluid(0.01, 1.0);

FluidDomain *slosh_problem = new FluidDomain(slosh_geo);

slosh_problem->applyMaterial({s0}, mat);
slosh_problem->generateMesh(T3, FRONT, "teste", "", true, false);
slosh_problem->setNumberOfSteps(500);
slosh_problem->setDeltat(0.05);
slosh_problem->setMaxNonlinearIterations(3);
slosh_problem->setNonlinearTolerance(1.0e-6);
slosh_problem->setGravity(0.0, -1.0, 0.0);
slosh_problem->setSpectralRadius(0.9);
slosh_problem->setReferenceConfiguration(ReferenceConfiguration::PAST);
slosh_problem->solveTransientProblem();