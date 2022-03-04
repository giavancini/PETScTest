Geometry *slosh_geo = new Geometry(0);

Point *p0 = slosh_geo->addPoint({-0.5, 0.0, 0.0});
Point *p1 = slosh_geo->addPoint({0.5, 0.0, 0.0});
Point *p2 = slosh_geo->addPoint({0.5, 1.01, 0.0});
Point *p3 = slosh_geo->addPoint({-0.5, 0.99, 0.0});
Point *p4 = slosh_geo->addPoint({-0.5, 0.0, -1.0});
Point *p5 = slosh_geo->addPoint({0.5, 0.0, -1.0});
Point *p6 = slosh_geo->addPoint({0.5, 1.01, -1.0});
Point *p7 = slosh_geo->addPoint({-0.5, 0.99, -1.0});

Line *l0 = slosh_geo->addLine({p0, p1});
Line *l1 = slosh_geo->addLine({p1, p2});
Line *l2 = slosh_geo->addSpline({p2, p3}, [](double x) { return (0.01 * std::sin(std::acos(-1.0) * x) + 1.0); }, 1000);
Line *l3 = slosh_geo->addLine({p3, p0});
Line *l4 = slosh_geo->addLine({p4, p5});
Line *l5 = slosh_geo->addLine({p5, p6});
Line *l6 = slosh_geo->addSpline({p6, p7}, [](double x) { return (0.01 * std::sin(std::acos(-1.0) * x) + 1.0); }, 1000);
Line *l7 = slosh_geo->addLine({p7, p4});
Line *l8 = slosh_geo->addLine({p0, p4});
Line *l9 = slosh_geo->addLine({p1, p5});
Line *l10 = slosh_geo->addLine({p3, p7});
Line *l11 = slosh_geo->addLine({p2, p6});

Surface *s0 = slosh_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = slosh_geo->addPlaneSurface({l4, l5, l6, l7});
Surface *s2 = slosh_geo->addPlaneSurface({l8, -(*l7), -(*l10), l3});
Surface *s3 = slosh_geo->addPlaneSurface({l9, l5, -(*l11), -(*l1)});
Surface *s4 = slosh_geo->addPlaneSurface({l0, l9, -(*l4), -(*l8)});
Surface *s5 = slosh_geo->addSurface({-(*l2), l11, l6, -(*l10)});

Volume *v0 = slosh_geo->addVolume({s0,s1,s2,s3,s4,s5});

slosh_geo->transfiniteLine({l0, l2, l4, l6, l8, l9, l10, l11}, 21);
slosh_geo->transfiniteLine({l1, l3, l5, l7}, 21);

slosh_geo->addDirichletBoundaryCondition({s2}, X, 0.0);
slosh_geo->addDirichletBoundaryCondition({s3}, X, 0.0);
slosh_geo->addDirichletBoundaryCondition({s4}, Y, 0.0);
slosh_geo->addDirichletBoundaryCondition({s0}, Z, 0.0);
slosh_geo->addDirichletBoundaryCondition({s1}, Z, 0.0);

Material *mat = new NewtonianFluid(0.01, 1.0);

FluidDomain *slosh_problem = new FluidDomain(slosh_geo);

slosh_problem->applyMaterial({v0}, mat);
slosh_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
slosh_problem->setNumberOfSteps(500);
slosh_problem->setDeltat(0.05);
slosh_problem->setMaxNonlinearIterations(6);
slosh_problem->setNonlinearTolerance(1.0e-6);
slosh_problem->setGravity(0.0, -1.0, 0.0);
slosh_problem->setSpectralRadius(0.9);
problem->setReferenceConfiguration(ReferenceConfiguration::PAST);
slosh_problem->solveTransientProblem();