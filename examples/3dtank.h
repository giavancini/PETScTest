Geometry *tank_geo = new Geometry(0);

Point *p0 = tank_geo->addPoint({0.0, 0.0, 0.0});
Point *p1 = tank_geo->addPoint({1.0, 0.0, 0.0});
Point *p2 = tank_geo->addPoint({1.0, 1.0, 0.0});
Point *p3 = tank_geo->addPoint({0.0, 1.0, 0.0});
Point *p4 = tank_geo->addPoint({0.0, 0.0, -1.0});
Point *p5 = tank_geo->addPoint({1.0, 0.0, -1.0});
Point *p6 = tank_geo->addPoint({1.0, 1.0, -1.0});
Point *p7 = tank_geo->addPoint({0.0, 1.0, -1.0});

Line *l0 = tank_geo->addLine({p0, p1});
Line *l1 = tank_geo->addLine({p1, p2});
Line *l2 = tank_geo->addLine({p2, p3});
Line *l3 = tank_geo->addLine({p3, p0});
Line *l4 = tank_geo->addLine({p4, p5});
Line *l5 = tank_geo->addLine({p5, p6});
Line *l6 = tank_geo->addLine({p6, p7});
Line *l7 = tank_geo->addLine({p7, p4});
Line *l8 = tank_geo->addLine({p0, p4});
Line *l9 = tank_geo->addLine({p1, p5});
Line *l10 = tank_geo->addLine({p3, p7});
Line *l11 = tank_geo->addLine({p2, p6});

Surface *s0 = tank_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = tank_geo->addPlaneSurface({l4, l5, l6, l7});
Surface *s2 = tank_geo->addPlaneSurface({l8, -(*l7), -(*l10), l3});
Surface *s3 = tank_geo->addPlaneSurface({l9, l5, -(*l11), -(*l1)});
Surface *s4 = tank_geo->addPlaneSurface({l0, l9, -(*l4), -(*l8)});
Surface *s5 = tank_geo->addPlaneSurface({-(*l2), l11, l6, -(*l10)});

Volume *v0 = tank_geo->addVolume({s0, s1, s2, s3, s4, s5});

tank_geo->transfiniteLine({l0, l2, l4, l6, l8, l9, l10, l11}, 11);
tank_geo->transfiniteLine({l1, l3, l5, l7}, 11);
//tank_geo->transfiniteSurface({s0, s1, s2, s3, s4, s5});
//tank_geo->transfiniteVolume({v0});

tank_geo->addDirichletBoundaryCondition({s2}, X, 0.0);
tank_geo->addDirichletBoundaryCondition({s3}, X, 0.0);
tank_geo->addDirichletBoundaryCondition({s4}, Y, 0.0);
tank_geo->addDirichletBoundaryCondition({s0}, Z, 0.0);
tank_geo->addDirichletBoundaryCondition({s1}, Z, 0.0);

Material *mat = new NewtonianFluid(0.01, 1000.0);

FluidDomain *tank_problem = new FluidDomain(tank_geo);

tank_problem->applyMaterial({v0}, mat);
tank_problem->generateMesh(TET4, FRONT, "", "", true, false);
tank_problem->setNumberOfSteps(100);
tank_problem->setDeltat(0.01);
tank_problem->setMaxNonlinearIterations(6);
tank_problem->setNonlinearTolerance(1.0e-4);
tank_problem->setGravity(0.0, -10.0, 0.0);
//tank_problem->setLumpedMass(false);
//tank_problem->setSpectralRadius(0.9);
tank_problem->solveTransientProblem();