Geometry *dam_geo = new Geometry(0);

Point *p0 = dam_geo->addPoint({0.0, 0.0, 0.0});
Point *p1 = dam_geo->addPoint({0.35, 0.0, 0.0});
Point *p2 = dam_geo->addPoint({0.35, 0.7, 0.0});
Point *p3 = dam_geo->addPoint({0.0, 0.7, 0.0});
Point *p4 = dam_geo->addPoint({0.0, 0.0, -0.35});
Point *p5 = dam_geo->addPoint({0.35, 0.0, -0.35});
Point *p6 = dam_geo->addPoint({0.35, 0.7, -0.35});
Point *p7 = dam_geo->addPoint({0.0, 0.7, -0.35});

Line *l0 = dam_geo->addLine({p0, p1});
Line *l1 = dam_geo->addLine({p1, p2});
Line *l2 = dam_geo->addLine({p2, p3});
Line *l3 = dam_geo->addLine({p3, p0});
Line *l4 = dam_geo->addLine({p4, p5});
Line *l5 = dam_geo->addLine({p5, p6});
Line *l6 = dam_geo->addLine({p6, p7});
Line *l7 = dam_geo->addLine({p7, p4});
Line *l8 = dam_geo->addLine({p0, p4});
Line *l9 = dam_geo->addLine({p1, p5});
Line *l10 = dam_geo->addLine({p3, p7});
Line *l11 = dam_geo->addLine({p2, p6});

Surface *s0 = dam_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = dam_geo->addPlaneSurface({l4, l5, l6, l7});
Surface *s2 = dam_geo->addPlaneSurface({l8, -(*l7), -(*l10), l3});
Surface *s3 = dam_geo->addPlaneSurface({l9, l5, -(*l11), -(*l1)});
Surface *s4 = dam_geo->addPlaneSurface({l0, l9, -(*l4), -(*l8)});
Surface *s5 = dam_geo->addPlaneSurface({-(*l2), l11, l6, -(*l10)});

Volume *v0 = dam_geo->addVolume({s0,s1,s2,s3,s4,s5});

dam_geo->transfiniteLine({l0, l2, l4, l6, l8, l9, l10, l11}, 21);
dam_geo->transfiniteLine({l1, l3, l5, l7}, 41);

dam_geo->addDirichletBoundaryCondition({s2}, X, 0.0);
dam_geo->addDirichletBoundaryCondition({s4}, Y, 0.0);
dam_geo->addDirichletBoundaryCondition({s0}, Z, 0.0);
dam_geo->addDirichletBoundaryCondition({s1}, Z, 0.0);

Material *mat = new NewtonianFluid(0.001, 1.0);

FluidDomain *dam_problem = new FluidDomain(dam_geo);

dam_problem->applyMaterial({v0}, mat);
dam_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
dam_problem->setNumberOfSteps(335);
dam_problem->setDeltat(0.005);
dam_problem->setMaxNonlinearIterations(6);
dam_problem->setNonlinearTolerance(1.0e-4);
dam_problem->setGravity(0.0, -1.0, 0.0);
dam_problem->setSpectralRadius(0.5);
dam_problem->solveTransientProblem();