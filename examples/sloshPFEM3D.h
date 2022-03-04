double h = 0.2;
double L = 1.0;

Geometry *slosh_geo = new Geometry(0);

Point *p0 = slosh_geo->addPoint({0.0        , 0.0,      0.0});
Point *p1 = slosh_geo->addPoint({10.0*L     , 0.0,      0.0});
Point *p2 = slosh_geo->addPoint({10.0*L     , 3.0*L,    0.0});
Point *p3 = slosh_geo->addPoint({0.0        , 7.0*L,    0.0});
Point *p4 = slosh_geo->addPoint({10.0*L     , 8.0*L,    0.0});
Point *p5 = slosh_geo->addPoint({0.0        , 8.0*L,    0.0});

Point *p6  = slosh_geo->addPoint({0.0        , 0.0,      -3*L});
Point *p7  = slosh_geo->addPoint({10.0*L     , 0.0,      -3*L});
Point *p8  = slosh_geo->addPoint({10.0*L     , 3.0*L,    -3*L});
Point *p9  = slosh_geo->addPoint({0.0        , 7.0*L,    -3*L});
Point *p10 = slosh_geo->addPoint({10.0*L     , 8.0*L,    -3*L});
Point *p11 = slosh_geo->addPoint({0.0        , 8.0*L,    -3*L});

Line *l0 = slosh_geo->addLine({p0, p1});
Line *l1 = slosh_geo->addLine({p1, p2});
Line *l2 = slosh_geo->addLine({p2, p3});
Line *l3 = slosh_geo->addLine({p3, p0});
Line *l4 = slosh_geo->addLine({p2, p4});
Line *l5 = slosh_geo->addLine({p4, p5});
Line *l6 = slosh_geo->addLine({p5, p3});

Line *l7  = slosh_geo->addLine({p6, p7});
Line *l8  = slosh_geo->addLine({p7, p8});
Line *l9  = slosh_geo->addLine({p8, p9});
Line *l10 = slosh_geo->addLine({p9, p6});
Line *l11 = slosh_geo->addLine({p8, p10});
Line *l12 = slosh_geo->addLine({p10, p11});
Line *l13 = slosh_geo->addLine({p11, p9});

Line *l14 = slosh_geo->addLine({p0, p6});
Line *l15 = slosh_geo->addLine({p1, p7});
Line *l16 = slosh_geo->addLine({p2, p8});
Line *l17 = slosh_geo->addLine({p3, p9});
Line *l18 = slosh_geo->addLine({p4, p10});
Line *l19 = slosh_geo->addLine({p5, p11});

Surface *s0 = slosh_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = slosh_geo->addPlaneSurface({l4, l5, l6, -(*l2)});
Surface *s2 = slosh_geo->addPlaneSurface({-(*l7), -(*l10), -(*l9), -(*l8)});
Surface *s3 = slosh_geo->addPlaneSurface({-(*l11), l9, -(*l13), -(*l12)});
Surface *s4 = slosh_geo->addPlaneSurface({l15, l8, -(*l16), -(*l1)});
Surface *s5 = slosh_geo->addPlaneSurface({l16, l11, -(*l18), -(*l4)});
Surface *s6 = slosh_geo->addPlaneSurface({-(*l14), -(*l3), l17, l10});
Surface *s7 = slosh_geo->addPlaneSurface({-(*l17), -(*l6), l19, l13});
Surface *s8 = slosh_geo->addPlaneSurface({-(*l0), l14, l7, -(*l15)});
Surface *s9 = slosh_geo->addPlaneSurface({-(*l2), l16, l9, -(*l17)});

Volume *v0 = slosh_geo->addVolume({s0, s2, s4, s6, s8, s9});

slosh_geo->transfiniteLine({l0}, 10*L / h + 1);
slosh_geo->transfiniteLine({l1}, 3*L / h + 1);
slosh_geo->transfiniteLine({l2}, sqrt(116*L*L) / h + 1);
slosh_geo->transfiniteLine({l3}, 7*L / h + 1);
slosh_geo->transfiniteLine({l4}, 5*L / h + 1);
slosh_geo->transfiniteLine({l5}, 10*L / h + 1);
slosh_geo->transfiniteLine({l6}, L / h + 1);

slosh_geo->transfiniteLine({l7}, 10*L / h + 1);
slosh_geo->transfiniteLine({l8}, 3*L / h + 1);
slosh_geo->transfiniteLine({l9}, sqrt(116*L*L) / h + 1);
slosh_geo->transfiniteLine({l10}, 7*L / h + 1);
slosh_geo->transfiniteLine({l11}, 5*L / h + 1);
slosh_geo->transfiniteLine({l12}, 10*L / h + 1);
slosh_geo->transfiniteLine({l13}, L / h + 1);

slosh_geo->transfiniteLine({l14}, 3*L / h + 1);
slosh_geo->transfiniteLine({l15}, 3*L / h + 1);
slosh_geo->transfiniteLine({l16}, 3*L / h + 1);
slosh_geo->transfiniteLine({l17}, 3*L / h + 1);
slosh_geo->transfiniteLine({l18}, 3*L / h + 1);
slosh_geo->transfiniteLine({l19}, 3*L / h + 1);

slosh_geo->addDirichletBoundaryCondition({s0}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s1}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s2}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s3}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s4}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s5}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s6}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s7}, ALL, 0.0);
slosh_geo->addDirichletBoundaryCondition({s8}, ALL, 0.0);

Material *mat = new NewtonianFluid(0.001, 1000.0);

FluidDomain *slosh_problem = new FluidDomain(slosh_geo);

slosh_problem->applyMaterial({v0}, mat);
slosh_problem->generateMesh(TET4, FRONT, "file", "", false, false);
slosh_problem->setNumberOfSteps(10000);
slosh_problem->setDeltat(0.001);
slosh_problem->setMaxNonlinearIterations(6);
slosh_problem->setNonlinearTolerance(1.0e-6);
slosh_problem->setGravity(0.0, -9.81, 0.0);
slosh_problem->setSpectralRadius(0.8);
slosh_problem->setMeshLength(h);
slosh_problem->solvePFEMProblem();