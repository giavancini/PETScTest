double h = 0.012/2.0;
double L = 0.146;

Geometry *dam_geo = new Geometry(0);

Point *p0 = dam_geo->addPoint({0.0  , 0.0  , 0.0});
Point *p1 = dam_geo->addPoint({L    , 0.0  , 0.0});
Point *p2 = dam_geo->addPoint({L    , 2.0*L, 0.0});
Point *p3 = dam_geo->addPoint({0.0  , 2.0*L, 0.0});
Point *p4 = dam_geo->addPoint({0.0  , 4.0*L, 0.0});
Point *p5 = dam_geo->addPoint({4.0*L, 0.0  , 0.0});
Point *p6 = dam_geo->addPoint({4.0*L, 4.0*L, 0.0});

Point *p7 = dam_geo->addPoint({0.0  , 0.0  ,  L});
Point *p8 = dam_geo->addPoint({L    , 0.0  ,  L});
Point *p9 = dam_geo->addPoint({L    , 2.0*L,  L});
Point *p10 = dam_geo->addPoint({0.0  , 2.0*L, L});
Point *p11 = dam_geo->addPoint({0.0  , 4.0*L, L});
Point *p12 = dam_geo->addPoint({4.0*L, 0.0  , L});
Point *p13 = dam_geo->addPoint({4.0*L, 4.0*L, L});

Line *l0 = dam_geo->addLine({p0, p1});
Line *l1 = dam_geo->addLine({p1, p2});
Line *l2 = dam_geo->addLine({p2, p3});
Line *l3 = dam_geo->addLine({p3, p0});
Line *l4 = dam_geo->addLine({p1, p5});
Line *l5 = dam_geo->addLine({p5, p6});
Line *l6 = dam_geo->addLine({p3, p4});
Line *l7 = dam_geo->addLine({p6, p4});

Line *l8 = dam_geo->addLine({p7, p8});
Line *l9 = dam_geo->addLine({p8, p9});
Line *l10 = dam_geo->addLine({p9, p10});
Line *l11 = dam_geo->addLine({p10, p7});
Line *l12 = dam_geo->addLine({p8, p12});
Line *l13 = dam_geo->addLine({p12, p13});
Line *l14 = dam_geo->addLine({p10, p11});
Line *l15 = dam_geo->addLine({p13, p11});

Line *l16 = dam_geo->addLine({p0, p7});
Line *l17 = dam_geo->addLine({p1, p8});
Line *l18 = dam_geo->addLine({p2, p9});
Line *l19 = dam_geo->addLine({p3, p10});
Line *l20 = dam_geo->addLine({p4, p11});
Line *l21 = dam_geo->addLine({p5, p12});
Line *l22 = dam_geo->addLine({p6, p13});

Surface *s0 = dam_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = dam_geo->addPlaneSurface({l8, l9, l10, l11});
Surface *s2 = dam_geo->addPlaneSurface({l4, l5, l7, -(*l6), -(*l2), -(*l1)});
Surface *s3 = dam_geo->addPlaneSurface({-(*l12), l9, l10, l14, -(*l15), -(*l13)});
Surface *s4 = dam_geo->addPlaneSurface({l21, l13, -(*l22), -(*l5)});
Surface *s5 = dam_geo->addPlaneSurface({l17, l9, -(*l18), -(*l1)});
Surface *s6 = dam_geo->addPlaneSurface({-(*l19), l6, l20, -(*l14)});
Surface *s7 = dam_geo->addPlaneSurface({-(*l16), -(*l3), l19, l11});
Surface *s8 = dam_geo->addPlaneSurface({-(*l0), l16, l8, -(*l17)});
Surface *s9 = dam_geo->addPlaneSurface({-(*l4), l17, l12, -(*l21)});
Surface *s10 = dam_geo->addPlaneSurface({-(*l2), l18, l10, -(*l19)});

Volume *v0 = dam_geo->addVolume({s0, s1, s5, s7, s8, s10});

dam_geo->transfiniteLine({l0}, L/h+1);
dam_geo->transfiniteLine({l1}, 2*L/h+1);
dam_geo->transfiniteLine({l2}, L/h+1);
dam_geo->transfiniteLine({l3}, 2*L/h+1);
dam_geo->transfiniteLine({l4}, 3*L/h+1);
dam_geo->transfiniteLine({l5}, 4*L/h+1);
dam_geo->transfiniteLine({l6}, 2*L/h+1);
dam_geo->transfiniteLine({l7}, 4*L/h+1);

dam_geo->transfiniteLine({l8}, L/h+1);
dam_geo->transfiniteLine({l9}, 2*L/h+1);
dam_geo->transfiniteLine({l10}, L/h+1);
dam_geo->transfiniteLine({l11}, 2*L/h+1);
dam_geo->transfiniteLine({l12}, 3*L/h+1);
dam_geo->transfiniteLine({l13}, 4*L/h+1);
dam_geo->transfiniteLine({l14}, 2*L/h+1);
dam_geo->transfiniteLine({l15}, 4*L/h+1);

dam_geo->transfiniteLine({l16}, L/h+1);
dam_geo->transfiniteLine({l17}, L/h+1);
dam_geo->transfiniteLine({l18}, L/h+1);
dam_geo->transfiniteLine({l19}, L/h+1);
dam_geo->transfiniteLine({l20}, L/h+1);
dam_geo->transfiniteLine({l21}, L/h+1);
dam_geo->transfiniteLine({l22}, L/h+1);

dam_geo->addDirichletBoundaryCondition({s0}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s1}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s2}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s3}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s4}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s6}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s7}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s8}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({s9}, ALL, 0.0);

Material *mat = new NewtonianFluid(0.001, 1000.0);

FluidDomain *dam_problem = new FluidDomain(dam_geo);

dam_problem->applyMaterial({v0}, mat);
dam_problem->generateMesh(TET4, FRONT, "file", "", false, false);
dam_problem->setNumberOfSteps(5000);
dam_problem->setDeltat(0.0002);
dam_problem->setMaxNonlinearIterations(6);
dam_problem->setNonlinearTolerance(1.0e-4);
dam_problem->setGravity(0.0, -9.81, 0.0);
dam_problem->setSpectralRadius(0.8);
dam_problem->setMeshLength(h);
dam_problem->setAlpha(1.35);
dam_problem->setExportFrequency(1);
dam_problem->solvePFEMProblem();
