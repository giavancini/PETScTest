Geometry *dam_geo = new Geometry(0);

double L = 0.5;
double h = 0.075;

Point *p0 = dam_geo->addPoint({0.0    , 0.0  });
Point *p1 = dam_geo->addPoint({4.0*L  , 0.0  });
Point *p2 = dam_geo->addPoint({4.0*L  , 5.0*L});
Point *p3 = dam_geo->addPoint({0.0    , 5.0*L});
Point *p4 = dam_geo->addPoint({0.0    , 8.0*L});
Point *p5 = dam_geo->addPoint({16.0*L , 0.0  });
Point *p6 = dam_geo->addPoint({20.0*L , 0.0  });
Point *p7 = dam_geo->addPoint({20.0*L , 5.0*L});
Point *p8 = dam_geo->addPoint({16.0*L , 5.0*L});
Point *p9 = dam_geo->addPoint({20.0*L , 8.0*L});

Line *l0 = dam_geo->addLine({p0, p1});
Line *l1 = dam_geo->addLine({p1, p2});
Line *l2 = dam_geo->addLine({p2, p3});
Line *l3 = dam_geo->addLine({p3, p0});
Line *l4 = dam_geo->addLine({p5, p6});
Line *l5 = dam_geo->addLine({p6, p7});
Line *l6 = dam_geo->addLine({p7, p8});
Line *l7 = dam_geo->addLine({p8, p5});
Line *l8 = dam_geo->addLine({p1, p5});
Line *l9 = dam_geo->addLine({p3, p4});
Line *l10 = dam_geo->addLine({p7, p9});

Surface *s0 = dam_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = dam_geo->addPlaneSurface({l4, l5, l6, l7});

dam_geo->transfiniteLine({l0}, 4.0*L/h+1);
dam_geo->transfiniteLine({l1}, 5.0*L/h+1);
dam_geo->transfiniteLine({l2}, 4.0*L/h+1);
dam_geo->transfiniteLine({l3}, 5.0*L/h+1);
dam_geo->transfiniteLine({l4}, 4.0*L/h+1);
dam_geo->transfiniteLine({l5}, 5.0*L/h+1);
dam_geo->transfiniteLine({l6}, 4.0*L/h+1);
dam_geo->transfiniteLine({l7}, 5.0*L/h+1);
dam_geo->transfiniteLine({l8}, 12.0*L/h+1);
dam_geo->transfiniteLine({l9}, 3.0*L/h+1);
dam_geo->transfiniteLine({l10}, 3.0*L/h+1);

dam_geo->addDirichletBoundaryCondition({l0}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l3}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l4}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l5}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l8}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l9}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l10}, ALL, 0.0);

Material *mat = new NewtonianFluid(0.001, 1000.0);

FluidDomain *dam_problem = new FluidDomain(dam_geo);

dam_problem->applyMaterial({s0}, mat);
dam_problem->applyMaterial({s1}, mat);
dam_problem->generateMesh(T3, FRONT, "teste", "", true, false);
dam_problem->setNumberOfSteps(16000);
dam_problem->setDeltat(0.0005);
dam_problem->setMaxNonlinearIterations(6);
dam_problem->setNonlinearTolerance(1.0e-6);
dam_problem->setGravity(0.0, -9.81, 0.0);
dam_problem->setSpectralRadius(0.8);
dam_problem->setMeshLength(h);
dam_problem->solvePFEMProblem();