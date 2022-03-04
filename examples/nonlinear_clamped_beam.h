Geometry *solid_geo = new Geometry(0);

Point *p0 = solid_geo->addPoint({0.0 , 0.0 , 0.0});
Point *p1 = solid_geo->addPoint({20.0, 0.0 , 0.0});
Point *p2 = solid_geo->addPoint({20.0, 0.125, 0.0});
Point *p3 = solid_geo->addPoint({10.0, 0.125, 0.0});
Point *p4 = solid_geo->addPoint({0.0 , 0.125, 0.0});

Line *l0 = solid_geo->addLine({p0, p1});
Line *l1 = solid_geo->addLine({p1, p2});
Line *l2 = solid_geo->addLine({p2, p3});
Line *l3 = solid_geo->addLine({p3, p4});
Line *l4 = solid_geo->addLine({p4, p0});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3, l4});

int n = 3*160;
solid_geo->transfiniteLine({l0}, n + 1);
solid_geo->transfiniteLine({l1}, n/160 + 1);
solid_geo->transfiniteLine({l2}, n/2 + 1);
solid_geo->transfiniteLine({l3}, n/2 + 1);
solid_geo->transfiniteLine({l4}, n/160 + 1);

// solid_geo->transfiniteLine({l0}, 201);
// solid_geo->transfiniteLine({l1}, 3);
// solid_geo->transfiniteLine({l2}, 101);
// solid_geo->transfiniteLine({l3}, 101);
// solid_geo->transfiniteLine({l4}, 3);
solid_geo->transfiniteSurface({s0}, "Right", {p0, p1, p2, p4});

solid_geo->addDirichletBoundaryCondition({l1}, ALL, 0.0);
solid_geo->addDirichletBoundaryCondition({l4}, ALL, 0.0);
solid_geo->addNeumannBoundaryCondition({p3}, 0.0, 640, 0.0);

Material* mat = new ElasticSolid(30000.0e3, 0.0, 2.537e-4);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({s0}, mat);
solid_problem->generateMesh(T3, FRONT, "teste", "", false, true);
solid_problem->setNumberOfSteps(2000);
solid_problem->setDeltat(2.50e-6);
solid_problem->setSpectralRadius(0.5);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->setGravity(0.0, 0.0, 0.0);
solid_problem->solveTransientProblem();