Geometry *solid_geo = new Geometry(0);

Point *p0 = solid_geo->addPoint({0.0, 0.0});
Point *p1 = solid_geo->addPoint({1.0, 0.0});
Point *p2 = solid_geo->addPoint({1.0, 1.0});
Point *p3 = solid_geo->addPoint({0.0, 1.0});

Line *l0 = solid_geo->addLine({p0, p1});
Line *l1 = solid_geo->addLine({p1, p2});
Line *l2 = solid_geo->addLine({p2, p3});
Line *l3 = solid_geo->addLine({p3, p0});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3});

solid_geo->transfiniteLine({l0}, 11);
solid_geo->transfiniteLine({l1}, 11);
solid_geo->transfiniteLine({l2}, 11);
solid_geo->transfiniteLine({l3}, 11);
//solid->transfiniteSurface({s0}, "Right", {p0, p1, p2, p3});

solid_geo->addDirichletBoundaryCondition({l3}, ALL, 0.0);
solid_geo->addNeumannBoundaryCondition({l1}, 1.0, 0.0, 0.0);

Material *mat = new ElasticSolid(1000.0, 0.0, 1.0);
mat->setPlaneAnalysis(PLANE_STRESS);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({s0}, mat);
solid_problem->generateMesh(T3, FRONT, "teste", "", true, false);
solid_problem->setNumberOfSteps(250);
solid_problem->setDeltat(8.e-4);
solid_problem->setSpectralRadius(1.0);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->setGravity(0.0, 0.0, 0.0);
solid_problem->solveTransientProblem();
