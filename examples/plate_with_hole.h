Geometry *solid_geo = new Geometry(0);

Point *p0 = solid_geo->addPoint({0.0, 0.0, 0.0}, 1.0, false);
Point *p1 = solid_geo->addPoint({10.0, 0.0, 0.0});
Point *p2 = solid_geo->addPoint({100.0, 0.0, 0.0});
Point *p3 = solid_geo->addPoint({100.0, 100.0, 0.0});
Point *p4 = solid_geo->addPoint({0.0, 100.0, 0.0});
Point *p5 = solid_geo->addPoint({0.0, 10.0, 0.0});

Line *l0 = solid_geo->addLine({p1, p2});
Line *l1 = solid_geo->addLine({p2, p3});
Line *l2 = solid_geo->addLine({p3, p4});
Line *l3 = solid_geo->addLine({p4, p5});
Line *l4 = solid_geo->addCircle({p5, p0, p1});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3, l4});

solid_geo->transfiniteLine({l0}, 41);
solid_geo->transfiniteLine({l1}, 51);
solid_geo->transfiniteLine({l2}, 51);
solid_geo->transfiniteLine({l3}, 41);
solid_geo->transfiniteLine({l4}, 21);
//solid->transfiniteSurface({s0}, "Right", {p0, p1, p2, p4});

solid_geo->addDirichletBoundaryCondition({l0}, Y, 0.0);
solid_geo->addDirichletBoundaryCondition({l3}, X, 0.0);
solid_geo->addNeumannBoundaryCondition({l1}, 10.0, 0.0, 0.0);

Material *mat = new ElasticSolid(1000000000.0, 0.0);
mat->setPlaneAnalysis(PLANE_STRESS);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({s0}, mat);
solid_problem->generateMesh(T10, FRONT, "teste", "", true, false);
solid_problem->setNumberOfSteps(10);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->solveStaticProblem();