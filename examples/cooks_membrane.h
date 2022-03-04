Geometry *solid_geo = new Geometry(0);

Point *p0 = solid_geo->addPoint({0.0  , 0.0  , 0.0});
Point *p1 = solid_geo->addPoint({48.0 , 44.0 , 0.0});
Point *p2 = solid_geo->addPoint({48.0 , 60.0 , 0.0});
Point *p3 = solid_geo->addPoint({0.0  , 44.0 , 0.0});

Line *l0 = solid_geo->addLine({p0, p1});
Line *l1 = solid_geo->addLine({p1, p2});
Line *l2 = solid_geo->addLine({p2, p3});
Line *l3 = solid_geo->addLine({p3, p0});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3});

int n = 10;
solid_geo->transfiniteLine({l0}, n+1);
solid_geo->transfiniteLine({l1}, n+1);
solid_geo->transfiniteLine({l2}, n+1);
solid_geo->transfiniteLine({l3}, n+1);
solid_geo->transfiniteSurface({s0}, "Right", {p0, p1, p2, p3});

solid_geo->addDirichletBoundaryCondition({l3}, ALL, 0.0);
solid_geo->addNeumannBoundaryCondition({l1}, 0.0, 1.0/16.0, 0.0);

Material *mat = new ElasticSolid(1.0, 0.33333333, 0.0);
mat->setPlaneAnalysis(PLANE_STRESS);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({s0}, mat);
solid_problem->generateMesh(T10, FRONT, "teste", "", true, false);
solid_problem->setNumberOfSteps(100);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->solveStaticProblem();