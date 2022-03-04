Geometry *solid_geo = new Geometry(0);

Point *p0 = solid_geo->addPoint({0.0, 0.0, 0.0});
Point *p1 = solid_geo->addPoint({1.0, 0.0, 0.0});
Point *p2 = solid_geo->addPoint({1.0, 1.0, 0.0});
Point *p3 = solid_geo->addPoint({0.0, 1.0, 0.0});
Point *p4 = solid_geo->addPoint({0.0, 0.0, -1.0});
Point *p5 = solid_geo->addPoint({1.0, 0.0, -1.0});
Point *p6 = solid_geo->addPoint({1.0, 1.0, -1.0});
Point *p7 = solid_geo->addPoint({0.0, 1.0, -1.0});

Line *l0 = solid_geo->addLine({p0, p1});
Line *l1 = solid_geo->addLine({p1, p2});
Line *l2 = solid_geo->addLine({p2, p3});
Line *l3 = solid_geo->addLine({p3, p0});
Line *l4 = solid_geo->addLine({p4, p5});
Line *l5 = solid_geo->addLine({p5, p6});
Line *l6 = solid_geo->addLine({p6, p7});
Line *l7 = solid_geo->addLine({p7, p4});
Line *l8 = solid_geo->addLine({p0, p4});
Line *l9 = solid_geo->addLine({p1, p5});
Line *l10 = solid_geo->addLine({p3, p7});
Line *l11 = solid_geo->addLine({p2, p6});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = solid_geo->addPlaneSurface({l4, l5, l6, l7});
Surface *s2 = solid_geo->addPlaneSurface({l8, -(*l7), -(*l10), l3});
Surface *s3 = solid_geo->addPlaneSurface({l9, l5, -(*l11), -(*l1)});
Surface *s4 = solid_geo->addPlaneSurface({l0, l9, -(*l4), -(*l8)});
Surface *s5 = solid_geo->addPlaneSurface({-(*l2), l11, l6, -(*l10)});

Volume *v0 = solid_geo->addVolume({s0,s1,s2,s3,s4,s5});

solid_geo->transfiniteLine({l0,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11}, 11);

solid_geo->addDirichletBoundaryCondition({s2}, ALL, 0.0);
solid_geo->addNeumannBoundaryCondition({s3}, 100000.0, 0.0, 0.0);

Material *mat = new ElasticSolid(210.e9, 0.0);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({v0}, mat);
solid_problem->generateMesh(TET4, DELAUNAY, "teste", "", true, false);
solid_problem->setNumberOfSteps(1);
solid_problem->setGravity(0.0, 0.0, 0.0);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->solveStaticProblem();