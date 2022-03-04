Geometry *solid_geo = new Geometry(0);

double L = 25.0;
double D = 4.0;
double w = 1.0;

Point *p0 = solid_geo->addPoint({0.0, 0.0, 0.0});
Point *p1 = solid_geo->addPoint({L  , 0.0, 0.0});
Point *p2 = solid_geo->addPoint({L  , D  , 0.0});
Point *p3 = solid_geo->addPoint({0.0, D  , 0.0});
Point *p4 = solid_geo->addPoint({0.0, 0.0, -w });
Point *p5 = solid_geo->addPoint({L  , 0.0, -w });
Point *p6 = solid_geo->addPoint({L  , D  , -w });
Point *p7 = solid_geo->addPoint({0.0, D  , -w });

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
Line *l10 = solid_geo->addLine({p2, p6});
Line *l11 = solid_geo->addLine({p3, p7});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = solid_geo->addPlaneSurface({l4, l5, l6, l7});
Surface *s2 = solid_geo->addPlaneSurface({l9, l5, -(*l10), -(*l1)});
Surface *s3 = solid_geo->addPlaneSurface({l8, -(*l7), -(*l11), l3});
Surface *s4 = solid_geo->addPlaneSurface({l0, l9, -(*l4), -(*l8)});
Surface *s5 = solid_geo->addPlaneSurface({-(*l2), l10, l6, -(*l11)});

Volume *v0 = solid_geo->addVolume({s0, s1, s2, s3, s4, s5});

double h = 0.25;
solid_geo->transfiniteLine({l0}, L/h + 1);
solid_geo->transfiniteLine({l1}, D/h + 1);
solid_geo->transfiniteLine({l2}, L/h + 1);
solid_geo->transfiniteLine({l3}, D/h + 1);
solid_geo->transfiniteLine({l4}, L/h + 1);
solid_geo->transfiniteLine({l5}, D/h + 1);
solid_geo->transfiniteLine({l6}, L/h + 1);
solid_geo->transfiniteLine({l7}, D/h + 1);
solid_geo->transfiniteLine({l8}, w/h + 1);
solid_geo->transfiniteLine({l9}, w/h + 1);
solid_geo->transfiniteLine({l10}, w/h + 1);
solid_geo->transfiniteLine({l11}, w/h + 1);
//solid_geo->transfiniteSurface({s0}, "Right", {p0, p1, p2, p3});

solid_geo->addDirichletBoundaryCondition({s3}, ALL, 0.0);
solid_geo->addDirichletBoundaryCondition({v0}, Z, 0.0);
solid_geo->addNeumannBoundaryCondition({s2}, 0.0, 10.0, 0.0);

Material *mat = new ElasticSolid(1.0e4, 0.25, 1.0);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({v0}, mat);
solid_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
solid_problem->setNumberOfSteps(1600);
solid_problem->setDeltat(0.01);
//solid_problem->setSpectralRadius(0.8);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->setGravity(0.0, 0.0, 0.0);
solid_problem->solveTransientProblem();