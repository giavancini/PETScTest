Geometry *solid_geo = new Geometry(0);

Point *p0 = solid_geo->addPoint({0.0  , 0.0   , 0.0});
Point *p1 = solid_geo->addPoint({20.0 , 0.0   , 0.0});
Point *p2 = solid_geo->addPoint({20.0 , 0.125 , 0.0});
Point *p3 = solid_geo->addPoint({10.0 , 0.125 , 0.0});
Point *p4 = solid_geo->addPoint({0.0  , 0.125 , 0.0});
Point *p5 = solid_geo->addPoint({0.0  , 0.0   , -1.0});
Point *p6 = solid_geo->addPoint({20.0 , 0.0   , -1.0});
Point *p7 = solid_geo->addPoint({20.0 , 0.125 , -1.0});
Point *p8 = solid_geo->addPoint({10.0 , 0.125 , -1.0});
Point *p9 = solid_geo->addPoint({0.0  , 0.125 , -1.0});

Line *l0 = solid_geo->addLine({p0, p1});
Line *l1 = solid_geo->addLine({p1, p2});
Line *l2 = solid_geo->addLine({p2, p3});
Line *l3 = solid_geo->addLine({p3, p4});
Line *l4 = solid_geo->addLine({p4, p0});
Line *l5 = solid_geo->addLine({p5, p6});
Line *l6 = solid_geo->addLine({p6, p7});
Line *l7 = solid_geo->addLine({p7, p8});
Line *l8 = solid_geo->addLine({p8, p9});
Line *l9 = solid_geo->addLine({p9, p5});
Line *l10 = solid_geo->addLine({p0, p5});
Line *l11 = solid_geo->addLine({p1, p6});
Line *l12 = solid_geo->addLine({p2, p7});
Line *l13 = solid_geo->addLine({p3, p8});
Line *l14 = solid_geo->addLine({p4, p9});

Surface *s0 = solid_geo->addPlaneSurface({l0, l1, l2, l3, l4});
Surface *s1 = solid_geo->addPlaneSurface({l5, l6, l7, l8, l9});
Surface *s2 = solid_geo->addPlaneSurface({l11, l6, -(*l12), -(*l1)});
Surface *s3 = solid_geo->addPlaneSurface({l10, -(*l9), -(*l14), l4});
Surface *s4 = solid_geo->addPlaneSurface({l0, l11, -(*l5), -(*l10)});
Surface *s5 = solid_geo->addPlaneSurface({-(*l2), l12, l7, -(*l13)});
Surface *s6 = solid_geo->addPlaneSurface({-(*l3), l13, l8, -(*l14)});

Volume *v0 = solid_geo->addVolume({s0, s1, s2, s3, s4, s5, s6});

int n = 3*160;
solid_geo->transfiniteLine({l0}, n + 1);
solid_geo->transfiniteLine({l1}, n/160 + 1);
solid_geo->transfiniteLine({l2}, n/2 + 1);
solid_geo->transfiniteLine({l3}, n/2 + 1);
solid_geo->transfiniteLine({l4}, n/160 + 1);
solid_geo->transfiniteLine({l5}, n + 1);
solid_geo->transfiniteLine({l6}, n/160 + 1);
solid_geo->transfiniteLine({l7}, n/2 + 1);
solid_geo->transfiniteLine({l8}, n/2 + 1);
solid_geo->transfiniteLine({l9}, n/160 + 1);
solid_geo->transfiniteLine({l10}, n/20 + 1);
solid_geo->transfiniteLine({l11}, n/20 + 1);
solid_geo->transfiniteLine({l12}, n/20 + 1);
solid_geo->transfiniteLine({l13}, n/20 + 1);
solid_geo->transfiniteLine({l14}, n/20 + 1);
//solid_geo->transfiniteSurface({s0}, "Right", {p0, p1, p2, p4});

solid_geo->addDirichletBoundaryCondition({s2}, ALL, 0.0);
solid_geo->addDirichletBoundaryCondition({s3}, ALL, 0.0);
solid_geo->addNeumannBoundaryCondition({l13}, 0.0, 640, 0.0);

Material* mat = new ElasticSolid(30000.0e3, 0.0, 2.537e-4);

SolidDomain *solid_problem = new SolidDomain(solid_geo);

solid_problem->applyMaterial({v0}, mat);
solid_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
solid_problem->setNumberOfSteps(2000);
solid_problem->setDeltat(2.50e-6);
solid_problem->setSpectralRadius(0.5);
solid_problem->setMaxNonlinearIterations(6);
solid_problem->setNonlinearTolerance(1.0e-6);
solid_problem->setGravity(0.0, 0.0, 0.0);
solid_problem->solveTransientProblem();