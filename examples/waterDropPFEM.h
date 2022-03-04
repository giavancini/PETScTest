Geometry *dam_geo = new Geometry(0);

Point *p0 = dam_geo->addPoint({-0.5,  0.0  });
Point *p1 = dam_geo->addPoint({0.5 ,  0.0  });
Point *p2 = dam_geo->addPoint({0.0 , -0.8  });
Point *p3 = dam_geo->addPoint({-5.0, -4.0});
Point *p4 = dam_geo->addPoint({5.0 , -4.0});

Line *l0 = dam_geo->addLine({p0, p2});
Line *l1 = dam_geo->addLine({p2, p1});
Line *l2 = dam_geo->addLine({p1, p0});
Line *l3 = dam_geo->addLine({p3, p4});

Surface *s0 = dam_geo->addPlaneSurface({l0, l1, l2});

int h = 20;
dam_geo->transfiniteLine({l0}, h+1);
dam_geo->transfiniteLine({l1}, h+1);
dam_geo->transfiniteLine({l2}, h+1);
dam_geo->transfiniteLine({l3}, 10*h+1);
//dam_geo->transfiniteSurface({s0});

dam_geo->addDirichletBoundaryCondition({l2}, ALL, 0.0);
dam_geo->addDirichletBoundaryCondition({l3}, ALL, 0.0);

Material *mat = new NewtonianFluid(0.001, 1.0);

FluidDomain *dam_problem = new FluidDomain(dam_geo);

dam_problem->applyMaterial({s0}, mat);
dam_problem->generateMesh(T3, FRONT, "teste", "", true, false);
dam_problem->setNumberOfSteps(1000);
dam_problem->setDeltat(0.01);
dam_problem->setMaxNonlinearIterations(6);
dam_problem->setNonlinearTolerance(1.0e-6);
dam_problem->setGravity(0.0, -1.0, 0.0);
dam_problem->setSpectralRadius(0.9);
dam_problem->setMeshLength(1.0/double(h));
dam_problem->solvePFEMProblem();