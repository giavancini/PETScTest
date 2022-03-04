Geometry *wave_geo = new Geometry(0);

double L = 80.0;
double d = 10.0;
double H = 2.0;
double w = 40.0;
double h = 2.0;

Point *p0 = wave_geo->addPoint({0.0 , 0.0   , 0.0});
Point *p1 = wave_geo->addPoint({L   , 0.0   , 0.0});
Point *p2 = wave_geo->addPoint({L   , d     , 0.0});
Point *p3 = wave_geo->addPoint({-L  , d     , 0.0});
Point *p4 = wave_geo->addPoint({-L  , 0.0   , 0.0});

Point *p5 = wave_geo->addPoint({0.0 , 0.0   , -w });
Point *p6 = wave_geo->addPoint({L   , 0.0   , -w });
Point *p7 = wave_geo->addPoint({L   , d     , -w });
Point *p8 = wave_geo->addPoint({-L  , d     , -w });
Point *p9 = wave_geo->addPoint({-L  , 0.0   , -w });

Line *l0 = wave_geo->addLine({p0, p1});
Line *l1 = wave_geo->addLine({p1, p2});
Line *l2 = wave_geo->addSpline({p2, p3}, [=](double x) 
           { return (d+H*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))); }, 1000);
Line *l3 = wave_geo->addLine({p3, p4});
Line *l4 = wave_geo->addLine({p4, p0});

Line *l5 = wave_geo->addLine({p5, p6});
Line *l6 = wave_geo->addLine({p6, p7});
Line *l7 = wave_geo->addSpline({p7, p8}, [=](double x) 
           { return (d+H*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))); }, 1000);
Line *l8 = wave_geo->addLine({p8, p9});
Line *l9 = wave_geo->addLine({p9, p5});

Line *l10 = wave_geo->addLine({p1, p6});
Line *l11 = wave_geo->addLine({p2, p7});
Line *l12 = wave_geo->addLine({p3, p8});
Line *l13 = wave_geo->addLine({p4, p9});

Surface *s0 = wave_geo->addPlaneSurface({l0, l1, l2, l3, l4});
Surface *s1 = wave_geo->addPlaneSurface({l5, l6, l7, l8, l9});
Surface *s2 = wave_geo->addPlaneSurface({l4, l0, l10, -(*l5), -(*l9), -(*l13)});
Surface *s3 = wave_geo->addSurface({-(*l2), l11, l7, -(*l12)});
Surface *s4 = wave_geo->addPlaneSurface({l10, l6, -(*l11), -(*l1)});
Surface *s5 = wave_geo->addPlaneSurface({l13, -(*l8), -(*l12), l3});

Volume *v0 = wave_geo->addVolume({s0, s1, s2, s3, s4, s5});

wave_geo->transfiniteLine({l0}, L / h + 1);
wave_geo->transfiniteLine({l1}, d / h + 1);
wave_geo->transfiniteLine({l2}, 2 * L / h + 1);
wave_geo->transfiniteLine({l3}, d / h + 1);
wave_geo->transfiniteLine({l4}, L / h + 1);

wave_geo->transfiniteLine({l5}, L / h + 1);
wave_geo->transfiniteLine({l6}, d / h + 1);
wave_geo->transfiniteLine({l7}, 2 * L / h + 1);
wave_geo->transfiniteLine({l8}, d / h + 1);
wave_geo->transfiniteLine({l9}, L / h + 1);

wave_geo->transfiniteLine({l10}, w / h + 1);
wave_geo->transfiniteLine({l11}, w / h + 1);
wave_geo->transfiniteLine({l12}, w / h + 1);
wave_geo->transfiniteLine({l13}, w / h + 1);
//wave_geo->transfiniteSurface({s0}, "RIGHT", {p1, p2, p3, p4});

wave_geo->addDirichletBoundaryCondition({s0}, Z, 0.0);
wave_geo->addDirichletBoundaryCondition({s1}, Z, 0.0);
wave_geo->addDirichletBoundaryCondition({s2}, Y, 0.0);
wave_geo->addDirichletBoundaryCondition({s4}, X, 0.0);
wave_geo->addDirichletBoundaryCondition({s5}, X, 0.0);

Material *mat = new NewtonianFluid(1.0, 1.0);

FluidDomain *wave_problem = new FluidDomain(wave_geo);

wave_problem->applyMaterial({v0}, mat);
wave_problem->generateMesh(TET4, FRONT, "teste", "", true, false);
wave_problem->setNumberOfSteps(5000);
wave_problem->setDeltat(0.01);
wave_problem->setMaxNonlinearIterations(6);
wave_problem->setNonlinearTolerance(1.0e-6);
wave_problem->setGravity(0.0, -9.81, 0.0);
wave_problem->setSpectralRadius(0.8);
wave_problem->setInitialVelocityX([=](double x, double y, double z) 
              {if (x != -L && x != L)
                return (sqrt(9.81*d)*(H/d)*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x)));
               else
                return 0.0;
              });
wave_problem->setInitialVelocityY([=](double x, double y, double z) 
              {if (y != 0.0)
                return (sqrt(3.0*9.81*d)*sqrt((H/d)*(H/d)*(H/d))*(y/d)*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*
                       (std::tanh(sqrt((3.0*H)/(4.0*d*d*d))*x))); 
               else
                return 0.0;
              });
wave_problem->setInitialPressure([=](double x, double y, double z) 
              {return -1.0 * 9.81 * ((d+H*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))) - y); });
wave_problem->solveTransientProblem();