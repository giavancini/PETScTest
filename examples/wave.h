Geometry *wave_geo = new Geometry(0);

double L = 80.0;
double d = 10.0;
double H = 2.0;
double h = 1.0;

Point *p0 = wave_geo->addPoint({0.0 , 0.0   , 0.0});
Point *p1 = wave_geo->addPoint({L   , 0.0   , 0.0});
Point *p2 = wave_geo->addPoint({L   , d     , 0.0});
Point *p3 = wave_geo->addPoint({-L  , d     , 0.0});
Point *p4 = wave_geo->addPoint({-L  , 0.0   , 0.0});

Line *l0 = wave_geo->addLine({p0, p1});
Line *l1 = wave_geo->addLine({p1, p2});
Line *l2 = wave_geo->addSpline({p2, p3}, [=](double x) 
           { return (d+H*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))*(1.0 / std::cosh(sqrt((3.0*H)/(4.0*d*d*d))*x))); }, 1000);
Line *l3 = wave_geo->addLine({p3, p4});
Line *l4 = wave_geo->addLine({p4, p0});

Surface *s0 = wave_geo->addPlaneSurface({l0, l1, l2, l3, l4});

wave_geo->transfiniteLine({l0}, L / h + 1);
wave_geo->transfiniteLine({l1}, d / h + 1);
wave_geo->transfiniteLine({l2}, 2 * L / h + 1);
wave_geo->transfiniteLine({l3}, d / h + 1);
wave_geo->transfiniteLine({l4}, L / h + 1);
//wave_geo->transfiniteSurface({s0}, "RIGHT", {p1, p2, p3, p4});

wave_geo->addDirichletBoundaryCondition({l0}, Y, 0.0);
wave_geo->addDirichletBoundaryCondition({l4}, Y, 0.0);
wave_geo->addDirichletBoundaryCondition({l1}, X, 0.0);
wave_geo->addDirichletBoundaryCondition({l3}, X, 0.0);

Material *mat = new NewtonianFluid(1.0, 1.0);

FluidDomain *wave_problem = new FluidDomain(wave_geo);

wave_problem->applyMaterial({s0}, mat);
wave_problem->generateMesh(T3, FRONT, "teste", "", true, false);
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