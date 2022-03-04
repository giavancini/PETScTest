//Fluid problem

double h = 0.012/2.0;
double L = 0.146;
double H = 0.08;
double D = 0.012;

Geometry *fluid_geo = new Geometry(0);

Point *p0  = fluid_geo->addPoint({0.0     , 0.0   , -L});
Point *p1  = fluid_geo->addPoint({L       , 0.0   , -L});
Point *p2  = fluid_geo->addPoint({L       , 2.0*L , -L});
Point *p3  = fluid_geo->addPoint({0.0     , 2.0*L , -L});
Point *p4  = fluid_geo->addPoint({0.0     , 4.0*L , -L});
Point *p5  = fluid_geo->addPoint({2.0*L   , 0.0   , -L});
Point *p6  = fluid_geo->addPoint({2.0*L+D , 0.0   , -L});
Point *p7  = fluid_geo->addPoint({4.0*L+D , 0.0   , -L});
Point *p8  = fluid_geo->addPoint({4.0*L+D , 4.0*L , -L});

Point *p9  = fluid_geo->addPoint({0.0     , 0.0   , 0.0});
Point *p10 = fluid_geo->addPoint({L       , 0.0   , 0.0});
Point *p11 = fluid_geo->addPoint({L       , 2.0*L , 0.0});
Point *p12 = fluid_geo->addPoint({0.0     , 2.0*L , 0.0});
Point *p13 = fluid_geo->addPoint({0.0     , 4.0*L , 0.0});
Point *p14 = fluid_geo->addPoint({2.0*L   , 0.0   , 0.0});
Point *p15 = fluid_geo->addPoint({2.0*L+D , 0.0   , 0.0});
Point *p16 = fluid_geo->addPoint({4.0*L+D , 0.0   , 0.0});
Point *p17 = fluid_geo->addPoint({4.0*L+D , 4.0*L , 0.0});

Point *p18 = fluid_geo->addPoint({2.0*L   , 0.0   , -0.99*L});
Point *p19 = fluid_geo->addPoint({2.0*L   , H     , -0.99*L});
Point *p20 = fluid_geo->addPoint({2.0*L+D , H     , -0.99*L});
Point *p21 = fluid_geo->addPoint({2.0*L+D , 0.0   , -0.99*L});
Point *p22 = fluid_geo->addPoint({2.0*L   , 0.0   , -0.01*L});
Point *p23 = fluid_geo->addPoint({2.0*L   , H     , -0.01*L});
Point *p24 = fluid_geo->addPoint({2.0*L+D , H     , -0.01*L});
Point *p25 = fluid_geo->addPoint({2.0*L+D , 0.0   , -0.01*L});

Line *l0  = fluid_geo->addLine({p0, p1});
Line *l1  = fluid_geo->addLine({p1, p2});
Line *l2  = fluid_geo->addLine({p2, p3});
Line *l3  = fluid_geo->addLine({p3, p0});
Line *l4  = fluid_geo->addLine({p3, p4});
Line *l5  = fluid_geo->addLine({p1, p5});
Line *l6  = fluid_geo->addLine({p5, p6});
Line *l7  = fluid_geo->addLine({p6, p7});
Line *l8  = fluid_geo->addLine({p7 , p8});
Line *l9  = fluid_geo->addLine({p8 , p4});
Line *l10 = fluid_geo->addLine({p9 , p10});
Line *l11 = fluid_geo->addLine({p10, p11});
Line *l12 = fluid_geo->addLine({p11, p12});
Line *l13 = fluid_geo->addLine({p12, p9});
Line *l14 = fluid_geo->addLine({p12, p13});
Line *l15 = fluid_geo->addLine({p10, p14});
Line *l16 = fluid_geo->addLine({p14, p15});
Line *l17 = fluid_geo->addLine({p15, p16});
Line *l18 = fluid_geo->addLine({p16, p17});
Line *l19 = fluid_geo->addLine({p17, p13});
Line *l20 = fluid_geo->addLine({p18, p19});
Line *l21 = fluid_geo->addLine({p19, p20});
Line *l22 = fluid_geo->addLine({p20, p21});
Line *l23 = fluid_geo->addLine({p21, p18});
Line *l24 = fluid_geo->addLine({p22, p23});
Line *l25 = fluid_geo->addLine({p23, p24});
Line *l26 = fluid_geo->addLine({p24, p25});
Line *l27 = fluid_geo->addLine({p25, p22 });
Line *l28 = fluid_geo->addLine({p0 , p9});
Line *l29 = fluid_geo->addLine({p3 , p12});
Line *l30 = fluid_geo->addLine({p4 , p13});
Line *l31 = fluid_geo->addLine({p1 , p10});
Line *l32 = fluid_geo->addLine({p2 , p11});
Line *l33 = fluid_geo->addLine({p5 , p18});
Line *l34 = fluid_geo->addLine({p18, p22});
Line *l35 = fluid_geo->addLine({p22, p14});
Line *l36 = fluid_geo->addLine({p19, p23});
Line *l37 = fluid_geo->addLine({p6 , p21});
Line *l38 = fluid_geo->addLine({p21, p25});
Line *l39 = fluid_geo->addLine({p25, p15});
Line *l40 = fluid_geo->addLine({p20, p24});
Line *l41 = fluid_geo->addLine({p7 , p16});
Line *l42 = fluid_geo->addLine({p8 , p17});

Surface *s0 = fluid_geo->addPlaneSurface({l0, l1, l2, l3});
Surface *s1 = fluid_geo->addPlaneSurface({l5, l6, l7, l8, l9, -(*l4), -(*l2), -(*l1)});
Surface *s2 = fluid_geo->addPlaneSurface({l10, l11, l12, l13});
Surface *s3 = fluid_geo->addPlaneSurface({l15, l16, l17, l18, l19, -(*l14), -(*l12), -(*l11)});
Surface *s4 = fluid_geo->addPlaneSurface({l23, l20, l21, l22});
Surface *s5 = fluid_geo->addPlaneSurface({l27, l24, l25, l26});
Surface *s6 = fluid_geo->addPlaneSurface({l28, -(*l13), -(*l29), l3});
Surface *s7 = fluid_geo->addPlaneSurface({l29, l14, -(*l30), -(*l4)});
Surface *s8 = fluid_geo->addPlaneSurface({l31, l11, -(*l32), -(*l1)});
Surface *s9 = fluid_geo->addPlaneSurface({l34, l24, -(*l36), -(*l20)});
Surface *s10 = fluid_geo->addPlaneSurface({l38, -(*l26), -(*l40), l22});
Surface *s11 = fluid_geo->addPlaneSurface({l41, l18, -(*l42), -(*l8)});
Surface *s12 = fluid_geo->addPlaneSurface({l10, -(*l31), -(*l0), l28});
Surface *s13 = fluid_geo->addPlaneSurface({l15, -(*l35), -(*l34), -(*l33), -(*l5), l31});
Surface *s14 = fluid_geo->addPlaneSurface({l33, -(*l23), -(*l37), -(*l6)});
Surface *s15 = fluid_geo->addPlaneSurface({l27, l35, l16, -(*l39)});
Surface *s16 = fluid_geo->addPlaneSurface({l7, l41, -(*l17), -(*l39), -(*l38), -(*l37)});
Surface *s17 = fluid_geo->addPlaneSurface({l21, l40, -(*l25), -(*l36)});
Surface *s18 = fluid_geo->addPlaneSurface({l2, l29, -(*l12), -(*l32)});
Surface *s19 = fluid_geo->addPlaneSurface({l9, l30, -(*l19), -(*l42)});

Volume *v0 = fluid_geo->addVolume({s2, s12, s0, s6, s18, s8});

fluid_geo->transfiniteLine({l0}, L / h + 1);
fluid_geo->transfiniteLine({l1}, 2*L / h + 1);
fluid_geo->transfiniteLine({l2}, L / h + 1);
fluid_geo->transfiniteLine({l3}, 2*L / h + 1);
fluid_geo->transfiniteLine({l4}, 2*L / h + 1);
fluid_geo->transfiniteLine({l5}, L / h + 1);
fluid_geo->transfiniteLine({l6}, D / h + 1);
fluid_geo->transfiniteLine({l7}, 2*L / h + 1);
fluid_geo->transfiniteLine({l8}, 4*L / h + 1);
fluid_geo->transfiniteLine({l9}, (4*L+D) / h + 1);

fluid_geo->transfiniteLine({l10}, L / h + 1);
fluid_geo->transfiniteLine({l11}, 2*L / h + 1);
fluid_geo->transfiniteLine({l12}, L / h + 1);
fluid_geo->transfiniteLine({l13}, 2*L / h + 1);
fluid_geo->transfiniteLine({l14}, 2*L / h + 1);
fluid_geo->transfiniteLine({l15}, L / h + 1);
fluid_geo->transfiniteLine({l16}, D / h + 1);
fluid_geo->transfiniteLine({l17}, 2*L / h + 1);
fluid_geo->transfiniteLine({l18}, 4*L / h + 1);
fluid_geo->transfiniteLine({l19}, (4*L+D) / h + 1);

fluid_geo->transfiniteLine({l20}, H / h + 1);
fluid_geo->transfiniteLine({l21}, D / h + 1);
fluid_geo->transfiniteLine({l22}, H / h + 1);
fluid_geo->transfiniteLine({l23}, D / h + 1);
fluid_geo->transfiniteLine({l24}, H / h + 1);
fluid_geo->transfiniteLine({l25}, D / h + 1);
fluid_geo->transfiniteLine({l26}, H / h + 1);
fluid_geo->transfiniteLine({l27}, D / h + 1);

fluid_geo->transfiniteLine({l28}, L / h + 1);
fluid_geo->transfiniteLine({l29}, L / h + 1);
fluid_geo->transfiniteLine({l30}, L / h + 1);
fluid_geo->transfiniteLine({l31}, L / h + 1);
fluid_geo->transfiniteLine({l32}, L / h + 1);

fluid_geo->transfiniteLine({l33}, 2);
fluid_geo->transfiniteLine({l34}, L / h + 1);
fluid_geo->transfiniteLine({l35}, 2);
fluid_geo->transfiniteLine({l36}, L / h + 1);
fluid_geo->transfiniteLine({l37}, 2);
fluid_geo->transfiniteLine({l38}, L / h + 1);
fluid_geo->transfiniteLine({l39}, 2);
fluid_geo->transfiniteLine({l40}, L / h + 1);
fluid_geo->transfiniteLine({l41}, L / h + 1);
fluid_geo->transfiniteLine({l42}, L / h + 1);

fluid_geo->transfiniteSurface({s9}, "Alternate");
fluid_geo->transfiniteSurface({s10}, "Alternate");
fluid_geo->transfiniteSurface({s4}, "Alternate");
fluid_geo->transfiniteSurface({s5}, "Alternate");
fluid_geo->transfiniteSurface({s17}, "Alternate");

fluid_geo->addDirichletBoundaryCondition({s0}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s1}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s2}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s3}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s6}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s7}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s11}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s12}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s13}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s14}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s15}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s16}, ALL, 0.0);
fluid_geo->addDirichletBoundaryCondition({s19}, ALL, 0.0);
fluid_geo->addInterfaceBoundaryCondition({s4});
fluid_geo->addInterfaceBoundaryCondition({s5});
fluid_geo->addInterfaceBoundaryCondition({s9});
fluid_geo->addInterfaceBoundaryCondition({s10});
fluid_geo->addInterfaceBoundaryCondition({s17});
Material *fluid_mat = new NewtonianFluid(0.001, 1000.0);
FluidDomain *fluid_problem = new FluidDomain(fluid_geo);
fluid_problem->applyMaterial({v0}, fluid_mat);
fluid_problem->generateMesh(TET4, FRONT, "fluid", "", true, true);
fluid_problem->setGravity(0.0, -9.81, 0.0);
fluid_problem->setMeshLength(h);
fluid_problem->setAlpha(1.35);

//Solid problem
Geometry *solid_geo = new Geometry(1);

Point *p26 = solid_geo->addPoint({2.0*L   , 0.0   , -0.99*L});
Point *p27 = solid_geo->addPoint({2.0*L   , H     , -0.99*L});
Point *p28 = solid_geo->addPoint({2.0*L+D , H     , -0.99*L});
Point *p29 = solid_geo->addPoint({2.0*L+D , 0.0   , -0.99*L});
Point *p30 = solid_geo->addPoint({2.0*L   , 0.0   , -0.01*L});
Point *p31 = solid_geo->addPoint({2.0*L   , H     , -0.01*L});
Point *p32 = solid_geo->addPoint({2.0*L+D , H     , -0.01*L});
Point *p33 = solid_geo->addPoint({2.0*L+D , 0.0   , -0.01*L});

Line *l43 =  solid_geo->addLine({p26, p27});
Line *l44 =  solid_geo->addLine({p27, p28});
Line *l45 =  solid_geo->addLine({p28, p29});
Line *l46 =  solid_geo->addLine({p29, p26});
Line *l47 =  solid_geo->addLine({p30, p31});
Line *l48 =  solid_geo->addLine({p31, p32});
Line *l49 =  solid_geo->addLine({p32, p33});
Line *l50 =  solid_geo->addLine({p33, p30});
Line *l51 =  solid_geo->addLine({p26, p30});
Line *l52 =  solid_geo->addLine({p27, p31});
Line *l53 =  solid_geo->addLine({p28, p32});
Line *l54 =  solid_geo->addLine({p29, p33});

Surface *s20 = solid_geo->addPlaneSurface({l43, l44, l45, l46});
Surface *s21 = solid_geo->addPlaneSurface({-(*l47), -(*l50), -(*l49), -(*l48)});
Surface *s22 = solid_geo->addPlaneSurface({l51, l47, -(*l52), -(*l43)});
Surface *s23 = solid_geo->addPlaneSurface({-(*l54), -(*l45), l53, l49});
Surface *s24 = solid_geo->addPlaneSurface({l50, -(*l51), -(*l46), l54});
Surface *s25 = solid_geo->addPlaneSurface({l48, -(*l53), -(*l44), l52});

Volume *v1 = solid_geo->addVolume({s20, s21, s22, s23, s24, s25});

solid_geo->transfiniteLine({l43}, H / h + 1);
solid_geo->transfiniteLine({l44}, D / h + 1);
solid_geo->transfiniteLine({l45}, H / h + 1);
solid_geo->transfiniteLine({l46}, D / h + 1);
solid_geo->transfiniteLine({l47}, H / h + 1);
solid_geo->transfiniteLine({l48}, D / h + 1);
solid_geo->transfiniteLine({l49}, H / h + 1);
solid_geo->transfiniteLine({l50}, D / h + 1);
solid_geo->transfiniteLine({l51}, L / h + 1);
solid_geo->transfiniteLine({l52}, L / h + 1);
solid_geo->transfiniteLine({l53}, L / h + 1);
solid_geo->transfiniteLine({l54}, L / h + 1);

solid_geo->transfiniteSurface({s20}, "Left");
solid_geo->transfiniteSurface({s21}, "Left");
solid_geo->transfiniteSurface({s22}, "Left");
solid_geo->transfiniteSurface({s23}, "Left");
solid_geo->transfiniteSurface({s24}, "Left");
solid_geo->transfiniteSurface({s25}, "Left");

solid_geo->addDirichletBoundaryCondition({s24}, ALL, 0.0);
solid_geo->addInterfaceBoundaryCondition({s20});
solid_geo->addInterfaceBoundaryCondition({s21});
solid_geo->addInterfaceBoundaryCondition({s22});
solid_geo->addInterfaceBoundaryCondition({s23});
solid_geo->addInterfaceBoundaryCondition({s25});

Material *solid_mat = new ElasticSolid(1.0e6, 0.0, 2500.0);
SolidDomain *solid_problem = new SolidDomain(solid_geo);
solid_problem->applyMaterial({v1}, solid_mat);
solid_problem->generateMesh(TET4, FRONT, "solid", "", true, true);
solid_problem->setGravity(0.0, 0.0, 0.0);

//Coupled problem
CoupledDomain *coupled_problem = new CoupledDomain(fluid_problem, solid_problem);
coupled_problem->setNumberOfSteps(5000);
coupled_problem->setDeltat(0.0002);
coupled_problem->setMaxNonlinearIterations(6);
coupled_problem->setNonlinearTolerance(1.0e-4);
coupled_problem->setSpectralRadius(0.8);
coupled_problem->setExportFrequency(10);
coupled_problem->solveMonolithicCoupledPFEMProblem();