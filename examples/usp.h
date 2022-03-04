Geometry *dam_geo = new Geometry(0);

    double L = 1.0;
    double h = 0.2;
    double pi = 3.14;

    Point *p0 = dam_geo->addPoint({0.0   , 7.0*L});
    Point *p1 = dam_geo->addPoint({0.0   , 2.0*L});
    Point *p2 = dam_geo->addPoint({4.0*L , 2.0*L});
    Point *p3 = dam_geo->addPoint({4.0*L , 7.0*L});
    Point *p4 = dam_geo->addPoint({3.0*L , 7.0*L});
    Point *p5 = dam_geo->addPoint({3.0*L , 2.0*L});
    Point *p6 = dam_geo->addPoint({2.0*L , 2.0*L}, 1.0, false);
    Point *p7 = dam_geo->addPoint({1.0*L , 2.0*L});
    Point *p8 = dam_geo->addPoint({1.0*L , 7.0*L});

    Point *p9 =  dam_geo->addPoint({5.0*L , 0.0});
    Point *p10 = dam_geo->addPoint({7.0*L , 0.0});
    Point *p11 = dam_geo->addPoint({7.0*L , 4.0*L});
    Point *p12 = dam_geo->addPoint({7.0*L , 5.0*L}, 1.0, false);
    Point *p13 = dam_geo->addPoint({7.0*L , 6.0*L});
    Point *p14 = dam_geo->addPoint({9.0*L , 6.0*L});
    Point *p15 = dam_geo->addPoint({9.0*L , 7.0*L});
    Point *p16 = dam_geo->addPoint({7.0*L , 7.0*L});
    Point *p17 = dam_geo->addPoint({7.0*L , 3.0*L});
    Point *p18 = dam_geo->addPoint({7.0*L , 2.0*L}, 1.0, false);
    Point *p19 = dam_geo->addPoint({7.0*L , 1.0*L});
    Point *p20 = dam_geo->addPoint({5.0*L , 1.0*L});

    Point *p21 = dam_geo->addPoint({10.0*L, 0.0});
    Point *p22 = dam_geo->addPoint({11.0*L, 0.0});
    Point *p23 = dam_geo->addPoint({11.0*L, 1.0*L});
    Point *p24 = dam_geo->addPoint({11.0*L, 7.0*L});
    Point *p25 = dam_geo->addPoint({10.0*L, 7.0*L});
    Point *p26 = dam_geo->addPoint({11.0*L, 2.0*L});
    Point *p27 = dam_geo->addPoint({11.0*L, 4.0*L});
    Point *p28 = dam_geo->addPoint({11.0*L, 6.0*L});

    Line *l0 = dam_geo->addLine({p0, p1});
    Line *l1 = dam_geo->addCircle({p1, p6, p2});
    Line *l2 = dam_geo->addLine({p2, p3});
    Line *l3 = dam_geo->addLine({p3, p4});
    Line *l4 = dam_geo->addLine({p4, p5});
    Line *l5 = dam_geo->addCircle({p7, p6, p5});
    Line *l6 = dam_geo->addLine({p7, p8});
    Line *l7 = dam_geo->addLine({p8, p0});

    Line *l8 = dam_geo->addLine({p9, p10});
    Line *l9 = dam_geo->addCircle({p10, p18, p11});
    Line *l10 = dam_geo->addCircle({p13, p12, p11});
    Line *l11 = dam_geo->addLine({p13, p14});
    Line *l12 = dam_geo->addLine({p14, p15});
    Line *l13 = dam_geo->addLine({p15, p16});
    Line *l14 = dam_geo->addCircle({p16, p12, p17});
    Line *l15 = dam_geo->addCircle({p19, p18, p17});
    Line *l16 = dam_geo->addLine({p19, p20});
    Line *l17 = dam_geo->addLine({p20, p9});

    Line *l18 = dam_geo->addLine({p21, p22});
    Line *l19 = dam_geo->addLine({p22, p23});
    Line *l20 = dam_geo->addCircle({p23, p27, p24});
    Line *l21 = dam_geo->addLine({p24, p25});
    Line *l22 = dam_geo->addLine({p25, p21});
    Line *l23 = dam_geo->addCircle({p26, p27, p28});
    Line *l24 = dam_geo->addLine({p28, p27});
    Line *l25 = dam_geo->addLine({p27, p26});

    Surface *s0 = dam_geo->addPlaneSurface({l0, l1, l2, l3, l4, -(*l5), l6, l7});
    Surface *s1 = dam_geo->addPlaneSurface({l8, l9, -(*l10), l11, l12, l13, l14, -(*l15), l16, l17});
    Surface *s2 = dam_geo->addPlaneSurface({l18, l19, l20, l21, l22, -(*l23), -(*l24), -(*l25)});

    dam_geo->transfiniteLine({l0}, 5.0*L / h + 1);
    dam_geo->transfiniteLine({l1}, pi*2.0*L / h + 1);
    dam_geo->transfiniteLine({l2}, 5.0*L / h + 1);
    dam_geo->transfiniteLine({l3}, 1.0*L / h + 1);
    dam_geo->transfiniteLine({l4}, 5.0*L / h + 1);
    dam_geo->transfiniteLine({l5}, pi*1.0*L / h + 1);
    dam_geo->transfiniteLine({l6}, 5.0*L / h + 1);
    dam_geo->transfiniteLine({l7}, 1.0*L / h + 1);

    dam_geo->transfiniteLine({l8}, 2.0*L / h + 1);
    dam_geo->transfiniteLine({l9}, pi*2.0*L / h + 1);
    dam_geo->transfiniteLine({l10}, pi*1.0*L / h + 1);
    dam_geo->transfiniteLine({l11}, 2.0*L / h + 1);
    dam_geo->transfiniteLine({l12}, 1.0*L / h + 1);
    dam_geo->transfiniteLine({l13}, 2.0*L / h + 1);
    dam_geo->transfiniteLine({l14}, pi*2.0*L / h + 1);
    dam_geo->transfiniteLine({l15}, pi*1.0*L / h + 1);
    dam_geo->transfiniteLine({l16}, 2.0*L / h + 1);
    dam_geo->transfiniteLine({l17}, 1.0*L / h + 1);

    dam_geo->transfiniteLine({l18}, 1.0*L / h + 1);
    dam_geo->transfiniteLine({l19}, 1.0*L / h + 1);
    dam_geo->transfiniteLine({l20}, pi*3.0*L / h + 1);
    dam_geo->transfiniteLine({l21}, 1.0*L / h + 1);
    dam_geo->transfiniteLine({l22}, 7.0*L / h + 1);
    dam_geo->transfiniteLine({l23}, pi*2.0*L / h + 1);
    dam_geo->transfiniteLine({l24}, 2.0*L / h + 1);
    dam_geo->transfiniteLine({l25}, 2.0*L / h + 1);

    Material *mat = new NewtonianFluid(0.001, 1000.0);

    FluidDomain *dam_problem = new FluidDomain(dam_geo);

    dam_problem->applyMaterial({s0}, mat);
    dam_problem->applyMaterial({s1}, mat);
    dam_problem->applyMaterial({s2}, mat);
    
    dam_problem->generateMesh(T3, FRONT, "teste", "", true, false);
