Func temporal_derivative(Func input,Expr temOrd) {
    // This program computes 0th,1st and 2nd order temporal derivative of an image sequence

    // Initialize temporal filters (0rd,1st,2nd) order and spatial filters ()
    uint8_t tn = 23; // temporal window
    uint8_t alpha = 10;
    float tau = 0.25;

    ///////////////////////////////
    // Generate temporal filters //
    ///////////////////////////////

    //@Cat: try to generate a temporal filter [1:1:23] mask with Halide
    //It is a over-complicated way to do such tasks.

    // // Pure definition: do nothing.
    // Func d_ltm;
    // d_ltm(t) = {undef<float>(),undef<float>(),undef<float>()};

    // // Update 0: set the top row of values to 0
    // d_ltm(0) = {0,0,0};

    // // Update 1: generate filter across the range [1:1:23]
    // RDom rt(1,tn-1);
    // Expr d_ltm0,d_ltm1,d_ltm2,SFILT;
    // d_ltm0 = exp(-pow(log(rt/alpha)/tau,2))/float(sqrt(3.141592)*alpha*tau*exp((pow(tau,2)/4)));
    // d_ltm1 = -2*((log(rt/alpha))/(float(pow(tau,2))*rt)) * d_ltm0;
    // d_ltm2 = -2*((log(rt/alpha))/(float(pow(tau,2))*rt)) * d_ltm1 - (2/(pow(tau*rt,2))) * (1-log(rt/alpha)) * d_ltm0;

    // d_ltm(rt) = {d_ltm0,d_ltm1,d_ltm2};

    // // Schedule it for running
    // d_ltm.compute_root();
    // Realization r = d_ltm.realize(tn);

    //@Cat: generate a mask by a normal for-loop

    Image<float> d_ltm0(22), d_ltm1(22), d_ltm2(22); // 23 is the size of temporay window
    d_ltm0(0) = float(0); d_ltm1(0) = float(0); d_ltm2(0) = float(0);
    for (int rt = 1; rt<tn ; rt++) {
        d_ltm0(rt) = exp(-pow(log(rt/alpha)/tau,2))/float(sqrt(M_PI)*alpha*tau*exp((pow(tau,2)/4)));
        d_ltm1(rt) = -2*((log(rt/alpha))/(float(pow(tau,2))*rt)) * d_ltm0(rt);
        d_ltm2(rt) = -2*((log(rt/alpha))/(float(pow(tau,2))*rt)) * d_ltm1(rt) - (2/(pow(tau*rt,2))) * (1-log(rt/alpha)) * d_ltm0(rt);
    }

    RDom rt(0,tn-1);

    //////////////////////////////////////////////////////////////////////
    // Apply derivative filter along temporal domain of input sequences //
    //////////////////////////////////////////////////////////////////////
    Func T("T"); // equivalent to T0,T1,T2 in ColorVideo_speed_angle

    // // Separate T
    // T(x,y,c,t) = {undef<float>(),undef<float>(),undef<float>()};
    // Expr T0("T0"),T1("T1"),T2("T2");
    // // blur2(x, y) = sum(tent(r.x, r.y) * input(x + r.x - 1, y + r.y - 1));
    // T0 = sum(rt,d_ltm0(rt.x)*input(x,y,c,t + rt.x - 1),"sum_T0");
    // T1 = sum(rt,d_ltm1(rt.x)*input(x,y,c,t + rt.x - 1),"sum_T1");
    // T2 = sum(rt,d_ltm2(rt.x)*input(x,y,c,t + rt.x - 1),"sum_T2");

    // T(x,y,c,t) = Tuple(T0,T1,T2);

    // Vector of T
    std::vector<Expr> T_expr(3,undef<float>());

    T_expr[0] = sum(rt,d_ltm0(rt.x)*input(x,y,c,t + rt.x - 1),"sum_T0");
    T_expr[1] = sum(rt,d_ltm1(rt.x)*input(x,y,c,t + rt.x - 1),"sum_T1");
    T_expr[2] = sum(rt,d_ltm2(rt.x)*input(x,y,c,t + rt.x - 1),"sum_T2");

    T(x,y,c,t) = Tuple(T_expr);

    // T.trace_stores();
    // T.reorder(t,x,y,c);
    return T;
}
