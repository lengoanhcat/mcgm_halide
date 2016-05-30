Func spatial_derivative(Func T, Expr spaOrd, Expr temOrd) {
// This program computers spatial gaussian derivative up to nfilt-oder

    uint8_t nfilt = 5;
    uint8_t numSTB = 63; // number of spatial-temporal basis
    uint8_t numTB = 3;
    uint8_t sn = 23; // spatial window
    float sigma = 1.5;

    //////////////////////////////
    // Generate spatial filters //
    //////////////////////////////

    Image<float> SFILT(sn,nfilt+1);
    for (int rs = 0; rs < sn; rs++) {
        float si = rs - (sn-1)/2;
        float gauss = exp(-(pow(si,2)/(4*sigma)))/sqrt(4*sigma*3.1415926);
        for (int ro = 0; ro <= nfilt; ro++) {
            float H = 0;
            for (int io = 0; io < ro; io++) {
                H = H + pow((-1),io) * (pow((2*si),(ro-2*io))/pow((4*sigma),(ro-io)))/(factorial(io)*factorial(ro-2*io));
            }
            H = factorial(ro)*H;
            SFILT(rs,ro) = H*gauss;
        }
    }

    // create a new range
   RDom rs(0,sn-1);

    /////////////////////////////////////////////////
    // Apply spatial derivative filter on T0,T1,T2 //
    /////////////////////////////////////////////////

    // // Version 1: does not generate anything
    Func tmp_col("tmp_col");
    std::vector<Expr> tmp_col_expr((nfilt+1)*numTB,undef<float>()); // define a vector of expression
    Func basis("basis");
    std::vector<Expr> basis_expr(numSTB,undef<float>());
    // basis(x,y,c,t) = Tuple(basis_expr);
    int iB = 0;

    // for (int iSf = 0; iSf < numSTB; iSf++)
    //     basis_expr[iSf] = Expr(1.0f);

    // FIR filter on horizontal axis
    for (int iTf = 0; iTf < numTB; iTf++ ) { // iTf: index of temporal filter, we only take the first two orders
        for (int iSf=0; iSf<=nfilt; iSf++){ // iSf: index of spatial filter
            tmp_col_expr[iTf*(nfilt+1)+iSf] = sum(rs,T(x + rs.x - 1,y,c,t)[iTf]*SFILT(rs.x,iSf),"sum_col");
        }
    }
    tmp_col(x,y,c,t) = Tuple(tmp_col_expr);

    // FIR filter on vertical axis
    for (int iTf = 0; iTf < numTB; iTf++)
        for (int iSf = 0; iSf<=nfilt; iSf++)
            for (int iSf1 = 0; iSf1<=iSf; iSf1++)
                for (int iSf2 = 0; iSf2<=iSf; iSf2++)
                    if ((iSf1+iSf2) == iSf) {
                        basis_expr[iB] = sum(rs,tmp_col(x,y + rs.x - 1,c,t)[iTf*(nfilt+1)+iSf1]*SFILT(rs.x,iSf2),"sum_basis");
                        iB++;
                    }

    basis(x,y,c,t) = Tuple(basis_expr);
    return basis;

    // Version 2:
    // Func tmp_col("tmp_col");
    // std::vector<Expr> tmp_col_expr(nfilt+1,undef<float>()); // define a vector of expression
    // tmp_col(x,y,c,t) = Tuple(tmp_col_expr);
    // Func basis("basis");
    // std::vector<Expr> basis_expr(numSTB,undef<float>());

    // int iB = 0;
    // for (int iTf = 0; iTf < 3; iTf++ ) { // iTf: index of temporal filter, we only take the first two orders
    //     for (int iSf=0; iSf<=nfilt; iSf++){ // iSf: index of spatial filter
    //         tmp_col_expr[iSf] = sum(rs,T(x + rs.x - 1,y,c,t)[iTf]*SFILT(rs.x,iSf),"sum_col");
    //     }
    //     for (int iSf = 0; iSf<=nfilt; iSf++)
    //         for (int iSf1 = 0; iSf1<=iSf; iSf1++)
    //             for (int iSf2 = 0; iSf2<=iSf; iSf2++)
    //                 if ((iSf1+iSf2) == nfilt) {
    //                     basis_expr[iB] = sum(rs,tmp_col_expr(x,y + rs.x - 1,c,t)[iSf1]*SFILT(rs.x,iSf2),"sum_basis");
    //                 }
    // }

    // basis(x,y,c,t) = Tuple(basis_expr);
    // return basis;
}
