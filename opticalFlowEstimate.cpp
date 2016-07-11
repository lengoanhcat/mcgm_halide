int Mgetfilterindex(uint8_t x_order,uint8_t y_order,uint8_t t_order,uint8_t numSTB, uint8_t numSB) {
// Figure out filter index with respect to filter order (x_order,y_order,t_order) in which (x,y) are spatial order and (t) is temporal orders.
    // numSTB: number of spatio-temporal basis

    int index = (x_order+y_order)*(x_order+y_order+1)/2 + y_order;
    // int numSB = numSTB / numTB; // numTB = 3 numSB = 21 (5th order in both x and y) -> numSTB = numTB*numSB = 63
    if (index >= numSB) index = -64;
    index = index + t_order * numSB;
    assert(index<numSTB);
    return index;
}

Func ColorMgetfilter(Func stBasis, float angle, uint8_t iXo, uint8_t iYo, uint8_t iTo, uint8_t iCo ) {
    // Compute a rotated basis at (iXo,iYo,iTo,iCo) order with angle value

    // temporary setting
    uint8_t numSTB = 63;
    uint8_t numSB = 21;

    angle = -1*angle - M_PI/2;
    float * weights;

    Func work; // work: rotated basis at a particular spatio-temporal order
    work(x,y,t) = cast<float>(0.0f);

    weights = (float *) calloc(iXo+iYo+1,sizeof(float));

    // compute weights for possible orders
    for (int i = 0; i <= iXo; i++)
        for (int j = 0; j <= iYo; j++)
            weights[iXo+iYo-i-j] += float(combination(iXo,i))*float(combination(iYo,j))*pow((-1.0f),float(i))*pow(cos(angle),float(iXo-i+j))*pow(sin(angle),float(iYo+i-j));

    // get filtered expression at paricular order and angle value
    // Func basis("basis");
    for (int k=0; k<=(iXo+iYo); k++) {
        int index = Mgetfilterindex(iXo+iYo-k,k,iTo,numSTB,numSB);
        // basis = spatial_temporal_derivative(T,iXo+iYo-k,k,iTo,iCo);
        if ((index > 0) && (weights[iXo+iYo-k] != 0))
            work(x,y,t) += weights[iXo+iYo-k]*stBasis(x,y,iCo,t)[index];
    }

    work.compute_root();

    free(weights);
    return work;
}

Func ColorMgather(Func stBasis, float angle, uint8_t * orders, Expr filterthreshold, Expr divisionthreshold, Expr divisionthreshold2) {
    uint8_t x_order = orders[0];
    uint8_t y_order = orders[1];
    uint8_t t_order = orders[2];
    uint8_t c_order = orders[3];

    Func X("X"),Y("Y"),T("T"),Xrg("Xrg"),Yrg("Yrg"),Trg("Trg");
    uint8_t max_order = x_order;
    // std::vector<Expr>Xk_expr (max_order,cast<float>(0.0f));
    // std::vector<Expr>Yk_expr (max_order,cast<float>(0.0f));
    // std::vector<Expr>Tk_expr (max_order,cast<float>(0.0f));

    uint8_t Xk_uI[max_order];
    uint8_t Yk_uI[max_order];
    uint8_t Tk_uI[max_order];
    Func Xk[max_order]; Func Yk[max_order]; Func Tk[max_order];

    // Expr Xk[max_order],Yk[max_order],Tk[max_order];

    for (int iO=0; iO < x_order; iO++) {
        Xk[iO](x,y,t) = Expr(0.0f);
        Yk[iO](x,y,t) = Expr(0.0f);
        Tk[iO](x,y,t) = Expr(0.0f);
        Xk_uI[iO] = 0;
        Yk_uI[iO] = 0;
        Tk_uI[iO] = 0;
    }

    int k = 0;
    for (int iXo = 0; iXo < x_order; iXo++) // x_order
        for (int iYo = 0; iYo < y_order; iYo++) // y_oder
            for (int iTo = 0; iTo < t_order; iTo++) // t_order
                for (int iCo = 0; iCo < c_order; iCo ++ ) // c_order: index of color channel
                {
                    if ((iYo+iTo+iCo == 0 || iYo+iTo+iCo == 1)
                        && ((iXo+iYo+iTo+iCo+1) < (x_order + 1))) {
                        X = ColorMgetfilter(stBasis, angle, iXo+1, iYo, iTo, iCo);
                        Y = ColorMgetfilter(stBasis, angle, iXo, iYo+1, iTo, iCo);
                        T = ColorMgetfilter(stBasis, angle, iXo, iYo, iTo+1, iCo);

                        Xrg = ColorMgetfilter(stBasis, angle, iXo+1, iYo, iTo, iCo+1);
                        Yrg = ColorMgetfilter(stBasis, angle, iXo, iYo+1, iTo, iCo+1);
                        Trg = ColorMgetfilter(stBasis, angle, iXo, iYo, iTo+1, iCo+1);

                        k = iXo + iYo + iTo + iCo;
                        Xk[k](x,y,t) += X(x,y,t) + Xrg(x,y,t);
                        Yk[k](x,y,t) += Y(x,y,t) + Yrg(x,y,t);
                        Tk[k](x,y,t) += T(x,y,t) + Trg(x,y,t);

                        Xk[k].update(Xk_uI[k]); Xk_uI[k]++;
                        Yk[k].update(Yk_uI[k]); Yk_uI[k]++;
                        Tk[k].update(Tk_uI[k]); Tk_uI[k]++;
                    }
                }

    // Scheduling
    for (int iO = 0; iO <= k; iO++) {
        Xk[iO].compute_root();
        Yk[iO].compute_root();
        Tk[iO].compute_root();
    }

    std::vector<Expr> st_expr(6,cast<float>(0.0f));
    for (int iK=0; iK <= k; iK++) {
        st_expr[0] += Xk[iK](x,y,t)*Tk[iK](x,y,t);
        st_expr[1] += Tk[iK](x,y,t)*Tk[iK](x,y,t);
        st_expr[2] += Xk[iK](x,y,t)*Xk[iK](x,y,t);
        st_expr[3] += Yk[iK](x,y,t)*Tk[iK](x,y,t);
        st_expr[4] += Yk[iK](x,y,t)*Yk[iK](x,y,t);
        st_expr[5] += Xk[iK](x,y,t)*Yk[iK](x,y,t);
    }

    Func st("st"); st(x,y,t) = Tuple(st_expr);
    st.compute_root();

    Expr x_clamped = clamp(x,0,width-1);
    Expr y_clamped = clamp(y,0,height-1);
    Func st_clamped("st_clamped"); st_clamped(x,y,t) = st(x_clamped,y_clamped,t);

    // float win = 7.0;
    // Image<float> meanfilter(7,7,"meanfilter_data");
    // meanfilter(x,y) = Expr(1.0f/(win*win));

    // RDom rMF(meanfilter);
    uint8_t win = 7;
    RDom rMF(0,win,0,win);

    Func st_filtered[6];
    for (uint8_t iPc=0; iPc<6; iPc++) {
    // iPc: index of product component
        // Apply average filter
        st_filtered[iPc](x,y,t) = sum(rMF,st_clamped(x + rMF.x,y + rMF.y,t)[iPc]/Expr(float(win*win)),"mean_filter");
        st_filtered[iPc].compute_root();
    }
    // Tuple st_tuple = Tuple(st_expr);

    // 4 debug
    // Func tmpOut("tmpOut"); tmpOut(x,y,t) = Tuple(st_filtered[0](x,y,t),st_filtered[1](x,y,t),st_filtered[2](x,y,t),st_filtered[3](x,y,t),st_filtered[4](x,y,t),st_filtered[5](x,y,t));
    // return tmpOut;

    Tuple pbx = Tuple(st_filtered[2](x,y,t),st_filtered[5](x,y,t),st_filtered[0](x,y,t));
    Tuple pby = Tuple(st_filtered[5](x,y,t),st_filtered[4](x,y,t),st_filtered[3](x,y,t));
    Tuple pbt = Tuple(st_filtered[0](x,y,t),st_filtered[3](x,y,t),st_filtered[1](x,y,t));

    Func pbxy("pbxy"); pbxy = cross(pby,pbx); pbxy.compute_root();
    Func pbxt("pbxt"); pbxt = cross(pbx,pbt); pbxt.compute_root();
    Func pbyt("pbyt"); pbyt = cross(pby,pbt); pbyt.compute_root();

    Func pbxyd("pbxyd"); pbxyd = dot(pby,pbx); pbxyd.compute_root();
    Func pbxtd("pbxtd"); pbxtd = dot(pbx,pbt); pbxtd.compute_root();
    Func pbytd("pbytd"); pbytd = dot(pby,pbt); pbytd.compute_root();

    // 4 debug
    // Func tmpOut("tmpOut"); tmpOut(x,y,t) = Tuple(pbxy(x,y,t)[0],pbxt(x,y,t)[0],pbyt(x,y,t)[0],pbxyd(x,y,t),pbxtd(x,y,t),pbytd(x,y,t));
    // return tmpOut;

    Func yt_xy("yt_xy"); yt_xy = dot(pbyt(x,y,t),pbxy(x,y,t)); yt_xy.compute_root();
    Func xt_yt("xt_yt"); xt_yt = dot(pbxt(x,y,t),pbyt(x,y,t)); xt_yt.compute_root();
    Func xt_xy("xt_xy"); xt_xy = dot(pbxt(x,y,t),pbxy(x,y,t)); xt_xy.compute_root();
    Func yt_yt("yt_yt"); yt_yt = dot(pbyt(x,y,t),pbyt(x,y,t)); yt_yt.compute_root();
    Func xt_xt("xt_xt"); xt_xt = dot(pbxt(x,y,t),pbxt(x,y,t)); xt_xt.compute_root();
    Func xy_xy("xy_xy"); xy_xy = dot(pbxy(x,y,t),pbxy(x,y,t)); xy_xy.compute_root();

    Tuple Tk_tuple = Tuple(Tk[0](x,y,t),Tk[1](x,y,t),Tk[2](x,y,t),
                           Tk[3](x,y,t),Tk[4](x,y,t));
    Func Tkd("Tkd"); Tkd = dot(Tk_tuple,Tk_tuple); Tkd.compute_root();

    // Expr Dimen = pbxyd/xy_xy;
    Expr kill(1.0f);

    Func Oxy; Oxy(x,y,t) = Mdefdiv(st_filtered[5](x,y,t) - Mdefdivang(yt_xy(x,y,t),yt_yt(x,y,t),pbxyd(x,y,t),divisionthreshold2)*st_filtered[3](x,y,t)*kill,st_filtered[4](x,y,t),divisionthreshold);
    Oxy.compute_root();

    Func Oyx; Oyx(x,y,t) = Mdefdiv(st_filtered[5](x,y,t) + Mdefdivang(xt_xy(x,y,t),xt_xt(x,y,t),pbxyd(x,y,t),divisionthreshold2)*st_filtered[0](x,y,t)*kill,st_filtered[2](x,y,t),divisionthreshold);
    Oyx.compute_root();

    Func C0; C0(x,y,t) = st_filtered[3](x,y,t) * Mdefdivang(Expr(-1.0f)*xt_yt(x,y,t),yt_yt(x,y,t),pbxyd(x,y,t),divisionthreshold2)*kill;
    C0.compute_root();

    Func M0; M0(x,y,t) = Mdefdiv(st_filtered[0](x,y,t) + C0(x,y,t), st_filtered[1](x,y,t)*pow(Mdefdivang(xt_yt(x,y,t),yt_yt(x,y,t),pbxyd(x,y,t),divisionthreshold2),Expr(2.0f)),divisionthreshold);
    M0.compute_root();

    Func C1; C1(x,y,t) = st_filtered[5](x,y,t) * Mdefdivang(Expr(-1.0f)*xt_xy(x,y,t),xy_xy(x,y,t),pbxyd(x,y,t),divisionthreshold2)*kill;
    C1.compute_root();

    Func P1; P1(x,y,t) = pow(Mdefdivang(xt_yt(x,y,t),xt_xt(x,y,t),pbxyd(x,y,t),divisionthreshold2),Expr(2.0f))*kill + 1.0f;
    P1.compute_root();

    // 4 debug
    // Func tmpOut("tmpOut"); tmpOut(x,y,t) = Tuple(Oxy(x,y,t),Oyx(x,y,t),C0(x,y,t),M0(x,y,t),C1(x,y,t),P1(x,y,t));
    // return tmpOut;


    Func Q1; Q1(x,y,t) = st_filtered[2](x,y,t) * (pow(Oyx(x,y,t),Expr(2.0f))+Expr(1.0f));
    Q1.compute_root();

    Func M1; M1(x,y,t) = Mdefdiv(((st_filtered[0](x,y,t)-C1(x,y,t))*P1(x,y,t)),Q1(x,y,t),divisionthreshold);
    M1.compute_root();

    Func C2; C2(x,y,t) = st_filtered[0](x,y,t) * Mdefdivang(Expr(-1.0f)*xt_yt(x,y,t),xt_xt(x,y,t),pbxyd(x,y,t),divisionthreshold2)*kill;
    C2.compute_root();

    Func M2; M2(x,y,t) = Mdefdiv(st_filtered[3](x,y,t)+C2(x,y,t),st_filtered[1](x,y,t)*(pow(Mdefdivang(xt_yt(x,y,t),xt_xt(x,y,t),pbxyd(x,y,t),divisionthreshold2),Expr(2.0f))*kill+Expr(1.0f)),divisionthreshold);
    M2.compute_root();

    Func C3; C3(x,y,t) = st_filtered[5](x,y,t) * Mdefdivang(yt_xy(x,y,t),xy_xy(x,y,t),pbxyd(x,y,t),divisionthreshold2)*kill;
    C3.compute_root();

    Func P3; P3(x,y,t) = pow(Mdefdivang(xt_yt(x,y,t),yt_yt(x,y,t),pbxyd(x,y,t),divisionthreshold2),Expr(2.0f))*kill + Expr(1.0f);
    P3.compute_root();

    Func Q3; Q3(x,y,t) = st_filtered[4](x,y,t) * (pow(Oxy(x,y,t),Expr(2.0f))+Expr(1.0f));
    Q3.compute_root();

    Func M3; M3(x,y,t) = Mdefdiv(((st_filtered[3](x,y,t)-C3(x,y,t))*P3(x,y,t)),Q3(x,y,t),divisionthreshold);
    M3.compute_root();

    Func basisAtAngle;
    basisAtAngle(x,y,t) = Tuple(M0(x,y,t),M1(x,y,t),M2(x,y,t),M3(x,y,t),Tkd(x,y,t));

    return basisAtAngle;

// Func hsv2rgb(Func colorImage) { // Took this function
//     Var x, y, c, t;
//     Func output;
//     output(x,y,c,t) = cast <float> (0.0f);
//     Expr fR, fG, fB; // R,G & B values
//     Expr fH = (colorImage(x,y,0,t)); //H value [0-360)
//     Expr fS = (colorImage(x,y,1,t)); //S value
//     Expr fV = (colorImage(x,y,2,t)); //V value

// //Conversion (I took the one on Wikipedia)
//     // https://fr.wikipedia.org/wiki/Teinte_Saturation_Valeur#Conversion_de_TSV_vers_RVB
//     Expr fHi = floor(fH / Expr(60.0f));
//     Expr fF = fH / 60.0f - fHi;
//     Expr fL = fV * (1 - fS);
//     Expr fM = fV * (1 - fF * fS) ;
//     Expr fN = fV * (1 - (1 - fF) * fS);

//     fR = select((0 == fHi),fV,
//                 (1 == fHi),fM,
//                 (2 == fHi),fL,
//                 (3 == fHi),fL,
//                 (4 == fHi),fN,
//                 (5 == fHi),fV,
//                 0.0f);

//     fG = select((0 == fHi),fN,
//                 (1 == fHi),fV,
//                 (2 == fHi),fV,
//                 (3 == fHi),fM,
//                 (4 == fHi),fL,
//                 (5 == fHi),fL,
//                 0.0f);

//     fB = select((0 == fHi),fL,
//                 (1 == fHi),fL,
//                 (2 == fHi),fN,
//                 (3 == fHi),fV,
//                 (4 == fHi),fV,
//                 (5 == fHi),fM,
//                 0.0f);

//     output(x,y,0,t) = fR;
//     output(x,y,1,t) = fG;
//     output(x,y,2,t) = fB;
//     return output;

// }

// Func angle2rgb (Func v) {
//     Var x, y, c, t;
//     Func ov, a;
//     ov(x,y,c,t) = cast <float> (0.0f);
//     Expr pi2(2*M_PI);
//     a(x,y,c,t) = v(x,y,c,t) / pi2;
//     ov(x,y,0,t) = a(x,y,c,t);
//     ov(x,y,1,t) = 1;
//     ov(x,y,2,t) = 1;
//     return ov;
// }

// Func outputvelocity(Func Blur, Func Speed, Func Angle, int border, Expr speedthreshold, Expr filterthreshold) {
//     extern Expr width;
//     extern Expr height;

//     Func Blur3, Speed3;
//     Blur3(x,y,c,t) = cast <float> (0.0f);
//     Speed3(x,y,c,t) = cast <float> (0.0f);

// //Scale the grey level images
//     Blur(x,y,0,t) = (Blur(x,y,0,t) - minimum(Blur(x,y,0,t))) / (maximum(Blur(x,y,0,t)) - minimum(Blur(x,y,0,t)));
//     //Concatenation along the third dimension
//     Blur3(x,y,0,t) = Blur(x,y,0,t);
//     Blur3(x,y,1,t) = Blur(x,y,0,t);
//     Blur3(x,y,2,t) = Blur(x,y,0,t);

// //Speed scaled to 1
//     //Concatenation along the third dimension
//     Speed3(x,y,1,t) = Speed(x,y,0,t);
//     Speed3(x,y,2,t) = Speed(x,y,0,t);

// //Use the log speed to visualise speed
//     Func LogSpeed;
//     LogSpeed(x,y,c,t) = fast_log(Speed3(x,y,c,t) + Expr(0.0000001f))/fast_log(Expr(10.0f));
//     LogSpeed(x,y,c,t) = (LogSpeed(x,y,c,t) - minimum(LogSpeed(x,y,c,t))) / (maximum(LogSpeed(x,y,c,t)) - minimum(LogSpeed(x,y,c,t)));

// //Make a colour image
//     // uint16_t rows = height;
//     // uint16_t cols = width;
//     // int depth = Angle.channels();

// //Do it the HSV way
//     Func colorImage;
//     colorImage(x,y,0,t) = Angle(x,y,0,t);

// //Do hsv to rgb
//     Func colorImage1;
//     colorImage1 = hsv2rgb(colorImage);

// // Assume the border equals to the size of spatial filter
// //Make the border
//     // int bir = rows + 2 * border;
//     // int bic = cols + 2 * border;
//     Expr orows = height / Expr(2);
//     Expr ocols = width / Expr(2);

// //Rotation matrix
//     int ph = 0;
//     Func mb, sb;
//    // if (rx < border - 1 || rx >= rows+border -1 || ry < border - 1 || ry >= cols+border - 1) {
//     Expr co1 = x - orows;
//     Expr co2 = - (y - ocols);
//     Expr cosPh(cos(ph));
//     Expr sinPh(sin(ph));
//     Expr rco1 = cosPh * co1 - sinPh * co2; //Using rotation matrix
//     Expr rco2 = sinPh * co1 + cosPh * co2;
//     // Expr justPi (M_PI);
//     mb(x,y,c,t) =
//         select (((x < (border - 1)) ||
//                   (x >= (height+border -1)) ||
//                   (y < (border - 1)) ||
//                   (y >= (width+border - 1))),
//                 atan2(rco1,rco2) + Expr(M_PI),mb(x,y,c,t));
//     sb(x,y,c,t) =
//          select (((x < (border - 1)) ||
//                   (x >= (height+border -1)) ||
//                   (y < (border - 1) ) ||
//                   (y >= (width+border - 1))),
//                   1, sb (x,y,c,t));

//     Func cb;
//     cb = angle2rgb(mb);

// //Get the old data
//     // Expr pi2(2*M_PI);
//     colorImage1(x,y,0,t)=colorImage(x,y,0,t) * Expr(2*M_PI);
//     colorImage1=angle2rgb(colorImage1);
//     colorImage1(x,y,c,t)=select(abs(Speed3(x,y,c,t))<speedthreshold,Expr(0.0f),colorImage1(x,y,c,t));
//     Func colorImage2;
//     colorImage2(x,y,c,t) = colorImage1(x,y,c,t) * Speed(x,y,c,t);

// //Put the data in the border
//     RDom bordx (border,rows + border);
//     RDom bordy (border,cols + border);
//     Func ang1, ang2;
//     ang1 (x,y,c,t) = cast <float> (0.0f);
//     ang2 (x,y,c,t) = cast <float> (0.0f);

//     cb(bordx, bordy,c,t) = colorImage1(x,y,c,t);
//     ang1 = cb;
//     cb(bordx, bordy,c,t) = colorImage2(x,y,c,t);
//     ang2 = cb;
//     sb(bordx, bordy,c,t) = Speed3(x,y,c,t);
//     Speed3 = sb;
//     sb(bordx, bordy,c,t) = Blur3(x,y,c,t);
//     Blur3 = sb;

//     // Func I;
//     // I (x,y,c,t) = Blur3(x,y,c,t) + Speed3(x,y - height,c,t) + ang1(x - width,y,c,t) + ang2(x - width,y - height,c,t);
//     //I = cat(2,cat(1,Blur,Speed),cat(1,ang1,ang2));
//     return I;
// }
}
Func opticalFlow_estimate (Func stBasis, uint8_t nAngle, uint8_t * orders, \
                           Expr filterthreshold, Expr divisionthreshold,\
                           Expr divisionthreshold2) {
// This function estimates components of optical flow fields as well as its speed and direction of movement
// basis: from spatio-temporal filters, angle: number of considered angles
    // orders: x (spatial index), y ( spatial index ), t (time index) and s ?
// ColorMgather function in MATLAB
// Pipeline:
// 1. Compute oriented filter basis at a particular angle
// {
//     basis -> X
//           -> Y
//           -> T
//           -> Xrg
//           -> Yrg
//           -> Trg
//           -> Xk
//           -> Yk
//           -> Tk
// }
    // Expr M0,M1,M2,M3,temp_a,temp_b,Tkd;
    // Expr M0 = cast<float>(0.0f);
    // Expr M1 = cast<float>(0.0f);
    // Expr M2 = cast<float>(0.0f);
    // Expr M3 = cast<float>(0.0f);
    // Expr temp_a = cast<float>(0.0f);
    // Expr temp_b = cast<float>(0.0f);
    // Expr Tkd = cast<float>(0.0f);

    // Expr D0 = cast<float>(0.0f);
    // Expr D1 = cast<float>(0.0f);
    // Expr D2 = cast<float>(0.0f);
    // Expr D3 = cast<float>(0.0f);
    // Expr N0 = cast<float>(0.0f);
    // Expr N1 = cast<float>(0.0f);
    // Expr N2 = cast<float>(0.0f);
    // Expr N3 = cast<float>(0.0f);
    // Expr A0 = cast<float>(0.0f);
    // Expr A1 = cast<float>(0.0f);
    // Expr A2 = cast<float>(0.0f);
    // Expr A3 = cast<float>(0.0f);
    Func Tkd;
    Func fA0; fA0(x,y,t) = Expr(0.0f);
    Func fA1; fA1(x,y,t) = Expr(0.0f);
    Func fA2; fA2(x,y,t) = Expr(0.0f);
    Func fA3; fA3(x,y,t) = Expr(0.0f);
    Func fD0; fD0(x,y,t) = Expr(0.0f);
    Func fD1; fD1(x,y,t) = Expr(0.0f);
    Func fD2; fD2(x,y,t) = Expr(0.0f);
    Func fD3; fD3(x,y,t) = Expr(0.0f);
    Func fN0; fN0(x,y,t) = Expr(0.0f);
    Func fN1; fN1(x,y,t) = Expr(0.0f);
    Func fN2; fN2(x,y,t) = Expr(0.0f);
    Func fN3; fN3(x,y,t) = Expr(0.0f);

    // std::vector<Expr> basisAtAngleExpr(5,cast<float>(0.0f));
    // Tuple basisAtAngle = Tuple(basisAtAngleExpr);
    Func basisAtAngle[nAngle/2];

// Compute spatial-temporal basis
    for (int iA = 0; iA <= nAngle / 2 - 1; iA++) {
        float aAngle = 2*iA*M_PI/nAngle;
        basisAtAngle[iA] = ColorMgather(stBasis, aAngle, orders, filterthreshold,
                                        divisionthreshold, divisionthreshold2);
        basisAtAngle[iA].compute_root();

        Expr M0,M1,M2,M3;
        M0 = basisAtAngle[iA](x,y,t)[0]; M1 = basisAtAngle[iA](x,y,t)[1];
        M2 = basisAtAngle[iA](x,y,t)[2]; M3 = basisAtAngle[iA](x,y,t)[3];

        fD0(x,y,t) += M0 * M0;
        fD1(x,y,t) += M0 * M2;
        fD2(x,y,t) += M2 * M0;
        fD3(x,y,t) += M2 * M2;
        fN0(x,y,t) += M1 * M0;
        fN1(x,y,t) += M1 * M2;
        fN2(x,y,t) += M3 * M0;
        fN3(x,y,t) += M3 * M2;
        // temp_a = abs(M0) * M1;
        // temp_b = abs(M2) * M3;
        Expr cosaAngle((float) cos(aAngle));
        Expr sinaAngle((float) sin(aAngle));
        // A0  += temp_a  * cosaAngle;
        // A1  += temp_b  * sinaAngle;
        // A2  += temp_a  * sinaAngle;
        // A3  += temp_b  * cosaAngle;
        // fA0(x,y,c,t) += abs(M0) * M1 * cosaAngle;
        fA0(x,y,t) += abs(M0) * M1 * cosaAngle;
        fA1(x,y,t) += abs(M2) * M3 * sinaAngle;
        fA2(x,y,t) += abs(M0) * M1 * sinaAngle;
        fA3(x,y,t) += abs(M2) * M3 * cosaAngle;
    }

    Tkd(x,y,t) = basisAtAngle[nAngle/2-1](x,y,t)[4];
    // return basisAtAngle[0]; // 4 debug
    // // Schedule basis
    // for (int iA = 0; iA <= nAngle / 2 - 1; iA++) {
    //     fD0.update(iA);
    //     fD1.update(iA);
    //     fD2.update(iA);
    //     fD3.update(iA);
    //     fN0.update(iA);
    //     fN1.update(iA);
    //     fN2.update(iA);
    //     fN3.update(iA);
    //     fA0.update(iA);
    //     fA1.update(iA);
    //     fA2.update(iA);
    //     fA3.update(iA);
    // }

    fD0.compute_root();
    fD1.compute_root();
    fD2.compute_root();
    fD3.compute_root();
    fN0.compute_root();
    fN1.compute_root();
    fN2.compute_root();
    fN3.compute_root();
    fA0.compute_root();
    fA1.compute_root();
    fA2.compute_root();
    fA3.compute_root();
    Tkd.compute_root();

    Func top_func;
    Func bottom_func;
    Expr speed0;
    Expr speed1;

// Polar fig ?
    top_func(x,y,t) = fN0(x,y,t) * fN3(x,y,t) - fN1(x,y,t) * fN2(x,y,t);
    top_func.compute_root();
    bottom_func(x,y,t) = fD0(x,y,t)  * fD3(x,y,t)  - fD1(x,y,t)  * fD2(x,y,t);
    bottom_func.compute_root();

    speed0 = sqrt(sqrt(abs(Mdefdiv(top_func(x,y,t) , bottom_func(x,y,t), Expr(0.0f)))));
    speed1 = Manglecalc(fA0(x,y,t) , fA1(x,y,t) , fA2(x,y,t) , fA3(x,y,t));
// Display the results ?
    speed0 = select(abs(Tkd(x,y,t)) > filterthreshold,speed0,Expr(0.0f));
    speed1 = select(abs(Tkd(x,y,t)) > filterthreshold,speed1,Expr(0.0f));

    Func speed("speed"); speed(x,y,t) = Tuple(speed0,speed1);
    return speed;

//Bimg = [T0n speed0 speed1] ?
//    Func img = outputvelocity(T0n,Func(speed0),Func(speed1),16, speedthreshold, filterthreshold);
}
