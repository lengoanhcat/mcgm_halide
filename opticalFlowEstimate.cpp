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

Expr ColorMgetfilter(Func basis, float angle, uint8_t iXo, uint8_t iYo, uint8_t iTo, uint8_t iCo ) {
    // Compute a rotated basis at (iXo,iYo,iTo,iCo) order with angle value

    // temporary setting
    uint8_t numSTB = 63;
    uint8_t numSB = 21;

    angle = -1*angle - M_PI/2;
    float * weights;

    Expr work; // work: rotated basis at a particular spatio-temporal order
    work = undef<float>();

    weights = (float *) calloc(iXo+iYo+1,sizeof(float));

    // compute weights for possible orders
    for (int i = 0; i <= iXo; i++)
        for (int j = 0; j <= iYo; j++)
            weights[iXo+iYo-i-j] += combination(iXo,iYo)*pow((-1.0f),i)*pow(cos(angle),(iXo-i+j))*pow(sin(angle),(iYo+i-j));

    // get filtered expression at paricular order and angle value
    for (int k=0; k<=(iXo+iYo); k++) {
        int index = Mgetfilterindex(iXo+iYo-k,k,iTo,numSTB,numSB);
        if ((index > 0) && (weights[iXo+iYo-k] != 0))
            work += weights[iXo+iYo-k]*basis(x,y,c,t)[index];
    }

    free(weights);
    return work;
}

Tuple ColorMgather(Func basis, float angle, uint8_t * orders, Expr filterthreshold, Expr divisionthreshold, Expr divisionthreshold2) {
    uint8_t x_order = orders[0];
    uint8_t y_order = orders[1];
    uint8_t t_order = orders[2];
    uint8_t c_order = orders[3];

    Expr X("X"),Y("Y"),T("T"),Xrg("Xrg"),Yrg("Yrg"),Trg("Trg");
    uint8_t max_order = x_order + 1;
    std::vector<Expr>Xk_expr (max_order,undef<float>());
    std::vector<Expr>Yk_expr (max_order,undef<float>());
    std::vector<Expr>Tk_expr (max_order,undef<float>());

    // Expr Xk[max_order],Yk[max_order],Tk[max_order];

    int k = 0;
    for (int iXo = 0; iXo < x_order; iXo++) // x_order
        for (int iYo = 0; iYo < y_order; iYo++) // y_oder
            for (int iTo = 0; iTo < t_order; iTo++) // t_order
                for (int iCo = 0; iCo < c_order; iCo ++ ) // c_order: index of color channel
                {
                    if ((iYo+iTo+iXo == 0 || iYo+iTo+iXo == 1) && ((iXo+iYo+iTo+1) < (x_order + 1))) {
                        X = ColorMgetfilter(basis, angle, iXo+1, iYo, iTo, iCo);
                        Y = ColorMgetfilter(basis, angle, iXo, iYo+1, iTo, iCo);
                        T = ColorMgetfilter(basis, angle, iXo, iYo, iTo+1, iCo);

                        Xrg = ColorMgetfilter(basis, angle, iXo+1, iYo, iTo, iCo+1);
                        Yrg = ColorMgetfilter(basis, angle, iXo, iYo+1, iTo, iCo+1);
                        Trg = ColorMgetfilter(basis, angle, iXo, iYo, iTo+1, iCo+1);

                        k = iXo + iYo + iTo + iCo;
                        Xk_expr[k] += X + Xrg;
                        Yk_expr[k] += Y + Yrg;
                        Tk_expr[k] += T + Trg;
                    }
                }

    Tuple Xk = Tuple(Xk_expr);
    Tuple Yk = Tuple(Yk_expr);
    Tuple Tk = Tuple(Tk_expr);

    std::vector<Expr> st_expr(6,undef<float>());
    for (int iK=0; iK < k; iK++) {
        st_expr[0] += Xk_expr[k]*Tk_expr[k];
        st_expr[1] += Tk_expr[iK]*Tk_expr[iK];
        st_expr[2] += Xk_expr[iK]*Xk_expr[iK];
        st_expr[3] += Yk_expr[iK]*Tk_expr[iK];
        st_expr[4] += Yk_expr[iK]*Yk_expr[iK];
        st_expr[5] += Xk_expr[iK]*Yk_expr[iK];
    }

    Func st("st"); st(x,y,c,t) = Tuple(st_expr);

    // Generate mean filter
    // float win = 7.0;
    // Func generateMeanFilter("meanfilter");
    // generateMeanFilter(x,y) = 1.0f/(win*win);
    // Image<float> meanfilter = generateMeanFilter.realize(win,win);
    float win = 7.0;
    Image<float> meanfilter(7,7,"meanfilter_data");
    meanfilter(x,y) = Expr(1.0f/(win*win));

    RDom rMF(meanfilter);

    for (uint8_t iPc=0; iPc<6; iPc++) // iPc: index of product component
        st_expr[iPc] = sum(rMF,meanfilter(rMF.x,rMF.y)*                 \
                           st(x + rMF.x - 1,y + rMF.y - 1,c,t)[iPc],"mean_filter");

    Tuple st_tuple = Tuple(st_expr);

    Tuple pbx = Tuple(st_tuple[2],st_tuple[5],st_tuple[0]);
    Tuple pby = Tuple(st_tuple[5],st_tuple[4],st_tuple[3]);
    Tuple pbt = Tuple(st_tuple[0],st_tuple[3],st_tuple[1]);

    Tuple pbxy = cross(pby,pbx);
    Tuple pbxt = cross(pbx,pbt);
    Tuple pbyt = cross(pby,pbt);

    Expr pbxyd("pbxyd"); pbxyd = dot(pby,pbx);
    Expr pbxtd("pbxtd"); pbxtd = dot(pbx,pbt);
    Expr pbytd("pbytd"); pbytd = dot(pby,pbt);

    Expr yt_xy("yt_xy"); yt_xy = dot(pbyt,pbxy);
    Expr xt_yt("xt_yt"); xt_yt = dot(pbxt,pbyt);
    Expr xt_xy("xt_xy"); xt_xy = dot(pbxt,pbxy);
    Expr yt_yt("yt_yt"); yt_yt = dot(pbyt,pbyt);
    Expr xt_xt("xt_xt"); xt_xt = dot(pbxt,pbxt);
    Expr xy_xy("xy_xy"); xy_xy = dot(pbxy,pbxy);

    Expr Tkd("Tkd"); Tkd = dot(Tk,Tk);

    // Expr Dimen = pbxyd/xy_xy;
    Expr kill(1.0f);

    Expr Oxy = Mdefdiv(st_tuple[5] - Mdefdivang(yt_xy,yt_yt,pbxyd,divisionthreshold2)*st_tuple[3]*kill,st_tuple[4],divisionthreshold);
    Expr Oyx = Mdefdiv(st_tuple[5] - Mdefdivang(xt_xy,xt_xt,pbxyd,divisionthreshold2)*st_tuple[0]*kill,st_tuple[2],divisionthreshold);
    Expr C0 = st_tuple[3] * Mdefdivang(Expr(-1.0f)*xt_yt,yt_yt,pbxyd,divisionthreshold2)*kill;
    Expr M0 = Mdefdiv(st_tuple[0] + C0, st_tuple[1]*pow(Mdefdivang(xt_yt,yt_yt,pbxyd,divisionthreshold2),2),divisionthreshold);

    Expr C1 = st_tuple[5] * Mdefdivang(Expr(-1.0f)*xt_xy,xy_xy,pbxyd,divisionthreshold2)*kill;
    Expr P1 = pow(Mdefdivang(xt_yt,xt_xt,pbxyd,divisionthreshold2),2)*kill + 1;
    Expr Q1 = st_tuple[2] * (pow(Oxy,2)+1);
    Expr M1 = Mdefdiv(((st_tuple[0]-C1)*P1),Q1,divisionthreshold);

    Expr C2 = st_tuple[0] * Mdefdivang(Expr(-1.0f)*xt_yt,xt_xt,pbxyd,divisionthreshold2)*kill;
    Expr M2 = Mdefdiv(st_tuple[3]+C2,st_tuple[1]*(pow(Mdefdivang(xt_yt,xt_xt,pbxyd,divisionthreshold2),2)*kill+1),divisionthreshold);

    Expr C3 = st_tuple[5] * Mdefdivang(yt_xy,xy_xy,pbxyd,divisionthreshold2)*kill;
    Expr P3 = pow(Mdefdivang(xt_yt,yt_yt,pbxyd,divisionthreshold2),2)*kill + 1;
    Expr Q3 = st_tuple[4] * (pow(Oxy,2)+1);
    Expr M3 = Mdefdiv(((st_tuple[3]-C3)*P3),Q3,divisionthreshold);

    return {M0,M1,M2,M3,Tkd};
}

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

Func opticalFlow_estimate(Func basis,uint8_t nAngle, uint8_t * orders, \
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

    Expr M0,M1,M2,M3,temp_a,temp_b,Tkd;
    Expr D0 = cast<float>(0.0f);
    Expr D1 = cast<float>(0.0f);
    Expr D2 = cast<float>(0.0f);
    Expr D3 = cast<float>(0.0f);
    Expr N0 = cast<float>(0.0f);
    Expr N1 = cast<float>(0.0f);
    Expr N2 = cast<float>(0.0f);
    Expr N3 = cast<float>(0.0f);
    Expr A0 = cast<float>(0.0f);
    Expr A1 = cast<float>(0.0f);
    Expr A2 = cast<float>(0.0f);
    Expr A3 = cast<float>(0.0f);

    Expr top, bottom, speed0, speed1;

    std::vector<Expr> basisAtAngleExpr(5,undef<float>());
    Tuple basisAtAngle = Tuple(basisAtAngleExpr);

    for (int iA = 0; iA <= nAngle / 2 - 1; iA++) {
        float aAngle = 2*iA*M_PI/nAngle;
        basisAtAngle = ColorMgather(basis, aAngle, orders, filterthreshold,
                                    divisionthreshold, divisionthreshold2);

        M0 = basisAtAngle[0]; M1 = basisAtAngle[1]; M2 = basisAtAngle[2];
        M3 = basisAtAngle[3]; Tkd = basisAtAngle[4];

        D0 += M0 * M0;
        D1 += M0 * M2;
        D2 += M2 * M0;
        D3 += M2 * M2;
        N0 += M1 * M0;
        N1 += M1 * M2;
        N2 += M3 * M0;
        N3 += M3 * M2;
        temp_a = abs(M0) * M1 ;
        temp_b  = abs(M2 ) * M3 ;
        Expr cosaAngle(cos(aAngle));
        Expr sinaAngle(sin(aAngle));
        A0  += temp_a  * cosaAngle;
        A1  += temp_b  * sinaAngle;
        A2  += temp_a  * sinaAngle;
        A3  += temp_b  * cosaAngle;
    }

//Polar fig ?

    top  = N0  * N3  - N1  * N2 ;
    bottom  = D0  * D3  - D1  * D2 ;
    speed0  = sqrt(sqrt(abs(Mdefdiv(top , bottom , Expr(0.0f)))));
    speed1  = Manglecalc(A0 , A1 , A2 , A3);
//Display the results ?

    speed0 = select(abs(Tkd) > filterthreshold,speed0,Expr(0.0f));

    // Between two lines for filtering speed 1, the first line generates a lot
    // of expression during code generation while the second allows neat compilation.
    // It is not clear how each line affect the accuracy.
    speed1 = select(abs(Tkd) > filterthreshold,speed0,Expr(0.0f));

    Func speed("speed"); speed(x,y,c,t) = Tuple(speed0,speed1);

    return speed;

//Bimg = [T0n speed0 speed1] ?
//    Func img = outputvelocity(T0n,Func(speed0),Func(speed1),16, speedthreshold, filterthreshold);
}
