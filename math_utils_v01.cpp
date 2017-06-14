int factorial(int n)
{
    // return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    int factorial_n = 1;
    for (int i = 2; i <= n; i++)
        factorial_n = factorial_n * i;

    return factorial_n;
}

int combination(int v, int k) {
// caculate combinations of v elements when choose k.
    if (v < k) return 0;
    else if (v == k) return 1;
    else return factorial(v)/(factorial(v-k)*factorial(k));
}

Func cross(Tuple a, Tuple b) {
//compute cross product of two vectors
    std::vector<Expr> cp_expr(3,cast<float>(0.0f));

    // Compute a cross-product of two vector
    cp_expr[0] = a[1]*b[2] - a[2]*b[1];
    cp_expr[1] = a[2]*b[0] - a[0]*b[2];
    cp_expr[2] = a[0]*b[1] - a[1]*b[0];

    Func cp_func; cp_func(x,y,t) = Tuple(cp_expr);

    return cp_func;
}

Func dot(Tuple a, Tuple b) {
//compute dot product of two vectors
    Expr dp("dp"); dp = cast<float>(0.0f);
    assert(a.size() == b.size());
    for (uint8_t iC = 0; iC < a.size(); iC++)
        dp += a[iC]*b[iC];

    Func dp_func; dp_func(x,y,t) = dp;
    return dp_func;
}

Expr Mdefdiv(Expr A, Expr B, Expr divth) {
// Avoid dividing by 0
    return select(abs(B)>divth,A/B,Expr(0.0f));
}

Expr Mdefdivang(Expr A, Expr B, Expr R, Expr divth) {
// Avoid dividing by 0
    // return select(abs(A/dot)>divth,A/B,Expr(0.0f));
    return select(R<divth,A/B,Expr(0.0f));
}

Expr Manglecalc(Expr A0, Expr A1, Expr A2, Expr A3) {
    // Warping angles into region [0-360]
    // Expr angle = cast<float>(0.0f);
    // angle = (atan2((A2 + A3), (A0 - A1)) * cast<uint8_t>(180) ) / Expr(M_PI);
    // angle = select  (angle >= Expr(360.0f), angle - 360.0f,
    //                  angle <= Expr(0.0f), angle + 360.0f,
    //                  angle);
    // return angle;

    return atan2((A2 + A3), (A0 - A1));
}
