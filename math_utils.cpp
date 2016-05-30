int factorial(int n)
{
    // return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
    int factorial_n = 1;
    for (int i = 1; i <= n; i++)
        factorial_n = factorial_n * i;

    return factorial_n;
}

int combination(int v, int k) {
// caculate combinations of v elements when choose k.
    if (v < k) return 0;
    else if (v == k) return 1;
    else return factorial(v)/(factorial(v-k)*factorial(k));
}

Tuple cross(Tuple a, Tuple b) {
//compute cross product of two vectors
    std::vector<Expr> cp_expr(3,undef<float>());

    // Compute a cross-product of two vector
    cp_expr[0] = a[1]*b[2] - a[2]*b[1];
    cp_expr[1] = a[2]*b[0] - a[0]*b[2];
    cp_expr[2] = a[0]*b[1] - a[1]*b[0];

    Tuple cp_tuple = Tuple(cp_expr);

    return cp_tuple;
}

Expr dot(Tuple a, Tuple b) {
//compute dot product of two vectors
    Expr dp("dp"); dp = undef<float>();
    assert(a.size() == b.size());
    for (uint8_t iC = 0; iC < a.size(); iC++)
        dp += a[iC]*b[iC];
    return dp;
}

Expr Mdefdiv(Expr A, Expr B, Expr divth) {
// Avoid dividing by 0
    return select(abs(B)>divth,A/B,0.0f);
}

Expr Mdefdivang(Expr A, Expr B, Expr dot, Expr divth) {
// Avoid dividing by 0
    return select(abs(A/dot)>divth,A/B,0.0f);
}

Expr Manglecalc(Expr A0, Expr A1, Expr A2, Expr A3) {
    // Warping angles into region [0-360]
    Expr angle = (atan2((A2 + A3), (A0 - A1)) * 180.0f) / Expr(M_PI);
    angle = select  (angle >= Expr(360.0f), angle - 360.0f,
                     angle <= Expr(0.0f), angle + 360.0f,
                     angle);

    return angle;
}
