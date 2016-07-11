Func color_derivative(Func input) {

        float m = 1.0f/3.0f;

        Func specDeri("specDeri"); //Multiply the image with the ColorRF matrix (see original Matlab code)
        // specDeri (x, y, c, t) = cast <float> (0.0f);
        // specDeri (x, y, 0, t) = input(x,y,0,t) * m       + input(x,y,1,t) * m      + input(x,y,2,t) * m;
        // specDeri (x, y, 1, t) = input(x,y,0,t) * 0.25f   + input(x,y,1,t) * 0.25f  - input(x,y,2,t) * 0.5f;
        // specDeri (x, y, 2, t) = input(x,y,0,t) * 0.5f    - input(x,y,1,t) * 0.5f;

        Image<float> colorRF(3,3,"colorRF");
        colorRF(0,0) = m; colorRF(1,0) = m; colorRF(2,0) = m;
        colorRF(0,1) = 0.25f; colorRF(1,1) = 0.25f; colorRF(2,1) = -0.5f;
        colorRF(0,2) = 0.5f; colorRF(1,2) = -0.5f; colorRF(2,2) = 0.0f;

        RDom rf_x(0,3);
        specDeri(x,y,c,t) = sum(rf_x,input(x,y,rf_x,t)*colorRF(rf_x,c),"sum_colorRF");

        return specDeri;
}
