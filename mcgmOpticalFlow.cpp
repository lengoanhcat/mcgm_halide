#include <Halide.h>
using namespace Halide;

#include <stdio.h>
#include <math.h>
// #include <mex.h>
// Include some support code for loading pngs.
// #include <halide_image_io.h>
// using namespace Halide::Tools;

// #define M_PI 3.141592
// Declare function variable
Var x("x"),y("y"),c("c"),t("t"),x_outer("x_outer"),
    y_outer("y_outer"),x_inner("x_inner"),y_inner("y_inner"),
    tile_index("tile_index"); // x,y: spatial indices, t and c are time and color indices
Expr width, height, noChn, noFrm;
uint8_t sn;

#include "halide_profile_utils.cpp"
#include "math_utils.cpp"
#include "colorDerivative.cpp"
#include "temporalDerivative.cpp"
#include "spatialDerivative.cpp"
// #include "spatialTemporalDerivative.cpp"
#include "opticalFlowEstimate.cpp"

namespace {

class McgmOpticalFlow : public Generator<McgmOpticalFlow> {
public:
    // This is the input image: a 4D sequence of (color) images with float32 pixels.
    ImageParam input{Float(32), 4, "input"};
    // The filter coefficient, alpha is the weight of the input to the
    // filter.
    // Param<uint8_t> spaOrd{"spatialOrder"};
    // Param<uint8_t> temOrd{"temporalOrder"};
    // Param<uint8_t> rotAng{"rotatingAngle"};
    // Param<uint8_t> spatialFilterSize{"spatialFilterSize"};
    // Param<uint8_t *> orders{"listOfOrders"};
    Param<float> filterthreshold{"filterthreshold"};
    Param<float> divisionthreshold{"divisionthreshold"};
    Param<float> divisionthreshold2{"divisionthreshold2"};
    Param<float> speedthreshold{"speedthreshold"};

    Func build() {

        width = input.width();
        height = input.height();
        noChn = input.channels();
        noFrm = input.extent(3);

        uint8_t orders[4] = {5,2,2,2};
        uint8_t angle = 24; // 24;
        // Expr filterthreshold((float) pow(10.0,-4.8));
        // Expr divisionthreshold((float) pow(10.0,-30.0));
        // Expr divisionthreshold2((float) pow(10.0,-25.0));
        // Expr speedthreshold(0.000001f);

        // Our input is an ImageParam, but blur_cols takes a Func, so
        // we define a trivial func to wrap the input.
        Func input_func("input_func");
        // input_func = BoundaryConditions::repeat_edge(input);
        input_func(x,y,c,t) = BoundaryConditions::constant_exterior(input,0)(x,y,c,t);
        input_func
            .parallel(t)
            .vectorize(x,4)
            .compute_root();

        // Scheduling is done inside each function separately
        // First, convert RGB channel to color derivative
        Func d_cspec("d_cspec");
        d_cspec = color_derivative(input_func);

        // Best scheduling strategy for first 2 stage so far
        d_cspec
            .parallel(t)
            .vectorize(x,4)
            .compute_root();

        // Do (0rd, 1st, 2nd) temporal filtering
        Func Tx("Tx"); Tx = temporal_derivative(d_cspec);
        Tx
            .parallel(t)
            .vectorize(x,4)
            .compute_root();

        // Do spatial filtering with differential gaussian kernel up to 5-th order
        Func stBasis("stBasis"); stBasis = spatial_derivative(Tx);
        stBasis.set_custom_print(&my_print);
        stBasis
            .parallel(t)
            .vectorize(x,4)
            .compute_root();

        // Compute velocity
        Func optFlw("optFlw");
        optFlw = opticalFlow_estimate(stBasis,angle,orders,filterthreshold,divisionthreshold,divisionthreshold2);
        optFlw
            .parallel(t)
            .vectorize(x,4)
            .compute_root();
        // return optFlw; // 4 debug

        Func outPut("outPut");
        outPut(x,y,t) = Tuple(Tx(x,y,0,t)[0],optFlw(x,y,t)[0],optFlw(x,y,t)[1]);

        // Computer angle and speed with speed model
        return outPut;
    }
};

// auto mcgmOpticalFlow = RegisterGenerator<McgmOpticalFlow>("McgmOpticalFlow");
    Halide::RegisterGenerator<McgmOpticalFlow> register_my_gen{"mcgmOpticalFlow"};

}
