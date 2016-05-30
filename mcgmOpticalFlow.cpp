#include <Halide.h>
using namespace Halide;

#include <stdio.h>
#include <math.h>
// Include some support code for loading pngs.
// #include <halide_image_io.h>
// using namespace Halide::Tools;

// #define M_PI 3.141592
// Declare function variable
Var x("x"),y("y"),c("c"),t("t"); // x,y: spatial indices, t and c are time and color indices
Expr width, height;
uint8_t sn;

#include "math_utils.cpp"
#include "colorDerivative.cpp"
#include "temporalDerivative.cpp"
#include "spatialDerivative.cpp"
#include "opticalFlowEstimate.cpp"

class McgmOpticalFlow : public Generator<McgmOpticalFlow> {
public:
    // This is the input image: a 4D sequence of (color) images with float32 pixels.
    ImageParam input{Float(32), 4, "input"};
    // The filter coefficient, alpha is the weight of the input to the
    // filter.
    Param<uint8_t> spaOrd{"spatialOrder"};
    Param<uint8_t> temOrd{"temporalOrder"};
    // Param<uint8_t> rotAng{"rotatingAngle"};
    // Param<uint8_t> spatialFilterSize{"spatialFilterSize"};
    // Param<uint8_t *> orders{"listOfOrders"};

    Func build() {

        Expr width = input.width();
        Expr height = input.height();

        uint8_t orders[4] = {5,2,2,2};
        uint8_t angle = 24;
        Expr filterthreshold((float) pow(10.0,-4.8));
        Expr divisionthreshold((float) pow(10.0,-30));
        Expr speedthreshold(0.000001f);

        // Our input is an ImageParam, but blur_cols takes a Func, so
        // we define a trivial func to wrap the input.
        Func input_func("input_func");
        input_func = BoundaryConditions::repeat_edge(input);

        // Scheduling is done inside each function separately
        // First, convert RGB channel to color derivative
        Func d_cspec("d_cspec");
        // d_cspec.trace_stores();
        d_cspec = color_derivative(input_func);

        // Do (0rd, 1st, 2nd) temporal filtering
        Func Tx("Tx"); Tx = temporal_derivative(d_cspec,temOrd);

        // Do spatial filtering with differential gaussian kernel up to 5-th order
        Func basis("basis"); basis = spatial_derivative(Tx,spaOrd,temOrd);

        // Compute velocity
        Func optFlw("optFlw");
        optFlw = opticalFlow_estimate(basis,angle,orders,filterthreshold,divisionthreshold,divisionthreshold);

        // std::vector<Argument> args;
        // T.compile_to_lowered_stmt("temporal_derivative.html", args, StmtOutputFormat::HTML);
        // T.compile_to_c("temporal_derivative_generated_by_halide.cpp",args,"temporal_derivative");


        // Computer angle and speed with speed model
         return optFlw;
    }
};

auto mcgmOpticalFlow = RegisterGenerator<McgmOpticalFlow>("McgmOpticalFlow");
