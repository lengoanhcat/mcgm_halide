#include "Halide.h"
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

#include "math_utils_v01.cpp"
#include "colorDerivative.cpp"
#include "temporalDerivative.cpp"
#include "spatialDerivative_autoscheduled.cpp"
// #include "spatialTemporalDerivative.cpp"
#include "opticalFlowEstimate_v02_autoscheduled.cpp" // using complex form of speed and inverse speed

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

        width = Expr((int) 960); // width = input.width();
        height = Expr((int) 1440); // height = input.height();
        noChn = Expr((int) 3); // noChn = input.channels();
        noFrm = Expr((int) 222);    // noFrm = input.extent(3);

        uint8_t orders[4] = {5,2,2,2};
        uint8_t angle = 24; // 24;
        // Expr filterthreshold((float) pow(10.0,-4.8));
        // Expr divisionthreshold((float) pow(10.0,-30.0));
        // Expr divisionthreshold2((float) pow(10.0,-25.0));
        // Expr speedthreshold(0.000001f);

        ////////////////////////////////////////////////////////////////
        // Our input is an ImageParam, but blur_cols takes a Func, so //
        // we define a trivial func to wrap the input.                //
        ////////////////////////////////////////////////////////////////
        Func input_func("input_func");
        // input_func = BoundaryConditions::repeat_edge(input);
        input_func(x,y,c,t) = BoundaryConditions::constant_exterior(input,0)(x,y,c,t);

        ////////////////////////////////////////////////////////
        // Scheduling is done inside each function separately //
        // First, convert RGB channel to color derivative     //
        ////////////////////////////////////////////////////////
        Func d_cspec("d_cspec");
        d_cspec = color_derivative(input_func);

        ///////////////////////////////////////////
        // Do (0rd, 1st, 2nd) temporal filtering //
        ///////////////////////////////////////////
        Func Tx("Tx"); Tx = temporal_derivative(d_cspec);

        /////////////////////////////////////////////////////////////////////////////
        // Do spatial filtering with differential gaussian kernel up to 5-th order //
        /////////////////////////////////////////////////////////////////////////////
        Func stBasis("stBasis"); stBasis = spatial_derivative(Tx);
        // stBasis.set_custom_print(&my_print);

        //////////////////////
        // Compute velocity //
        //////////////////////
        Func optFlw("optFlw");
        optFlw = opticalFlow_estimate(stBasis,angle,orders,filterthreshold,divisionthreshold,divisionthreshold2);

        //////////////////////////////////////
        // Auto-schedule the whole pipeline //
        //////////////////////////////////////

        int num_stencils = 2;
        std::vector<Func> stencils;
        for(int i = 0; i < num_stencils; i ++) {
            Func s("stencil_" + std::to_string(i));
            stencils.push_back(s);
        }

        for (int i = 0; i < num_stencils; i++) {
            stencils[i](x,y,t) = optFlw(x,y,t)[i];
            // Provide estimates on the output
            stencils[i].estimate(x,0,width).estimate(y,0,height).estimate(t,0,noFrm);
        }

        // Provide estimates on input
        input.dim(0).set_bounds_estimate(0,960);
        input.dim(1).set_bounds_estimate(0,1440);
        input.dim(2).set_bounds_estimate(0,3);
        input.dim(3).set_bounds_estimate(0,222);

        // Provide estimates on the parameters
        filterthreshold.set_estimate((float) pow(10.0,-4.8));
        divisionthreshold.set_estimate((float) pow(10.0,-30.0));
        divisionthreshold2.set_estimate((float) pow(10.0,-25.0));
        speedthreshold.set_estimate(0.000001f);

        std::vector<Target::Feature> lomond_features;
        lomond_features.push_back(Target::SSE41);
        lomond_features.push_back(Target::Matlab);
        lomond_features.push_back(Target::LargeBuffers);

        Target target;
        target.os = Target::Linux;
        target.arch = Target::X86;
        target.bits = 64;
        target.set_features(lomond_features);
        // Target target = get_target_from_environment();
        Pipeline pipeline(stencils);

        std::cout << "\n\n******************************************\nSCHEDULE:\n"
                  << "******************************************\n"
                  << pipeline.auto_schedule(target)
                  << "\n******************************************\n\n";

        Func outPut("outPut");
        outPut(x,y,t) = Tuple(stencils[0](x,y,t),stencils[1](x,y,t));

        // Computer angle and speed with speed model
        // return pipeline.get_func(84);
        return outPut;
    };
};

// auto mcgmOpticalFlow = RegisterGenerator<McgmOpticalFlow>("McgmOpticalFlow");
    Halide::RegisterGenerator<McgmOpticalFlow> register_my_gen{"mcgmOpticalFlow"};
}
