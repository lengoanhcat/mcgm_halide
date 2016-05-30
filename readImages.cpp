#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <Halide.h>
#include <halide_image_io.h>
#include <ctime>
#include <FreeImagePlus.h>
#include <stdio.h>
#include <cstdlib>
#include <stdio.h>
#include <math.h>

using namespace Halide;
using namespace Halide::BoundaryConditions;
using namespace Halide::Tools;
using namespace std;

unsigned int bufferSize = 23;
Halide::Image<float> images[23]; //Buffer of all the frames that we want

Halide::Image<float> * readImages (int noFrames);
Halide::Image<float> * readImages (int noFrames) { //noFrames is the number of the last frame that we want to extract
    
    if (noFrames < bufferSize || noFrames > 200)
    {
        printf("The number of the last frame has to be greater than the buffer size (here %d) and less than the total number of frames (here 200)", bufferSize);
        return NULL;
    }
    
    else   
    {   //Reading data into the buffer a buffer called images
        for (int iFrm = noFrames - bufferSize; iFrm < noFrames; iFrm++) {
            char buffer[50];
            sprintf(buffer,"./Stimulus/newCR3%.3d.png",iFrm + 1);
            
            //Circular shift + Reading the image and storing it in the buffer
            images [iFrm-noFrames+bufferSize] = load_image(buffer);
        }
        return images;
    }
}

//Main to test

// int main(int argc, char **argv) {
//     int lastFrame = 150;
//     time_t tstart, tend;
//     tstart = time(0);
//     Halide::Image<float> * test = readImages (lastFrame);
//     tend = time(0);
//     cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;
//     for (int i = 0; i < bufferSize; i++)
//     {
//         char buffer[60];
//         sprintf(buffer,"./testDeriResults/Result%.3d.png", lastFrame - bufferSize + i + 1);
//         save_image(test[i],buffer);
//     }
// }

