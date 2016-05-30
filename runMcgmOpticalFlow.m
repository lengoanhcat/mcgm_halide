%% This script is written to replicate Alan New Model for computing
%% optical flow with Multichannel Gaussian Model

%% (re)Compile mcgmOpticalFlow.mexa64 with Halide

% Add the path to mex_halide.m.
setenv('HALIDE_PATH','/usr/local/halide');
addpath(fullfile(getenv('HALIDE_PATH'), 'tools'));

% Build the mex library from the blur generator.
% setenv('HL_DEBUG_CODEGEN','1');
mex_halide('./mcgmOpticalFlow.cpp');

% Add path for displaying color result
addpath('/home/lpzal1/Alan_new_model');

close all
% bufferData: store frame data of predefined temporal buffer
% frmHeight: frame height, frmWidth: frame width,
% noColorChan: number of color channels, noFrm: number of frames
vidDirName = 'Stimulus/'; noFrm = 50;

frmWidth = 320;
frmHeight = 240;
noColorChan = 3;
bufferSize = 23;
bufferData = single(zeros(frmHeight,frmWidth,noColorChan,bufferSize));
speedData0 = single(zeros(frmHeight,frmWidth,noColorChan,bufferSize));
speedData1 = single(zeros(frmHeight,frmWidth,noColorChan,bufferSize));
% tic
figure(1);
%% Read data into the buffer
for iFrm = 1:noFrm % run for noFrm
    %% Read data from each of first (offset + bufferSize - 1) frames
    frameName = [vidDirName, num2str(iFrm,'newCR3%.3d'), '.png'];
    frameData = im2single(imread(frameName));    
    deriData = ColorDeri(frameData);
    if iFrm < bufferSize
        bufferData(:,:,:,iFrm) = frameData;
        continue;
    else
        bufferData = circshift(bufferData,1,4);
        bufferData(:,:,:,1) = frameData;
    end

    % Do nothing just a funny test for now    
    mcgmOpticalFlow(bufferData,5,3,speedData0,speedData1);
    subplot(1,2,1); imshow(speedData0(:,:,:,ceil(bufferSize/2)));
    subplot(1,2,2); imshow(speedData0(:,:,:,ceil(bufferSize/2)));
    print('test');
% $$$     img=outputvelocity(T0n,speed0,speed1,16, speedthreshold, ...
% $$$                        filterthreshold);    
end
% toc
%% Creating a video to verify that it has worked

% aviobj = VideoWriter ('test1.avi');
% open(aviobj)
% 
% figure()
% 
% for k = 1:bufferSize
%     imshow(bufferData(:,:,:,k))
%     frame = getframe;
%     writeVideo(aviobj,frame);
% end
% close (aviobj);
