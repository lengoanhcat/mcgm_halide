%% This script is written to replicate Alan New Model for computing
%% optical flow with Multichannel Gaussian Model

% Add path for displaying color result
% addpath('/home/lpzal1/Alan_new_model');

close all
setenv('HL_NUM_THREADS','320');
addpath('/storage/cat/mexopencv');
% bufferData: store frame data of predefined temporal buffer
% frmHeight: frame height, frmWidth: frame width,
% noColorChan: number of color channels, noFrm: number of frames

% $$$ vidDirName = ['Stimulus/' testCase '/']; %  'Stimulus/plaid/'

% $$$ testCase = 'newCR3'; noFrm = 200;
% $$$ frmWidth = 320; frmHeight = 240; noColorChan = 3;

% testCase = 'bbaf2n_'; noFrm = 75;
% frmWidth = 360; frmHeight = 288; noColorChan = 3;
% testCase = 'ep0201f_'; noFrm = 44;
% frmWidth = 640; frmHeight = 480; noColorChan = 3;
% testCase = 'h2n2a_'; noFrm = 200;
% frmWidth = 1280/4; frmHeight = 960/3; noColorChan = 3;
% testCase = 'ns01t01_'; noFrm = 66;
% frmWidth = 1280; frmHeight = 720; noColorChan = 3;
% $$$ testCase = 'grating'; noFrm = 128;
% $$$ frmWidth = 128; frmHeight = 128; noColorChan = 3;


testCase = 'Cart_GRAPPA_2_PF_6_8_9fps'
% testCase = 'Cart_GRAPPA_3_PF_5_8_17fps'
% testCase = 'Cart_GRAPPA_3_PF_6_8_12fps' % could not run
% testCase = 'Cart_GRAPPA_4_PF_5_8_22fps'
% testCase = 'Cart_TPAT_2_PF_6_8_9fps'
% testCase = 'Cart_TPAT_3_PF_5_8_16fps'
% testCase = 'Cart_TPAT_3_PF_6_8_11fps'
% testCase = 'Cart_TPAT_4_PF_4_8_16fps'
% testCase = 'Cart_TPAT_4_PF_5_8_15fps'
% testCase = 'Cart_TPAT_4_PF_5_8_15fps_int'
% testCase = 'Cart_TPAT_4_PF_5_8_21fps'
% testCase = 'Cart_TPAT_4_PF_6_8_14fps'
% testCase = 'Cart_TPAT_5_PF_6_8_15fps'
% testCase = 'Cart_TPAT_3_PF_6_8_11fps';
inVidName = ['../FredMRI/' testCase '.mov'];
% frmWidth = 683; frmHeight = 858; noFrm = 100; noColorChan = 3;
% frmWidth = 320; frmHeight = 480; noColorChan = 3;

% inVidName = '../mcgm_halide/Video_stimulus/NINJA2_S001_S001_T001.MOV';

%% Interpolated by Gaussian Blur
% $$$ testCase = 'gratingInterp'; noFrm = 30;
% $$$ frmWidth = 128; frmHeight = 128; noColorChan = 1;
% $$$ testCase = 'newCR3Interp'; noFrm = 32;
% $$$ frmWidth = 320; frmHeight = 240; noColorChan = 1;

vidDirName = ['Stimulus/' testCase '/'];
outVidName = ['results/' testCase '_complexInterpretation_motion.avi'];
vidFileHdl = VideoWriter(outVidName);
vidFileHdl.FrameRate = 10;

% $$$ bufferSize = 23;
% $$$ bufferData = single(zeros(frmHeight,frmWidth,noColorChan,bufferSize));
% $$$ speedData0 = single(zeros(frmHeight,frmWidth,noColorChan,bufferSize));
% $$$ speedData1 = single(zeros(frmHeight,frmWidth,noColorChan, ...
% $$$                           bufferSize));
% $$$ speedData2 = single(zeros(frmHeight,frmWidth,noColorChan, ...
% $$$                           bufferSize));

%% Process streaming video online
% $$$ % tic
% $$$ figure(1);
% $$$ %% Read data into the buffer
% $$$ for iFrm = 1:noFrm % run for noFrm
% $$$     %% Read data from each of first (offset + bufferSize - 1) frames
% $$$     frameName = [vidDirName, num2str(iFrm,'newCR3%.3d'), '.png'];
% $$$     frameData = im2single(imread(frameName));    
% $$$     if iFrm < bufferSize
% $$$         bufferData(:,:,:,iFrm) = frameData;
% $$$         continue;
% $$$     else
% $$$         bufferData = circshift(bufferData,1,4);
% $$$         bufferData(:,:,:,1) = frameData;
% $$$     end
% $$$   
% $$$     mcgmOpticalFlow(bufferData,5,3,speedData0,speedData1);
% $$$     subplot(1,2,1); imshow(speedData0(:,:,:,ceil(bufferSize/2)));
% $$$     subplot(1,2,2); imshow(speedData1(:,:,:,ceil(bufferSize/2)));
% $$$     
% $$$     % print('test');
% $$$ % $$$     img=outputvelocity(T0n,speed0,speed1,16, speedthreshold, ...
% $$$ % $$$                        filterthreshold);    
% $$$ end

%% --------------------------------------------------------------------------------------------------

% $$$ %% Displaying with 4-D output
% $$$ inData = single(zeros(frmHeight,frmWidth,noColorChan,noFrm));
% $$$ o0 = single(zeros(frmHeight,frmWidth,noColorChan,noFrm));
% $$$ outVar = 'o0';
% $$$ numOut = 63;
% $$$ 
% $$$ %% IMPORTANT NOTE: create separate memory places for output
% $$$ %% variables to avoid overwriting on same places
% $$$ for iO = 1:numOut-1
% $$$     eval(['o',num2str(iO),'=single(zeros(frmHeight,frmWidth,noColorChan,noFrm));']);
% $$$     outVar = [outVar,',o',num2str(iO)];
% $$$ end
% $$$ 
% $$$ for iFrm = 1:noFrm
% $$$      %% Read data from each of first (offset + bufferSize - 1) frames
% $$$      frameName = [vidDirName, num2str(iFrm,'newCR3%.3d'), '.png'];
% $$$      inData(:,:,:,iFrm) = im2single(imread(frameName));
% $$$ end
% $$$ 
% $$$ eval(['mcgmOpticalFlow(inData,filterthreshold,divisionthreshold,divisionthreshold2,speedthreshold,',outVar,');']);
% $$$ 
% $$$ %% Display colour channels
% $$$ % figure('units','pixels','position',[0 0 320 480],'Resize','off');
% $$$ figure(1);
% $$$ load('../Alan_new_model/basis_1stframe_mm.mat');
% $$$ for iO = 0:numOut-1
% $$$     eval(['output=o',num2str(iO),';']);
% $$$     title(['Output ',num2str(iO)]);
% $$$     subplot(1,2,1); imshow(tmpOut(:,:,iO+1),[]); title('MATLAB');
% $$$     for iFrm = 1:1   
% $$$         subplot(1,2,2); imshow(output(:,:,1,iFrm),[]); ...
% $$$             title('HALIDE');
% $$$     end    
% $$$     % saveas(1,['./results/frm_' num2str(iFrm) '.bmp']);
% $$$     pause(1.0);
% $$$ end

%% Display only 1st colour channels
% $$$ figure(1);
% $$$ for iFrm = 1:noFrm        
% $$$     for iO = 0:numOut-1    
% $$$         subplot(2,2,1); imshow(o0(:,:,1,iFrm),[]); title('Chn1');
% $$$         subplot(2,2,2); imshow(o1(:,:,1,iFrm),[]); title('Chn2');
% $$$         subplot(2,2,3); imshow(o2(:,:,1,iFrm),[]); title('Chn3');
% $$$         pause(1.0/24.0);
% $$$     end
% $$$ end

%% --------------------------------------------------------------------------------------------------

% % Continuous input
% for iFrm = 1:noFrm
%      %% Read data from each of first (offset + bufferSize - 1) frames
%      % frameName = [vidDirName, num2str(iFrm,'newCR3%.3d'),
%      % '.png'];
%      frameName = [vidDirName, num2str(iFrm,[testCase '%.3d']), ...
%                   '.png'];
%      %% Use whole frame
%      inData(:,:,:,iFrm) = im2single(imread(frameName));
%      %% Crop central face
%      % tmp = im2single(imread(frameName));
%      % inData(:,:,:,iFrm) = tmp((960/2-frmHeight/2):(960/2+frmHeight/2-1),(1280/2-frmWidth/2):(1280/2+frmWidth/2-1),:);
% end

% Continuous Video Input
% v = VideoReader(inVidName);
% iFrm = 0;
% while hasFrame(v)
%     iFrm = iFrm + 1;
%     inData(:,:,:,iFrm) = im2single(imresize(readFrame(v),[frmHeight frmWidth]));
% end

v = cv.VideoCapture(inVidName);
noFrm = v.FrameCount;
% frmWidth = v.FrameWidth;
% frmHeight = v.FrameHeight;
noColorChan = 3;
inData = single(zeros(frmHeight,frmWidth,noColorChan,noFrm));
for iFrm=1:noFrm
    inData(:,:,:,iFrm) = im2single(imresize(v.read(),[frmHeight frmWidth]));
    % inData(:,:,:,iFrm) = single(imresize(v.read(),[frmHeight frmWidth]))/single(255);    
end

%% Displaying 3-D result
o0 = single(zeros(frmHeight,frmWidth,noFrm));
outVar = 'o0';
numOut = 2;
for iO = 1:numOut-1
    eval(['o',num2str(iO),'=single(zeros(frmHeight,frmWidth,noFrm));']);
    outVar = [outVar,',o',num2str(iO)];
end

% Jump and smooth input by Gaussian
% $$$ for iFrm = 6:4:noFrm
% $$$      %% Read data from each of first (offset + bufferSize - 1) frames
% $$$      % frameName = [vidDirName, num2str(iFrm,'newCR3%.3d'),
% $$$      % '.png'];
% $$$      frameName = [vidDirName, num2str(iFrm,[testCase '%.3d']), '.png'];
% $$$      inData(:,:,:,iFrm) = im2single(imread(frameName));
% $$$ end
% $$$ 
% $$$ % Create bluring effect
% $$$ windowWidth = int16(23);
% $$$ gaussFilter = gausswin(windowWidth,1);
% $$$ gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
% $$$                                               % inData = inData * 6;
% $$$ 
% $$$ for indRow = 1:frmHeight
% $$$     for indCol = 1:frmWidth
% $$$         for indChn = 1:noColorChan
% $$$             inData(indRow,indCol,indChn,:) = ...
% $$$                 conv(squeeze(inData(indRow,indCol,indChn,:)),gaussFilter,'same');
% $$$         end
% $$$     end
% $$$ end

% $$$ gjVidFileHdl = VideoWriter('./Stimulus/generated_jump.avi');
% $$$ gjVidFileHdl.FrameRate = 5;
% $$$ 
% $$$ open(gjVidFileHdl);
% $$$ inData = mat2gray(inData);
% $$$ for iFrm = 1:noFrm
% $$$     writeVideo(gjVidFileHdl,inData(:,:,:,iFrm));
% $$$ end
% $$$ close(gjVidFileHdl);

%% Model parameter for motion model
filterthreshold = single(10^(-5)); % 1e-4.8
divisionthreshold = single(1e-30); % 1e-30
divisionthreshold2 = single(0.99); % 1e-25
speedthreshold = single(1e-6); % 1e-6

eval(['mcgmOpticalFlow_v02_autoscheduled(inData,filterthreshold,divisionthreshold,divisionthreshold2,speedthreshold,',outVar,');']);

%% Model parameter for stereo model
% $$$ filterthreshold = single(1e-5); % 1e-4.8
% $$$ divisionthreshold = single(1e-30); % 1e-30
% $$$ divisionthreshold2 = single(0.99); % 1e-25
% $$$ speedthreshold = single(1e-6); % 1e-6
% $$$ 
% $$$ eval(['mcgmOpticalFlow_v02_stereo(inData,filterthreshold,divisionthreshold,divisionthreshold2,speedthreshold,',outVar,');']);

%% Output with colour conversion

% Alan's wrapping orientation
speed1 = o1; % temporary fix
% $$$ speed1 = speed1 + pi;
speed1 = speed1 + (speed1 >= 2*pi)*(-2*pi);
speed1 = speed1 + (speed1 <= 0)*2*pi;
speed1 = speed1*180/pi;
% figure('units','pixels','position',[0 0 595 479],'Resize','off');
open(vidFileHdl);
for iFrm = 1:noFrm
    T0n = rgb2gray(inData(:,:,:,iFrm)); % temporary fix
    speed0 = o0(:,:,iFrm);
    
    outImg=outputvelocity(T0n,speed0,speed1(:,:,iFrm),16, speedthreshold, ...
                       filterthreshold);
    % imshow(outImg); F = getframe;    
    % writeVideo(vidFileHdl,F.cdata);
    writeVideo(vidFileHdl,uint8(round(im2single(outImg)*255)));
end
close(vidFileHdl);

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
