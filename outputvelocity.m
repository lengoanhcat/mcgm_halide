function [ I ] = outputvelocity(Blur, Speed, Angle, border, speedthreshold, filterthreshold)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Scale the grey level images
Blur=(Blur-min(Blur(:)))./(max(Blur(:))-min(Blur(:)));
Blur=cat(3, Blur,Blur,Blur);

%Speed scaled to 1
Speed=1*Speed;
Speed=cat(3, Speed,Speed,Speed);
% Speed(:,:,:) = Speed(:,:,:) .* (abs(Speed(:,:,:)) > speedthreshold);
% Speed(:,:,:) = Speed(:,:,:) .* (abs(Speed(:,:,:)) > filterthreshold);

%Use the log speed to visualise speed
LogSpeed=log10(Speed+0.0000001);
LogSpeed=(LogSpeed-min(LogSpeed(:)))./(max(LogSpeed(:))-min(LogSpeed(:)));



%Make a colour image
[rows, cols, depth] = size(Angle);

%Do it the HSV way
c=zeros(rows,cols,3);
Angle=Angle/360;

c(:,:,1)=Angle;
c(:,:,2)=1;
c(:,:,3)=1;
% c(:,:,:) = c(:,:,:) .* (abs(Speed(:,:,:)) > speedthreshold);
% c(:,:,:) = c(:,:,:) .* (abs(Speed(:,:,:)) > filterthreshold);

%Do hsv to rgb
c1=hsv2rgb(c);

%Make the border
bir=rows+2*border;
bic=cols+2*border;
orows=bir/2;
ocols=bir/2;

%rotation matrix
ph=0;
rot=[cos(ph), -sin(ph);sin(ph), cos(ph)];

cb=zeros(bir,bic,3);
sb=cb;
mb=zeros(bir,bic,1);

for i=1:bir
    for j=1:bic
        if (i<border || i>=rows+border || j<border || j>=cols+border)
            co(1,1)=(i-orows); co(2,1)=-(j-ocols);
            rco=(rot*co);
            mb(i,j)=(atan2(rco(1,1),rco(2,1))+pi);
            a=mb(i,j)/(2*pi);
            cb(i,j,1)=a;
            cb(i,j,2)=1;
            cb(i,j,3)=1;
            
            sb(i,j,:)=1;
        end
    end
end

% %COMMENT IN HERE FOR HSV
% %make the border rgb
% cb=hsv2rgb(cb);
% 
% %get the data for the brightness modified map
% c2=c;
% c2(:,:,3)=LogSpeed(:,:,1);
% c2=hsv2rgb(c2);
% %COMMENT IN HERE FOR HSV

%COMMENT IN HERE FOR OLD SCHOOL
%make the old (school) border
cb=angle2rgb(mb);

%get the old data

c1=angle2rgb(c(:,:,1)*2*pi);
c1(:,:,:) = c1(:,:,:) .* (abs(Speed(:,:,:)) > speedthreshold);
% c1(:,:,:) = c1(:,:,:) .* (abs(Speed(:,:,:)) > filterthreshold);
c2=c1(:,:,:).*Speed(:,:,:);
%COMMENT IN HERE FOR OLD SCHOOL

%put the data in the border

cb(border:rows+border-1, border:cols+border-1,:)=c1(:,:,:); ang1=cb;
cb(border:rows+border-1, border:cols+border-1,:)=c2(:,:,:); ang2=cb;
sb(border:rows+border-1, border:cols+border-1,:)=Speed(:,:,:); Speed=sb;
sb(border:rows+border-1, border:cols+border-1,:)=Blur(:,:,:); Blur=sb;


I = cat(2,cat(1,Blur,Speed),cat(1,ang1,ang2));
% imtool(I, [min(I(:)) max(I(:))] );

end
