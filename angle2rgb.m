function [ ov] = angle2rgb( v )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[rows, cols, depth]=size(v);

for i=1:rows
    for j=1:cols       
        a=v(i,j,1);
        if a>=0 && a<pi/2
            ov(i,j,1)=1;
            ov(i,j,2)=sin(a);
            ov(i,j,3)=0;
        elseif a>=pi/2 && a<pi
            ov(i,j,1)=cos(a-pi/2);
            ov(i,j,2)=1;
            ov(i,j,3)=0;
        elseif a>=pi && a<3*pi/2
            ov(i,j,1)=0;
            ov(i,j,2)=cos(a-pi);
            ov(i,j,3)=sin(a-pi);
        elseif a>=3*pi/2 && a<=2*pi
            ov(i,j,1)=sin(a-(3*pi/2));
            ov(i,j,2)=0;
            ov(i,j,3)=cos(a-(3*pi/2));      
    end 
end
   
end

