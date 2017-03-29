% Implemetation of the Constrained Simulating Annealing method presented
% in the paper "Phase retrieval for the binary signal£ºuniqueness and algorithm"
% By Ziyang Yuan and Hongxia Wang
% This simulating experiments are aim for 2D David star
% Code written by Ziyang Yuan
% Date:2017.3.22
clc;
clear;
load 'David'
Pic=I;
[m,n]=size(I);
for i=1:m
    for j=1:n
        if(Pic(i,j)==1)
            Pic(i,j)=0;
        else
            Pic(i,j)=1;
        end
    end
end
imshow(Pic)
c1=abs(fft2(Pic)).^2; 
geshu=sqrt(round(c1(1,1)));
[x,objectbest]=csa(Pic,m*m,geshu,m);
xBest=bestMatch2D(x,Pic);
imagesc(xBest);
title('Recoverd');
figure;
imagesc(Pic);
title('Real');