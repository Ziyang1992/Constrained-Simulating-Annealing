% Implemetation of the Constrained Simulating Annealing method presented
% in the paper "Phase retrieval for the binary signal£ºtheorem and algorithm"
% By Ziyang Yuan and Hongxia Wang
% This simulating experiments are aim for 2D sparse image
% Code written by Ziyang Yuan
% Date:2017.1.16
clc;
clear;
sparsity=30;
jilu=zeros(sparsity,1);
for k=1:sparsity
    tt=0;
    for cishu=1:50
Pic=zeros(64,64);
%% randomly create the image
S=randperm(21*21);
locs=S(1:k);
for i=0:20
    for j=0:20
        xiabiao=j*21+i+1;
        ind= find(locs==xiabiao);
      if(isempty(ind)==0)
for x=1:64
    for y=1:64
        if (x-3-3*i)^2+(y-3-3*j)^2<=1^2                              
            Pic(x,y)=1;
        end    
    end
end
      end
    end
end
c1=abs(fft2(Pic)).^2; 
c2=ifft2(c1);
geshu=sqrt(round(c1(1,1)));
 [x6,objectbest]=csa(Pic,64*64,c2,64);
 if objectbest<10^-3
     tt=tt+1;
 end
    end
    jilu(k)=tt/cishu;
end