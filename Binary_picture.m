% Implemetation of the Constrained Simulating Annealing method presented
% in the paper "Phase retrieval for the binary signal£ºuniqueness and algorithm"
% By Ziyang Yuan and Hongxia Wang
% This simulating experiments are aim for 2D sparse image
% Code written by Ziyang Yuan
% Date:2017.1.16
% 
clc;
clear;
number=30; %the number of crosses
replication=50; % the replication time
record=zeros(number,1);
for k=1:number
    tt=0;
    for times=1:replication
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
quantity=sqrt(round(c1(1,1)));
 [x6,objectbest]=csa(Pic,64*64,quantity,64);
 if objectbest<10^-3
     tt=tt+1;
 end
    end
    record(k)=tt/times;
end
h=plot(1:number,record,'*');
 xlabel('number of crosses'), ylabel('Mean recovery rate'), ...
     title('2D constrained SA')
 str='Obj';
legend(h,str);