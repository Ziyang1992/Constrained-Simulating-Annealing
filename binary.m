% Implemetation of the Constrained Simulating Annealing method presented
% in the paper "Phase retrieval for the binary signal£ºtheorem and algorithm"
% By Ziyang Yuan and Hongxia Wang
% This simulating experiments are aim for 1D sparse signal
% Code written by Ziyang Yuan
% Date:2017.1.16
% modified: 2017.3.27
clc;
clear;
%% set the initial parameters
tic;
n=1000;
supp=n/2;
value=30;   % the range of the sparsity k
replication=50; % the replicate times
record1=zeros(1,value);
record2=zeros(1,value);
for k=1:value   % The range of sparsity k 
    tt=0;
    tt1=0;
    for times=1:replication  % The test times of each sparsity k
locs=randperm(n/2);
x=zeros(n/2,1);
x(locs(1:k))=1;
x1=[x;zeros(n/2,1)];
c1=abs(fft(x1)).^2;
G=@(x)fft(x);
Ginv=@(x)ifft(x);
Gtrans=@(x)ifft(x)*length(x);
Gm=fft(eye(n));
loss=@(x)(norm(abs(G(x)).*abs(G(x))-c1)/norm(c1));
%% start the binary phase retrieval 
xco=ifft(c1);
geshu=round(xco(1));
x2=xco(2:n/2);
x2(x2<0.00001)=0;
x2=round(x2);
[x3,objectbest]=csa(x,x2,geshu,n/2);
 if objectbest<10^-3   % one criterion: calculate the loss function
     tt=tt+1;
 end
xBest=bestMatch(x3,x); % another stopping criterion: find whether it is equal to orignal signal
  if norm(xBest-x)<10^-3
     tt1=tt1+1;
   end
    end
   record1(k)=tt/times;
   record2(k)=tt1/times;
end
h=plot(1:value,record1,'o',1:value,record2,'*');
 xlabel('Sparsity'), ylabel('Mean recovery rate'), ...
     title('1D constrained SA')
 str=['Obj';'MSE'];
legend(h,str);