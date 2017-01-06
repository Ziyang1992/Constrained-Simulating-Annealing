% Implemetation of the Constrained Simulating Annealing method presented
% in the paper "Phase retrieval for the binary signal£ºtheorem and algorithm"
% By Ziyang Yuan and Hongxia Wang
% This simulating experiments are aim for 1D sparse signal
% Code written by Ziyang Yuan
% Date:2017.1.16
clc;
clear;
%% set the initial parameters
tic;
n=1000;
supp=n/2;
sparsity=30;
jilu=zeros(sparsity,1);
for k=1:sparsity
    tt=0;
    for cishu=1:50
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
geshu=xco(1);
x2=xco(2:n/2);
x2(x2<0.00001)=0;
x2=round(x2);
[x3,objectbest]=csa(x,x2,k,n/2);
 if objectbest<10^-3
     tt=tt+1;
 end
    end
   jilu(k)=tt/cishu;
end
