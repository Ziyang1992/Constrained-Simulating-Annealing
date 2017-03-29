function [ xbest,objectbest] = csa(xbb1,xlength,xsparsity,supp)
% 2D constrained Simulating Annealing method 
% Details can refer to the article "Phase retrieval for the binary signal£ºtheorem and algorithm"
% By Ziyang Yuan and Hongxia Wang 
G=@(x)fft2(x);
 Gtrans=@(x)ifft2(x)*length(x);
 y=abs(G(xbb1)).*abs(G(xbb1));
 loss=@(x)(norm(abs(G(x)).*abs(G(x))-y)/norm(y));
countrestart=0;
MaxTime=20;
%debug
tic;
r = 0.99;
T = xlength/2;

k=xsparsity;
xlen=xlength;
setArea=supp*supp;
S=randperm(setArea);
%initialize S
S=S(1:k);
xiter1=zeros(xlen,1);
%initialize x_iter
xiter1(S)=1; 
xbest=xiter1;
xiter=reshape(xiter1,supp,supp);
objectbest=loss(xiter);
objectlast=loss(xiter);
while 1    
    T=r*T;
    xiter=reshape(xiter1,supp,supp);
    Gx=G(xiter);
    grad1 = Gtrans(Gx.*(abs(Gx).*abs(Gx)-y)); 
    I1=grad1.*xiter.*(grad1>0); 
    I1=reshape(I1,supp*supp,1);
    table=cumsum(I1);
    table=table/table(end);  
    i1=(find(rand<table,1,'first')); 
    
    xMidIter=xiter1;
    xMidIter(i1)=0;
    xMidIter1=reshape(xMidIter,supp,supp);
    Gx=G(xMidIter1);
    grad2=Gtrans(Gx.*(abs(Gx).*abs(Gx)-y));  
    I2=-grad2.*(1-xMidIter1); 
    I2=reshape(I2,supp*supp,1);
    [~,index2]=max(I2);
    
    xnew=xMidIter;
    xnew(index2)=1;
    xnew1=reshape(xnew,supp,supp);
    lossnew=loss(xnew1);
    if lossnew<objectlast
       xiter1=xnew;
       objectlast=lossnew;
       if objectlast<objectbest
           xbest=xiter1;
           objectbest=objectlast;
       end
    elseif rand(1)<exp((objectlast-lossnew/T))
        xiter1=xnew;
        objectlast=lossnew;
    end
    if objectbest<1e-2
        break; 
    end
    if T<20
        T=xlength/2;
        S=randperm(setArea);
        S=S(1:k); 
        xiter1=zeros(xlen,1);
        xiter1(S)=1; 
        countrestart=countrestart+1;
        objectlast=Inf;
        if toc>MaxTime
            break
        end
    end
end
end

