 function [xbest,objectbest]=csa(xbb,x2,k,n)
 % 1D constrained Simulating Annealing method 
% Details can refer to the article "Phase retrieval for the binary signal��theorem and algorithm"
% By Ziyang Yuan and Hongxia Wang 
 tic;
 r=0.95;
 L=1;
 T=n;
 xbb1=[xbb;zeros(n,1)];
 G=@(x)fft(x);
 Gtrans=@(x)ifft(x)*length(x);
 y=abs(G(xbb1)).*abs(G(xbb1));
 loss=@(x)(norm(abs(G(x)).*abs(G(x))-y));
 maxtime=60;
 if(k==1)
     xbest=zeros(n,1);
     xbest(1)=1;
     objectbest=loss([xbest;zeros(n,1)]); 
     return;
 end
 if(k==2)
     xbest=zeros(n,1);
     xbest(1)=1;
     mm=find(x2);
     xbest(1+mm)=1;
     objectbest=loss([xbest;zeros(n,1)]); 
     return;
 end
 %% construct the feasible set
  ge=k;
  x3=zeros(n,1);
  [yin1,~]=find(x2);
  d1=x2(yin1);
  yin=[d1';yin1'];
  sou=zeros(1,ge+1);
  sou(1)=1;
  x3(1)=1;
  i=1;
  [~,mn]=size(yin);
  %%%%%%%%%%%%%
  x3(1+yin(2,mn))=1;
  sou(2)=1+yin(2,mn);
  yin(:,mn)=[];
  h=1;
  xiabiao=search(yin,sou,mn-1); %% get the feasible set
  xiabiao1=xiabiao;
  he=sou(1)+sou(2);
  [~,si1]=size(xiabiao);
 while(i<=si1)
      ind= find(xiabiao==he-xiabiao(i));
      if(isempty(ind)==0)
      fenzu(1,h)=xiabiao(i);
      fenzu(2,h)=xiabiao(ind);
      h=h+1;
      xiabiao([i,ind])=[];
      i=1;
      [~,si1]=size(xiabiao);
      else
          i=i+1;
      end
 end
 [~,si]=size(fenzu);
 xchu=x3;
 [x3,loc1]=initialization(x3,si,k,fenzu); %% estimate evaluation
 xbest=x3;
 objectbest=loss([x3;zeros(n,1)]);
 objectprevious=objectbest;
 %% Constrained Simulating Annealing method
 while(1)
 T=r*T;
 x4=[x3;zeros(n,1)];
 grad1=Gtrans(G(x4).*(abs(G(x4)).*abs(G(x4))-y));
 I1=grad1.*x4.*(grad1>0);
[~,i1]=max(I1);
 xMidIter=x4;
xMidIter(i1)=0;
  grad2=-Gtrans(G(xMidIter).*(abs(G(xMidIter)).*abs(G(xMidIter))-y));
  xiabiao2=xiabiao1;
  [~,~,k2]=intersect(loc1,xiabiao2);
  xiabiao2(k2)=[];
 guodu=zeros(2*n,1);
 guodu(xiabiao2)=1;
I2=grad2.*guodu;
[~,i2]=max(I2);
 xMidIter(i2)=1;
 objectlast=loss(xMidIter); 
 if(objectlast<objectprevious)
     x3=xMidIter(1:n);
     xbest=x3;
     objectbest=objectlast;
     objectprevious=objectlast;
 elseif rand(1)<exp((objectprevious-objectlast/T))
        x3=xMidIter(1:n);
        objectprevious=objectlast;
 end
 if(T<20)
    L=L+1;
    T=n;
    [x3,~]=initialization(xchu,si,k,fenzu);
    objectprevious=loss([x3;zeros(n,1)]);
 end
   if(objectbest<=0.001)
       return
   end
   if(toc>maxtime)
     return
   end
   loc1=find(x3);
 end
 end


