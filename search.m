% Through the support information to find the feasible set
function xiabiao=chazhao(yin,sou,k)
i=1;
m=2;
kk=1;
while(i<=k)
        a1=abs(1+yin(2,i)-sou(1:m));
        pan1=prod(double(ismember(a1,yin(2,:))));
        if(pan1==1)
            xiabiao(1,kk)=1+yin(2,i);
            kk=kk+1;
        end
       i=i+1;
end
