%%% find the initialization through the feasible set
function [x3,loc1]=initialization(x3,si,k,fenzu)
 if(si>=k-2)
      S1=randperm(si);
      S1=S1(1:k-2);
      for i=1:k-2
          kk=randperm(2);
          x3(fenzu(kk(1),S1(i)))=1;
      end
      loc1=find(x3);
  else
      S1=randperm(si);
      S2=S1(1:k-2-si);
      for i=1:si
          if(any(S2==i)==1)
              x3(fenzu(1,i))=1;
              x3(fenzu(2,i))=1;
          else
             kk=randperm(2);
            x3(fenzu(kk(1),S1(i)))=1;
          end
      end
      loc1=find(x3);
  end