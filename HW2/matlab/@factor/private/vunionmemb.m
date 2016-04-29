function [u,m1,m2]=vunionmemb(a,b)
% union & membership info for sorted uint32s

 %u=vunion(a,b);
 %m1=ismember(u,a);
 %m2=ismember(u,b);
 %return;

 if (isempty(a)) u=b; return; end;
 if (isempty(b)) u=a; return; end;

 r=false(1,max(max(a),max(b))+1);
 a=a+1; b=b+1;
 r(a)=true;
 r(b)=true;
 u=uint32(find(r));
 r(b)=false;
 r(a)=true;
 m1 = r(u);
 r(a)=false;
 r(b)=true;
 m2 = r(u);
 u = u-1;

