function u = vunion(a,b)
% union of two sorted uint32s

 %u = unique([a b]);

 if (isempty(a)) u=b; return; end;
 if (isempty(b)) u=a; return; end;

 r=false(1,max(max(a),max(b))+1);
 r(a+1)=true;
 r(b+1)=true;
 u=uint32(find(r))-1;

