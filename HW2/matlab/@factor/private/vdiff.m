function d = vdiff(a,b)
% set difference of two sorted uint32s

 %d=setdiff(a,b);

 if (isempty(a)) d=a; return; end;
 if (isempty(b)) d=a; return; end;

 r=false(1,max(max(a),max(b))+1);
 r(a+1)=true;
 r(b+1)=false;
 d=uint32(find(r))-1;
