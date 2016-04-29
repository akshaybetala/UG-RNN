function m = vmember(a,b)
% membership binary vector for sorted uint32s

 %m=ismember(a,b);

 if (isempty(a)||isempty(b)) m=false(size(a)); return; end;

 r=false(1,max(max(a),max(b))+1);
 r(b+1)=true;
 m=r(a+1);

