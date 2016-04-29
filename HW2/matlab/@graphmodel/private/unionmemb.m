function [u,m1,m2]=unionmemb(a,b)
% union & membership info for sorted uint32s
 u=union(a,b);
 m1=ismember(u,a);
 m2=ismember(u,b);
