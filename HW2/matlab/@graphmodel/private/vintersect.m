function s = vintersect(s1,s2)
% intersection of two sorted uint32s

r=vmember(s1,s2);
s=s1(r);

