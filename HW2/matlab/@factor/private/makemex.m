% Script to compile mex (c++) versions of various functions

if (exist('vdiff.cpp'))
	  mex vdiff.cpp 
end;
if (exist('vmember.cpp'))
	  mex vmember.cpp 
end;
if (exist('vunion.cpp'))
	  mex vunion.cpp 
end;
if (exist('vunionmemb.cpp'))
	  mex vunionmemb.cpp 
end;
if (exist('vintersect.cpp'))
	  mex vintersect.cpp 
end;
%if (exist('repmat.c')&&exist('mexutil.c'))
%	  mex repmat.c mexutil.c
%end;

