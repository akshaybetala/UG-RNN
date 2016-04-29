function [flist,name]=readWCSP(placeholder,filename)
% readWCSP : read in a collection of factors in WCSP (weighted constraint satisfaction) format
% [flist] = readWCSP(factor(),filename)
% Read 'filename' in WCSP format (see e.g. http://graphmod.ics.uci.edu/group/WCSP_file_format)

% (c) Alexander Ihler 2010

fp=fopen(filename);
if (~fp) error(sprintf('readWCSP: File %s does not exist',filename)); end;

str=getNextRegexp(fp,'^\S*');    % get first line
spaces=strfind(str,' ');
name=str(1:spaces(1));
prologue=sscanf(str(spaces(1)+1:end),'%d %d %d %ld');
nVariables=prologue(1); maxdim=prologue(2); nConstraints=prologue(3); upperbound=prologue(4);

str=getNextRegexp(fp,'^[0-9]');    % get next numeric line
dims=sscanf(str,'%d',nVariables); dims=dims'; 

% Now read the constraints
for c=1:nConstraints,
  str=getNextRegexp(fp,'^\s*[0-9]');    % get constraint info
  tmp=sscanf(str,'%d'); tmp=tmp';      % involved variables, base cost, and # of deviations
  nvar=tmp(1); vars=tmp(1+(1:nvar))+1; defaultCost=tmp(nvar+2); nTuples=tmp(nvar+3);
  varsc=num2cell(dims(vars));
  table=zeros(varsc{:})+defaultCost;              % currently store as a big table
  for i=1:nTuples,
    str=getNextRegexp(fp,'^\s*[0-9]');    % 
    tuple=sscanf(str,'%d',nvar+1);                % read tuple and
    val=tuple(nvar+1); tuple=tuple(1:nvar)';      %  its associated weight
    tupleC=num2cell(tuple+1); 
    table(tupleC{:})=val; 
  end;
  flist{c}=factor(vars,table);
end;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str=getNextRegexp(fp,pattern)
while ( fp )
  str=fgets(fp); 
  if (length(regexp(str,pattern)) > 0) break; end;
end;
if (~fp) error('readWCSP: Premature end of file encountered'); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



