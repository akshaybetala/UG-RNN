function [flist,names,labels,evid]=readErgo(placeholder,filename)
% readErgo : read in a collection of factors from a file in ERGO format
% [flist,names,labels] = readErgo(factor(),filename)
% Read 'filename' in Ergo format (see e.g. http://graphmod.ics.uci.edu/group/Ergo_file_format)

% (c) Alexander Ihler 2010

fp=fopen(filename,'r');
if (fp==-1) error(sprintf('readErgo: File %s does not exist',filename)); end;

str=getNextRegexp(fp,'^[0-9]');    % get next numeric line
nvar=sscanf(str,'%d',1);
str=getNextRegexp(fp,'^[0-9]');    % get next numeric line
dims=sscanf(str,'%d'); dims=dims'; 
if (length(dims)~=nvar) error(sprintf('readErgo: dimension list does not match number of variables in %s',filename)); end;
dimCell=cell(1,length(dims)); for i=1:length(dims), dimCell{i}=dims(i); end;
npar=zeros(1,nvar); par=cell(1,nvar);
for i=1:nvar,
  str=getNextRegexp(fp,'^[0-9]');    % get next numeric line
  tmp=sscanf(str,'%d'); tmp=tmp';
  npar(i)=tmp(1); par{i}=tmp(2:end)+1;
  if (length(par{i})~=npar(i)) error(sprintf('readErgo: wrong number of parents for var %d in %s',i,filename)); end;
end;

str=getNextRegexp(fp,'/\*\s*Probabilities\s*\*/');    % advance to next region
flist=cell(1,nvar);
for i=1:nvar,
  %fprintf('.');
  str=getNextRegexp(fp,'^[0-9]');    % get next numeric line
  nvals=sscanf(str,'%d',1);
  par{i} = par{i}( dims(par{i})~=1 );
  if (prod(dims([i,par{i}]))~=nvals) error(sprintf('readErgo: wrong number of values for factor %d in %s',i,filename)); end;
  if (isempty(par{i})) table=zeros(1,dimCell{i});
  else                 table=zeros(dimCell{[i,par{i}(end:-1:1)]});
  end;
  v=1; while(v<=nvals),
    str=getNextRegexp(fp,'^\s*[0-9]');    % get next numeric line (inital space ok)
    vals=sscanf(str,'%f',dims(i));
    table(v:v+dims(i)-1)=vals; v=v+dims(i);
  end;
  flist{i}=factor([i,par{i}(end:-1:1)],table);
end;
fprintf('Read %d conditional probabilities.\n',nvar);

str=getNextRegexp(fp,'/\*\s*Names\s*\*/');    % advance to next region (???)
%str=getNextRegexpDebug(fp,'/\**\sNames*\s\*/');    % advance to next region
names=cell(1,nvar);
for i=1:nvar,
  %fprintf('n');
  str=getNextRegexp(fp,'\w*');        % get next word-line
  tmp=regexp(str,'\w*','match');
  names{i}=tmp{1};  %  and get variable name
end;
%fprintf('Read %d variable names.\n',nvar);

str=getNextRegexp(fp,'/\*\s*Labels\s*\*/');    % advance to next region
labels=cell(1,nvar);
for i=1:nvar,
  %fprintf('l');
  str=getNextRegexp(fp,'\w*');        % get next word-lines
  labels{i}=regexp(str,'\w*','match');      %  and parse them into variable value labels
end;
%fprintf('Read %d value labels.\n',nvar);

fclose(fp);

evid={};
fp = fopen([filename '.evid']);
if (fp~=-1)      % if there's an evidence file
  fprintf('Reading evidence file %s\n',[filename '.evid']);
  str=getNextRegexp(fp,'/\*\s*Evidence\s*\*/');    % load preamble
  str=getNextRegexp(fp,'^[0-9]');      % get number of evidence values
  nevid=sscanf(str,'%d',1);
  evid=cell(1,nevid);
  for i=1:nevid,
    str=getNextRegexp(fp,'^\s*[0-9]');    % get (variable,value) pair
  tmp=sscanf(str,'%d',2); tmp=tmp'+1;
  table=zeros(1,dims(tmp(1))); table(tmp(2))=1;
  evid{i}=factor(tmp(1),table);      % convert to a factor
  end;
  fclose(fp);
end;
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str=getNextRegexp(fp,pattern)
str='';
while ( fp && length(regexp(str,pattern))==0)
  str=fgets(fp); 
end;
if (~fp) error('readErgo: Premature end of file encountered'); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



