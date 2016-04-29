function [flist,names,labels,evid]=readUai10(placeholder,filename)
% readUai10 : read in a collection of factors from a file in UAI-2010 format
% [flist,names,labels,evidence] = readUAI10(factor(),filename)
% Read 'filename' in UAI-2010 competition format (see e.g. http://www.cs.huji.ac.il/project/UAI10/fileFormat.php)
% FIX: readUAI*: evidence read incorrectly; hybrid between UAI08 and UAI10

% (c) Alexander Ihler 2010

fp=fopen(filename);
if (fp==-1) error(sprintf('readUAI10: File %s does not exist',filename)); end;

str=getNextRegexp(fp,'^\S*');    % get type string
%str=str(1:end-1);  % strip newline
str=strtrim(str);
switch(lower(str))
% case 'bayes',
%  error('Bayes not yet supported?');
% case 'markov',
  case {'bayes', 'markov'},
  str=getNextRegexp(fp,'^\s*[0-9]');    % get number of variables
  nvars=sscanf(str,'%d',1);
  str=getNextRegexp(fp,'^\s*[0-9]');    % get dimensions
  dims=sscanf(str,'%d'); dims=dims'; 
  if (length(dims)~=nvars) error(sprintf('readUAI10: dimension list does not match number of variables in %s',filename)); end;
  dimCell=cell(1,length(dims)); for i=1:length(dims), dimCell{i}=dims(i); end;

  str=getNextRegexp(fp,'^[0-9]');    % get number of cliques
  ncliques=sscanf(str,'%d',1);
  perClique=zeros(1,ncliques); varClique=cell(1,ncliques); nclique=zeros(1,ncliques);
  for i=1:ncliques,
    str=getNextRegexp(fp,'^\s*[0-9]');    % get clique size and dependencies
    tmp=sscanf(str,'%d'); tmp=tmp';
    nclique(i)=tmp(1); varClique{i}=tmp(2:end)+1;
    if (length(varClique{i})~=nclique(i)) error(sprintf('readUAI10: wrong number of variables for clique %d in %s',i,filename)); end;
  end;
end;

%str=getNextRegexp(fp,'/\*\s*Probabilities\s*\*/');    % advance to next region
flist=cell(1,ncliques);
for i=1:ncliques,
  fprintf('.');
  str=getNextRegexp(fp,'^[0-9]');    % get number of values listed
  [nvals,nread,emsg,pos]=sscanf(str,'%d',1);                                 % Octave has more limited sscanf:
  if (pos == 0 && nread ~= 0) [st,en]=regexp(str,'^[0-9]'); pos=en+1; end;   % fix for Octave
  varClique{i} = varClique{i}( dims(varClique{i})~=1 );
  if (prod(dims(unique(varClique{i})))~=nvals), error(sprintf('readUAI10: wrong number of values for factor %d in %s',i,filename)); end;
  %fprintf('Read clique %d, %d values\n',i,nvals); 
	table=zeros(dimCell{[varClique{i}(end:-1:1)]},1);
	v=1;
	[vals, nread]=sscanf(str(pos:end),'%f',nvals);
  table(v:v+nread-1)=vals; v=v+nread;
  while(v<=nvals),
    str2=getNextRegexp(fp,'^\s*[0-9]');    % get next numeric line (inital space ok)
    [vals, nread]=sscanf(str2,'%f',nvals);
    table(v:v+nread-1)=vals; v=v+nread;
  end;
	flist{i}=factor([varClique{i}(end:-1:1)],table);    % Matlab is little-endian; UAI is big-endian
end;
fprintf('\n');

% %%% No name region (?) %%%
%str=getNextRegexp(fp,'/\*\s*Names\s*\*/');    % advance to next region (???)
%names=cell(1,nvars);
%for i=1:nvars,
%  fprintf('n');
%  str=getNextRegexp(fp,'\w*');        % get next word-line
%  tmp=regexp(str,'\w*','match');
%  names{i}=tmp{1};  %  and get variable name
%end;
%fprintf('\n');
names={};

% %%% No value list region (?) %%%
%str=getNextRegexp(fp,'/\*\s*Labels\s*\*/');    % advance to next region
%labels=cell(1,nvars);
%for i=1:nvars,
%  fprintf('l');
%  str=getNextRegexp(fp,'\w*');        % get next word-lines
%  labels{i}=regexp(str,'\w*','match');      %  and parse them into variable value labels
%end;
%fprintf('\n');
labels={};

fclose(fp);

evid={};
if (length(dir([filename '.evid'])))      % if there's an evidence file
  fprintf('Reading evidence file %s\n',[filename '.evid']);
  fp = fopen([filename '.evid']);
  %str=getNextRegexp(fp,'/\*\s*Evidence\s*\*/');    % load preamble
%  str=getNextRegexp(fp,'^[0-9]');      % get number of evidence values
%  nevid=sscanf(str,'%d',1);
  nevid = 1;
  evid=cell(1,nevid);
  for i=1:nevid,
    str=getNextRegexp(fp,'^\s*[0-9]');    % get (variable,value) pair
    tmp=sscanf(str,'%d'); tmp=tmp'+1;
		if (tmp(1)*2-2 ~= length(tmp)-1) fprintf('Warning: incorrect number of evidence values?\n'); tmp, end;
		tmp=tmp(2:end);
    evid{i}=cell(1,length(tmp)/2);      % # of var/value pairs in this sample
    for j=0:(length(tmp)/2-1),
      table=zeros(1,dims(tmp(2*j+1))); table(tmp(2*j+2))=1;
      evid{i}{j+1}=factor(tmp(2*j+1),table);    % convert to a factor
    end;
  end;
  fclose(fp);
end;
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function str=getNextRegexp(fp,pattern)
%str=' ';
while ( fp )
  str=fgets(fp); 
  if (length(regexp(str,pattern)) > 0) break; end;
end;
if (~fp) error('readUai: Premature end of file encountered'); end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



