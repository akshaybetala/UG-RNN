function writeWCSP(placeholder,filename,flist)
% writeWSCP : write out a collection of factors to a file in weighted CSP format
% writeWCSP(factor(),filename,flist)
% Write 'filename' in weighted CSP format (see e.g. http://graphmod.ics.uci.edu/group/WCSP_file_format)
% TODO !!!: sparsity / use most common value in output

% (c) Alexander Ihler 2012

fp=fopen(filename,'w');
if (fp==-1) error(sprintf('writeWCSP: Failed to open file %s for write',filename)); end;

nvar=0; dmax=0; ub=-inf;
for i=1:length(flist), 
  nvar = max(nvar, max(vars(flist{i}))); 
  dmax = max(dmax, max(dims(flist{i})));
  ub   = max(ub,   table(maxmarginal(flist{i},[])) );
end;

dim=zeros(1,nvar);
for i=1:length(flist), dim(vars(flist{i}))=dims(flist{i}); end;

% Preamble: problem name, # variables, max # dimensions, # constraints, global value upper bound
fprintf(fp,'%s %d %d %d %d\n',filename,nvar,dmax,length(flist),ub);

% Dimensions of each variable
fprintf(fp,' %d',dim); fprintf(fp,'\n');

% Print each factor as a list of weights
for i=1:length(flist), 
	v=vars(flist{i}); 
  fprintf(fp,'%d ',length(v));            % # variables in scope
	fprintf(fp,' %d',v-1);                  % variable IDs in order
  fprintf(fp,'  0  %d\n',numel(flist{i}));  % default value 0, all tuples different (!!!)

  tab = table(flist{i});                  % now all tuples
  for t=1:numel(flist{i}),
    sub = ind2subv(flist{i},t);           % their state configuration
    fprintf(fp,'%d ',sub-1); 
    fprintf(fp,'%d\n',tab(t));            % and value
  end;
end;

fclose(fp);


