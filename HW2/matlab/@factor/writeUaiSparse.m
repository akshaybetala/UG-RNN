function writeUai10(placeholder,filename,flist)
% writeUaiSparse : write out a collection of factors to a file in UAI-SPARSE format
% writeUaiSparse(factor(),filename,flist)
% Write 'filename' in UAI-SPARSE competition format (see TODO)
% TODO FIX: writeUAI*: no evidence support

% (c) Alexander Ihler 2010, 2015

fp=fopen(filename,'w');
if (fp==-1) error(sprintf('writeUaiSparse: Failed to open file %s for write',filename)); end;

% Factor lists are Markov networks
fprintf(fp,'SPARSE\n');

% Print # of variables and dimensions
nv=0;
for i=1:length(flist), nv=max(nv,max(vars(flist{i}))); end;
fprintf(fp,'%d\n',nv);
dim=zeros(1,nv);
for i=1:length(flist), dim(vars(flist{i}))=dims(flist{i}); end;
fprintf(fp,' %d',dim); fprintf(fp,'\n');

fprintf(fp,'\n');

% Print # of factors and each factor's variables
fprintf(fp,'%d\n',length(flist));
for i=1:length(flist), 
	v=vars(flist{i}); fprintf(fp,'%d ',length(v)); 
	fprintf(fp,' %d',v(end:-1:1)-1); 
	fprintf(fp,'\n'); 
end;

fprintf(fp,'\n');

for i=1:length(flist),
	t=table(flist{i}); fprintf(fp,'%d ',numel(t));
  tstart  = [1;diff(t(:))]~=0; 
  tval    = t(tstart);
  tlength = diff([find(tstart);length(t(:))+1])';
  for j=1:length(tval),
    if (tlength(j)==1) fprintf(fp,' %f',tval(j));
    else               fprintf(fp,' (%d:%f)',tlength(j),tval(j));
    end;
  end;
  %fprintf(fp,' %f',t(:));
	fprintf(fp,'\n');
end;

fclose(fp);


