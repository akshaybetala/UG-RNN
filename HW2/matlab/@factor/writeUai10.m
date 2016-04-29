function writeUai10(placeholder,filename,flist)
% writeUai10 : write out a collection of factors to a file in UAI-2010 format
% writeUAI10(factor(),filename,flist)
% Write 'filename' in UAI-2010 competition format (see e.g. http://www.cs.huji.ac.il/project/UAI10/fileFormat.php)
% FIX: writeUAI*: no evidence support

% (c) Alexander Ihler 2010

fp=fopen(filename,'w');
if (fp==-1) error(sprintf('writeUAI10: Failed to open file %s for write',filename)); end;

% Factor lists are Markov networks
fprintf(fp,'MARKOV\n');

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
  fprintf(fp,' %f',t(:));
	fprintf(fp,'\n');
end;

fclose(fp);


