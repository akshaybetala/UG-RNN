function display(f)
% Display factor object information
% display(f): prints the variable set, and uses matlab internal disp() to show the table or its size.
%

% (c) Alexander Ihler 2010

%if (~isempty(f))
fprintf('Variables: '); fprintf('%d ',double(f.v));
fprintf('\nTable: \n' ); disp(f.t);
%else fprintf('Empty factor\n');
%end;

