function display(gm)
% display(gm): Display graphmodel object
%  Prints the variable set, and uses matlab internal disp() to show the table or its size.
%
fprintf('Graphical model over %d variables, containing %d factors\n',gm.nVar,length(gm.factors)); %-gm.nVacant)
if (~isempty(gm.Alg.name))
	fprintf(['  Specialized to algorithm ' gm.Alg.name '\n']);
end;

