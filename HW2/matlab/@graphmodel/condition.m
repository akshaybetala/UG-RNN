function gm = condition(gm, vars, val)
% gm = condition(GM,vars,tuple) : create a "slice" graphical model with "vars" fixed to value "tuple"
%   Conditions each factor in the graphical model, returning a new graphical model over the other variables

% (c) Alexander Ihler 2013

for f=withVariable(gm,vars)
  gm.factors{f} = condition(gm.factors{f},vars,val);
end;

%fs = getFactor(GM);
%for f=1:length(fs), if (~isempty(fs{f})),
%  fs{f} = condition(fs{f}, vars,val);
%end; end;
%gm = graphmodel(fs);

