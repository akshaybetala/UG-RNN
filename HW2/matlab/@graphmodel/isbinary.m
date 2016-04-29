function b = isbinary(gm)
% isbinary(fg) : returns true if the model is over only binary variables

% (c) Alexander Ihler 2010

b=~any( gm.dims > 2 );
