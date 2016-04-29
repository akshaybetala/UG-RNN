function b=isfinite(f)
% isnan(f): returns true if factor f has any NAN entries

% (c) Alexander Ihler 2010

b=any(isnan(f.t(:)));

