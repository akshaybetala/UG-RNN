function b=isfinite(f)
% isfinite(f) : returns false if factor f has any infinite/nan values

% (c) Alexander Ihler 2010

b=all(isfinite(f.t(:)));

