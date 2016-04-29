function H=entropy(f)
% Compute the entropy in nats of (the normalized version of) a factor f.
% H=entropy(f)
%

% TODO: add conditional entropy support?

% (c) Alexander Ihler 2010

if(0)
t=table(f);
L = log(t);
Z = sum(t(:));
H = t(:).*L(:);
H = -sum(H(~isinf(L(:))));
H = H/Z + log(Z);

else

fnorm=table(normalize(f));
T=fnorm(:).*log(fnorm(:));
%T=fnorm(:).*log2(fnorm(:));
T(fnorm==0)=0;
H = -sum( T );
%H = -sum( fnorm.t(:).*log2(fnorm.t(:)) );
end;
