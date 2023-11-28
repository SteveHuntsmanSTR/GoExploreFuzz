function t = strongCutoff(d)

% Minimal t such that exp(-u*d) is positive semidefinite and admits a
% nonnegative weighting for any u > t Here d is a symmetric extended real
% matrix with zero diagonal.
%
% Last modified 20210319 by Steve Huntsman
%
% Copyright (c) 2021, Systems & Technology Research. All rights reserved.

%% Check d matrix, square, nonnegative extended real, zero diag
if ~ismatrix(d), error('d not matrix'); end
[m,n] = size(d);
if m ~= n, error('d not square'); end
isExtendedReal = (isinf(d)|isfinite(d))&isreal(d);
if any(any(~isExtendedReal|d<0)), error('d not nonnegative extended real'); end
if any(diag(d)~=0), error('d diagonal not zero'); end
if max(max(abs(d'./d-1))) > sqrt(eps), error('d not symmetric'); end

%%
t = log(n-1)/min(min(d+diag(inf(1,n))));
lower = 0;
upper = t;
while 1-lower/upper>sqrt(eps)
    t = (lower+upper)/2;
    Z = exp(-t*d);
    spec = eig(Z);
    if min(spec)>=0
        w = Z\ones(size(Z,1),1);
        if all(w>0)
            upper = t;
        else
            lower = t;
        end
    else
        lower = t;
    end
end