function D = functionalEditDistance(s,t,c,d)

% Weighted edit distance between arrays s and t using function handles c
% and d to respectively compute the costs for deletion/insertion and
% substitution. Symmetry of D requires symmetry of d, as well as a single
% cost for deletion/inserion.
%
% s and t can be arrays of any class, but must be treatable as
% one-dimensional. c and d must be compatible with s and t.
%
% This uses the Wagner-Fischer algorithm. In particular, the dynamic
% programming matrix is never actually populated: instead, we simply update
% the current column using the previous one.
%
% Example: 
%     % Unit costs all around yields Levenshtein distance
%     c = @(s) numel(s);
%     d = @(s,t) ~isequal(s,t);
%     D = functionalEditDistance('SUNDAY','SATURDAY',c,d) % yields 3
%
% Last modified 20210702 by Steve Huntsman
%
% Copyright (c) 2021, Systems & Technology Research. All rights reserved.

%% Ensure arrays are nominally compatible and make them 1D
% Error checking actual compatibility of arrays would have to be recursive;
% incorporating functions would only be more involved. We avoid these
% tarpits. Caveat emptor.
if ~ismatrix(s), error('s not an array'); end
if ~ismatrix(t), error('t not an array'); end
if ~isequal(class(s),class(t)), error('s and t not of the same class'); end
s = s(:)';  % row better for introspection/debugging
t = t(:)';  % ibid.

%% Initializations
m = numel(s);
n = numel(t);
% Precompute deletion/insertion costs for both inputs
del = nan(m,1);
for j = 1:m
    del(j) = c(s(j));
end
ins = nan(1,n);
for k = 1:n
    ins(k) = c(t(k));
end
% At most two columns of the dynamic programming matrix are ever in memory
prevCol = [0;cumsum(del)];

%% Main loop
for k = 1:n
    currCol = zeros(size(prevCol));
    currCol(1) = prevCol(1)+ins(k);
    for j = 1:m
        sub = d(s(j),t(k)); % substitution cost
        if isequal(s(j),t(k)) && sub
            warning('sub ~= 0: fixing');
            sub = 0;
        end
        if sub < 0, error('sub < 0'); end
        % Update current column of dynamic programmming matrix
        currCol(j+1) = min([currCol(j)+del(j),...
            prevCol(j+1)+ins(k),...
            prevCol(j)+sub]);
    end
    prevCol = currCol;
end

%% Output
D = currCol(end);