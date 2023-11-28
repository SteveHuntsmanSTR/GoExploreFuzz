function cellId = stateCell(dis,landmarks,K,x)

% Produce cell identifiers for states
%
% Inputs:
%     dis,  function handle for dissimilarity (assumed symmetric)
%     landmarks,  
%           cell array of landmarks (versus a larger cell array of states 
%           and landmark indices)
%     K,    rank cutoff <= numel(lan) = L
%     x,    state(s) to map to cell
%
% Output: 
%     cellId,	array of size [numel(x),K] whose ith row is a sorted tuple
%               of the K closest landmarks (as indices, i.e., a subindex of
%               landmark indices)
%
% NB. To work with states and landmarkIndex as produced by
%     generateLandmarks, simply preclude executing this function by the
%     assignment
%         landmarks = states(landmarkIndex);
%
% Last modified 20220425 by Steve Huntsman
%
% Copyright (c) 2022, Systems & Technology Research. All rights reserved.

%% Check scalar input
% Other inputs are too tricky to check meaningfully here
if ~isscalar(K), error('K not scalar'); end
if ~isfinite(K), error('K not finite'); end
if ~isreal(K), error('K not real'); end
if K ~= round(K), error('K not integral'); end
if K < 1, error('K < 1'); end
L = numel(landmarks);
if K > L, error('K > L'); end

%%
dissimilarityFromLandmarks = nan(numel(x),L);
for i = 1:numel(x)
    for ell = 1:L
        dissimilarityFromLandmarks(i,ell) = dis(x{i},landmarks{ell});
    end
end
[~,ind] = sort(dissimilarityFromLandmarks,2);
cellId = ind(:,1:K);