function [initialStates_X,initialStates_Y,numEvals] = ...
    generateInitialFuzz(gamma,globalGenerator,dis_Y,n)

% Generate states/points (we use the term state by reference to a
% dissimilarity framework for Go-Explore-type algorithms, which this was
% developed for) with distinct corresponding traces (as measured by dis_Y).
% Differs from generateLandmarksFuzz.m in that we only seek distinctness,
% not maximal diversity.
%
% Inputs:
%     gamma,    
%           Function handle for (e.g.) a computationally expensive "state
%           to trace" map. The motivating example for greybox *zzing is the
%           map from a CPS input to the corresponding (possibly truncated)
%           path in a control flow graph.
%     globalGenerator,  
%           function handle for global state generator
%     dis_Y,  
%           function handle for dissimilarity on range of gamma (assumed
%           symmetric)
%     n,    number of points to generate
%
% Output: 
%     initialStates_X,  cell array of all states produced by
%                       globalGenerator (for reuse) 
%     initialStates_Y,  image of initialStates_X under gamma: entries
%                       guaranteed pairwise distinct by construction if
%                       dis_Y is zero on diagonal
%     numEvals,         number of evals of gamma
%
% Last modified 20230515 by Steve Huntsman
%
% Copyright (c) 2023, Systems & Technology Research. All rights reserved.

%% Check scalar inputs
% Function handles are too tricky to check meaningfully
if ~isscalar(n), error('n not scalar'); end
if ~isfinite(n), error('n not finite'); end
if ~isreal(n), error('n not real'); end
if n ~= round(n), error('n not integral'); end
if n < 1, error('n < 1'); end

%% First state
initialStates_X = cell(1,n);
initialStates_Y = cell(1,n);
initialStates_X{1} = globalGenerator();
initialStates_Y{1} = gamma(initialStates_X{1});

%% Form subsequent initial states, removing dissimilarity collisions ASAP
numStates = 1;
numEvals = 1;
increment = ceil(sqrt(n));  % for tracking progress in stdout
while numStates < n
    if mod(numEvals,increment) == 0
        disp([num2str(numEvals),' gamma evals; ',num2str(numStates),'/',...
            num2str(n),' desired unique initial landmarks']);
    end
    %%
    newState_X = globalGenerator();
    newState_Y = gamma(newState_X);
    numEvals = numEvals+1;
    disInitialVsNew = zeros(1,numStates);
    for j = 1:numStates
        disInitialVsNew(j) = dis_Y(initialStates_Y{j},newState_Y);
    end
    if any(disInitialVsNew==0)
        continue;
    else
        numStates = numStates+1;
        initialStates_X{numStates} = newState_X;
        initialStates_Y{numStates} = newState_Y;
    end
end