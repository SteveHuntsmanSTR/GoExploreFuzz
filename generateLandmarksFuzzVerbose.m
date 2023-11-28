function [initialStates_X,initialStates_Y,allInitial_X,allInitial_Y,...
    landmarkIndex,magnitude,counter] = ...
    generateLandmarksFuzzVerbose(gamma,globalGenerator,dis_Y,L,T)

% Variant of generateLandmarksFuzz used for benchmarking. Key addition is
% allInitial_*, which contains all the initial states generated, regardless
% of collisions.
%
% Generate diverse landmark states/points (we use the term state by
% reference to a dissimilarity framework for Go-Explore-type algorithms,
% which this was developed for) and corresponding traces. Differs from
% generateLandmarks.m in that the possibility of collisions with respect to
% dis_Y is accounted for, and traces themselves are computed/retained. 
% 
% Such collisions can happen when (e.g.) the composition
% gamma(globalGenerator()) can produce collisions, as when trying to
% produce program paths for fuzzing purposes, in which case path collisions
% are almost inevitable. This is also possible when globalGenerator itself
% can produce collisions. Even in this case, setting gamma = @(x) x and
% dis_Y = dis_X (= a dissimilarity on the range of globalGenerator) can be
% useful: note however that initialStates_Y will then be redundant.
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
%     L,    number of landmarks (<= T)
%     T,    number of points to generate (including landmarks, so >= L)
%
% Output: 
%     initialStates_X,  cell array of all states produced by
%                       globalGenerator (for reuse) 
%     initialStates_Y,  image of initialStates_X under gamma: entries
%                       guaranteed pairwise distinct by construction if
%                       dis_Y is zero on diagonal
%     allInitial_X,     cell array of all states in X produced without 
%                       collisions removed 
%     allInitial_Y,     cell array of all states in Y produced without 
%                       collisions removed 
%     landmarkIndex,    indices of states that are landmarks 
%     magnitude,        magnitude of landmarks after given generations
%     counter,          number of evals of gamma
%
% Last modified 20230518 by Steve Huntsman
%
% Copyright (c) 2023, Systems & Technology Research. All rights reserved.

%% Check scalar inputs
% Function handles are too tricky to check meaningfully
if ~isscalar(L), error('L not scalar'); end
if ~isfinite(L), error('L not finite'); end
if ~isreal(L), error('L not real'); end
if L ~= round(L), error('L not integral'); end
if L < 1, error('L < 1'); end
if ~isscalar(T), error('T not scalar'); end
if ~isfinite(T), error('T not finite'); end
if ~isreal(T), error('T not real'); end
if T ~= round(T), error('T not integral'); end
if T < L, error('T < L'); end

%% First state
initialStates_X = cell(1,T);
initialStates_Y = cell(1,T);
initialStates_X{1} = globalGenerator();
initialStates_Y{1} = gamma(initialStates_X{1});
allInitial_X{1} = initialStates_X{1};
allInitial_Y{1} = initialStates_Y{1};

%% Form subsequent initial states, removing dissimilarity collisions ASAP
ell = 1;
counter = 1;
while ell < L
    if mod(counter,L) == 0 && counter > 0
        disp([num2str(counter),' gamma evals; ',...
            num2str(ell),'/',num2str(L),' desired unique initial landmarks']);
    end
    %%
    newState_X = globalGenerator();
    newState_Y = gamma(newState_X);
    counter = counter+1;
    allInitial_X{counter} = newState_X; %#ok<AGROW> 
    allInitial_Y{counter} = newState_Y; %#ok<AGROW> 
    disInitialVsNew = zeros(1,ell);
    for j = 1:ell
        disInitialVsNew(j) = dis_Y(initialStates_Y{j},newState_Y);
    end
    if any(disInitialVsNew==0)
        continue;
    else
        ell = ell+1;
        initialStates_X{ell} = newState_X;
        initialStates_Y{ell} = newState_Y;
    end
end
landmarkIndex = 1:L;

%% Initial dissimilarity matrix, scale, weighting, etc.
d0 = zeros(L);
for j = 2:L
    for k = 1:(j-1)
        % There is some superfluous indexing here for illustration
        d0(j,k) = dis_Y(initialStates_Y{landmarkIndex(j)},...
            initialStates_Y{landmarkIndex(k)});
    end
end
d0 = max(d0,d0');
% The t = 0 limit is not worth considering here, as it yields either unit
% magnitude or naughty behavior. For the sake of generality (and because
% it's a one-time cost) we invoke the strong cutoff scale.
t = strongCutoff(d0)*(1+sqrt(eps));
Z0 = exp(-t*d0);
w0 = Z0\ones(size(Z0,1),1);
magnitude = nan(1,T);
magnitude(L) = sum(w0);

%% Main loop
for i = (L+1):T
    %% Propose new state to replace the one with least weighting component
    % But only if there is no collision
    [~,ind] = min(w0);
    disInitialVsNew = 0;
    while any(disInitialVsNew==0)
        if mod(counter,L) == 0 && counter > 0
            disp([num2str(counter),' globalGenerator evals; ',...
                num2str(i-1),'/',num2str(T),...
                ' desired unique initial states, including ',...
                num2str(L),' landmarks']);
        end
        %%
        newState_X = globalGenerator();
        newState_Y = gamma(newState_X);
        allInitial_X{counter} = newState_X; %#ok<AGROW> 
        allInitial_Y{counter} = newState_Y; %#ok<AGROW> 
        counter = counter+1;
        disInitialVsNew = zeros(1,i-1);
        for j = 1:(i-1)
            disInitialVsNew(j) = dis_Y(initialStates_Y{j},newState_Y);
        end
        if any(disInitialVsNew==0)
            continue;
        else
            initialStates_X{i} = newState_X;
            initialStates_Y{i} = newState_Y;
        end
    end
    
    %% Gauge impact on magnitude at original (strong cutoff) scale
    newRow = zeros(1,L);
    for ell = 1:L
        if ell ~= ind
            newRow(ell) = dis_Y(newState_Y,initialStates_Y{landmarkIndex(ell)});
        end
    end
    d1 = d0;
    d1(ind,:) = newRow;
    d1(:,ind) = newRow';
    Z1 = exp(-t*d1);
    w1 = Z1\ones(size(Z1,1),1);
    
    %% Update
    if sum(w1) > magnitude(i-1)
        d0 = d1;
        w0 = w1;
        landmarkIndex(ind) = i;
        magnitude(i) = sum(w1);
    else
        magnitude(i) = magnitude(i-1);
    end
end

%%
counter = counter-1;