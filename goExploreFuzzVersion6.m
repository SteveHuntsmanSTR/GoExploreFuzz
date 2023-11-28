function history = goExploreFuzzVersion6(initialStates_X,initialStates_Y,...
    landmarkIndex,gamma,phi,dis_X,dis_data_Y,update_Y,...
    K,globalGenerator,evalBudget,localGenerator,maxExploreEffort,varargin)

% Variant of Go-Explore suitable for (e.g.) greybox *uzzing time-consuming
% things like programs that make expensive function calls, cyber-physical
% systems that incur delays due to interactions, etc. This is intended to
% basically _be_ a general purpose greybox *zzer for this regime that
% nicely separates concerns and leverages all of the information that seems
% like it ought to be both generically available and generically useful.
%
% Cf. https://arxiv.org/abs/2211.12337 and its corresponding core function
% goExploreDissimilarity.m for additional details (including examples of
% most inputs not tied to fuzzing) and context.
%
% We suggest passing the flags "entropic" and "noPareto" in conjunction
% through varargin: see below.
%
% Unlike goExploreDissimilarity, this does not build a surrogate for an
% objective, and it also requires a representation of the objective as a
% composition of two functions, along with dissimilarities on the domains
% of both of these two functions. This is because a greybox *zzer already
% (implicitly) requires this. Also in contrast to goExploreDissimilarity.m,
% this uses a pure exploration mechanism based on magnitude(s) alone. The
% reasoning here is that if we already have the sort of decent
% understanding about how a program behaves in a "neighborhood" of some
% input that a surrogate-based exploitation mechanism would use, then we
% don't need to fuzz there in the first place.
%
% Another difference with goExploreDissimilarity is that the initial states
% and landmarks are taken as arguments rather than internally produced.
% This separation of concerns also means that the arguments L and T in
% goExploreDissimilarity are dropped. One key reason for this is to allow
% the output of one run to be used to start a subsequent run.
%
% Inputs:
%     initialStates_X,
%               Cell array of initial states s.t. dis_Y composed with gamma
%               is nonzero on pairs of distinct states, as produced by
%                     [initialStates_X,initialStates_Y,...
%                         landmarkIndex,~,~] = generateLandmarksFuzz(...
%                         gamma,globalGenerator,dis_Y,L,T);
%               Note that to produce T initial states it may be necessary
%               to run globalGenerator >> T times in this process.
%     initialStates_Y,
%               cellfun(gamma,initialStates_X). Since gamma is expensive,
%               we memorialize this. 
%     landmarkIndex,
%               Index into initialStates_* that yields landmarks.
%     gamma,    Function handle for (e.g.) a computationally expensive
%               map from a space X of "genotypes" to a space Y of
%               "phenotypes" or "traces". The motivating example for
%               greybox *zzing is the map from a cyber(-physical system)
%               input to the corresponding (possibly truncated) path in a
%               control flow graph.
%     phi,      Function handle for (e.g.) a computationally cheap
%               "trace to fitness" objective. Two examples for *zzing are
%               the y-coordinate of a "layered" drawing of the CFG (to get
%               "deep" within the code) or a constant (to get coverage in
%               the sense of diverse vertices, paths, etc. via dis_Y). 
%               TO MAKE phi CONSTANT IN A WAY THAT THE CODE CAN UNDERSTAND,
%               USE phi such that func2str(phi) matches the pattern
%                     "@("+alphanumericsPattern+")"+digitsPattern
%               as with, e.g., @(x) 0, @(abc123) 456, etc. Note also that
%               we can encourage small inputs or short paths in this
%               construct.
%     NB. f in goExploreDissimilarity.m corresponds to @(x) phi(gamma(x)).
%     dis_X,    Function handle for dissimilarity on X (this
%               should be cheap). An example for *zzing is a dissimilarity
%               that is based on a grammar for generating inputs. This can
%               be produced by representing the production rules as a tree
%               and building a corresponding (ultra)metric, possibly with
%               augmentations for continuous parameters. 
%     dis_data_Y,    
%               Function handle for parametrized dissimilarity on Y of the
%               form
%                     dis_data_Y(data_Y,y1,y2)
%               where data_Y is updated per-epoch using update_Y (see just
%               below). The motivating example for *zzing is (a lift to
%               sets of vertices [if not to arcs, paths, etc.] of) a metric
%               on a Markov chain corresponding to the empirical random
%               walk on (a strongly connected modification of) the CFG:
%               here, data_Y would include the observed transitions.
%     update_Y, Function handle for updating "dynamical" data on Y a la
%                     data_Y = update_Y(data_Y,newData)
%               where update_Y encodes some ambient "structure" such as
%               (e.g.) a CFG, and newData is of the form history(j).state_Y
%               or []. Initially, both data_Y and newData MUST be [], so
%               update_Y must handle this initialization requirement
%               gracefully. The archetypal example of arc traversals takes
%               ROUGHLY the following form for D a (suitable strongly
%               connected digraph related to the) CFG and newData
%               corresponding to vertex indices for a path:
%                     update_Y = @(data_Y,newData) data_Y+...
%                         sparse(newData(1:end-1),newData(2:end),1,...
%                         size(D.Nodes,1),size(D.Nodes,1));
%               Note that this requires setting data = 0 initially vs,
%               e.g., data = [], and it assumes numel(newData) > 1, which
%               together mandate using a non-anonymous function that
%               gracefully handles all of this.
%     K,        Rank cutoff <= numel(landmarks).
%     globalGenerator,  
%               Function handle for global state generator.
%     evalBudget, 
%               Number of evaluations to perform.
%     localGenerator,  
%               Function handle for exploration subroutine.
%     maxExploreEffort,  
%               Number bounding exploration effort, as measured in number
%               of evaluations of f per "expedition."
% Optional inputs:
%     Flags passed in as "modern" MATLAB strings (i.e., using double
%     quotes). Possible flags include
%         "entropic"
%             for a power schedule along the lines of
%             https://doi.org/10.1145/3368089.3409748 instead of a
%             maximum-diversity (for phi constant) or heuristic improvement
%             (for phi nonconstant) approach as in goExploreDissimilarity.m
%         "exploitToExplore"
%             trade from objective to diversity contributions for "going"
%         "exploreToExploit"
%             trade from diversity to objective contributions for "going"
%         "goUniform"
%             to use the uniform distribution for "going"
%         "maxDiversity"
%             to use the maximum diversity distribution for "going"
%         "noPareto"
%             to avoid downselection of probes using biobjective Pareto
%             domination
%         "parallel"
%             to execute the gamma evalutions in parallel
%         "pureExploit"
%             to use the objective contributions alone for "going"
%         "pureExplore"
%             to use the diversity contributions alone for "going"
%         "simtropic"
%             for a similarity-sensitive variant of a power schedule along
%             the lines of https://doi.org/10.1145/3368089.3409748 instead
%             of a maximum-diversity (for phi constant) or heuristic
%             improvement (for phi nonconstant) approach as in
%             goExploreDissimilarity.m
%         "unitBandwidth"
%             to enforce unit bandwidth during probing
%         "scaleZero"
%             to compute approximations of maximum diversity constructions
%             at scale zero (note that this is generally the best we can do
%             without devolving into NP-hard constructions).
%     At present, with the obvious exception of "parallel", many of these
%     flags are intended/were used for benchmarking purposes. 
% 
%     Our experiments suggest using "entropic," which should not be
%     surprising ("simtropic" actually does seem less justified than
%     "entropic" when unpacking the latter's rationale). The
%     "unitBandwidth" flag in particular should be avoided absent some good
%     reason. Surprisingly, however, "noPareto" seems to have effects that
%     range from salutary to negligible (and never significantly
%     deleterious), so we suggest using this flag as well.
%
% Output: 
%     history,	
%               Struct with fields
%                   * state_X     i.e., an input to gamma
%                   * cell        the cell containing the state
%                   * birth       the "epoch" in which the state was "born"
%                   * reign       last epoch in which the state was elite
%                   * state_Y     cellfun(@gamma,state_X)
%                   * objective   the objective value.
%                   * source      the cell from which the state was reached
%
% NB. Most of the expensive operations in this code (function evaluations
% and perhaps dis_Y evaluations, though not cutoff scales) could be easily
% performed in parallel with minor modifications. See in particular the
% code comments "*** CAN USE parfor HERE ***" which are not exhaustive.
%
% Last modified 20230912 by Steve Huntsman
%
% Copyright (c) 2023, Systems & Technology Research. All rights reserved.

%% For checking flags
stringIn = @(str,b) ...
    cell2mat(cellfun(@(a)ismember(a,str),b,'UniformOutput',false));

%%
landmarks_X = initialStates_X(landmarkIndex);
landmarks_Y = initialStates_Y(landmarkIndex);

%% Initialize data_Y and dis_Y
data_Y = update_Y([],[]);
for j = 1:numel(initialStates_Y)
    data_Y = update_Y(data_Y,initialStates_Y{j});
end
dis_Y = @(y1,y2) dis_data_Y(data_Y,y1,y2);  % this will vary
% Memorialize an immutable copy for use with stateCell
dis_Y_cell = dis_Y;                         % this will not vary

%% Assign states to cells
% Note that we are sort of overloading MATLAB terminology for cellArray:
% this is a byproduct of sticking with Go-Explore's "cell" terminology
cellArray = stateCell(dis_Y_cell,landmarks_Y,K,initialStates_Y); % matrix
% Get unique cells to encode a set per se of inhabited cells--all
% implicitly, to save the expense of working with the cell representation.
% Use 'stable' option just so we can append efficiently later if any future
% code refactoring suggests it
[~,~,cellNumber] = unique(cellArray,'rows','stable');

%% Initialize history with states used to generate landmarks
% The history (to be instantiated shortly) will detail all the states for
% which the objective is ever evaluated, their corresponding cells, the
% "epoch" in which they were "born" (i.e., the outer loop iteration during
% which the objective was evaluated on them), the last epoch in which they
% were elite (0 if never elite), and their objective values
%
% Note that we can use this information to (cheaply) reconstruct, e.g.,
%     * the number K of closest landmarks used to define cells: this is
%         size([history.cell],1)
%     * the set of inhabited cells: up to ordering/transpose, this is
%         [history([history.reign]==max([history.reign])).cell]
epoch = 1;
% objective *** CAN USE parfor HERE *** BUT there may be no need
objective = nan(size(initialStates_X));
for j = 1:numel(objective)
    objective(j) = phi(initialStates_Y{j});
end
% Reign
reign = zeros(size(initialStates_Y));
for j = 1:max(cellNumber)
    inCell = find(cellNumber==j);  
    [~,indArgmin] = min(objective(inCell));
    reign(inCell(indArgmin)) = epoch;
end
% The instantiation of history as a struct has the practical benefit
% that--with the possible exception of states, which might not be (column)
% vectors in general, but often will be in practice--all of the fields of
% history can be realized as matrices a la [history.cell], etc. In the
% general case, even the states can be realized as cell arrays using
% braces
initialCells = ...
    mat2cell(cellArray',size(cellArray,2),ones(1,size(cellArray,1)));
initialBirth = num2cell(epoch*ones(size(initialStates_X)));
initialReign = num2cell(reign);
initialObjective = num2cell(objective);
initialSource = initialCells;
initialHistory = [initialStates_X;initialCells;initialBirth;...
    initialReign;initialStates_Y;initialObjective;initialSource];
history = cell2struct(initialHistory,...
    ["state_X","cell","birth","reign","state_Y","objective","source"],1);

%% Boolean indicating if phi is (specified as) constant
% I.e., of the form 
%     @(foo) bar 
% for foo an alphanumeric variable and bar a constant made of just digits,
% e.g. @(x) 0, @(abc123) 456, etc. 
pat = "@("+alphanumericsPattern+")"+digitsPattern;
phiIsConstant = matches(func2str(phi),pat);

%% Main loop
evalCount = numel(history);
% Use evalCount in lieu of numel(history) per se in order to better respect
% the evaluation budget
while evalCount < evalBudget
        
    %% Identify elites (a/k/a the "archive" in Go-Explore paper-speak)
    % We keep the entire history to enhance exploration and in mind of the
    % expense of evaluating objectives for our applications
    isElite = [history.reign]==epoch;
    elite_X = {history(isElite).state_X};
    elite_Y = {history(isElite).state_Y};
    elite_cell = {history(isElite).cell};   % don't use []: this is simpler

    %% Begin a new epoch
    epoch = epoch+1;

    %% <GO>
    
    %% Compute dissimilarity matrices for elites
    dissimilarityMatrix_X = zeros(nnz(isElite),nnz(isElite));
    dissimilarityMatrix_Y = zeros(nnz(isElite),nnz(isElite));
    % Symmetry of dis_* is assumed (as otherwise the connection between
    % diversity saturation and weightings breaks down); hence symmetry is
    % built into dissimilarityMatrix_*
    for j1 = 1:nnz(isElite)
        for j2 = (j1+1):nnz(isElite)
            dissimilarityMatrix_X(j1,j2) = dis_X(elite_X{j1},elite_X{j2});
            dissimilarityMatrix_Y(j1,j2) = dis_Y(elite_Y{j1},elite_Y{j2});
        end
    end
    dissimilarityMatrix_X = ...
        max(dissimilarityMatrix_X,dissimilarityMatrix_X');
    dissimilarityMatrix_Y = ...
        max(dissimilarityMatrix_Y,dissimilarityMatrix_Y');

    %% Make educated guess about quantization of dis_X
    % For adjusting bandwidth during probing
    isQuantized_X = ...
        all(dissimilarityMatrix_X==ceil(dissimilarityMatrix_X),'all');  
    
    %% Compute diversity-maximizing distribution (or proxy) on elites
    if any(stringIn("scaleZero",varargin))
        w = dissimilarityMatrix_Y\ones(size(dissimilarityMatrix_Y,1),1);
        w = max(w,0);
    else
        % Compute an appropriate cutoff scale to ensure a bona fide
        % diversity-maximizing distribution on elites. Don't make any
        % assumptions about positive definiteness.
        t = strongCutoff(dissimilarityMatrix_Y)*(1+sqrt(eps));
        Z = exp(-t*dissimilarityMatrix_Y);
        % Solve for weighting while handling degeneracies 
        if max(abs(Z-ones(size(Z))),[],'all') < eps^.75 % sqrt(eps) too big
            w = ones(size(Z,1),1);
        else
            w = Z\ones(size(Z,1),1);
        end
        if any(w<0)
            warning(['min(w) = ',num2str(min(w)),' < 0: adjusting post hoc']); 
            w = w-min(w); 
        end
    end
    maxDiversityPdf = w/sum(w);
    
    %% Construct "go distribution" goPdf balancing diversity and objective
    % The improvement of the objective values should drive progress instead
    % of (e.g.) a regularization coefficient a la temperature in simulated
    % annealing: i.e., this is only implicitly dynamic
    eliteObjectiveValue = [history(isElite).objective];
    % Match top and middle quantiles of logarithm of maximum diversity PDF
    % and elite objective values en route to producing the "go
    % distribution." Note that matching ranges or moments is complicated by
    % the fact that maxDivPdf has one entry that is basically zero, so its
    % logarithm is hard to normalize otherwise
    denom = max(log(maxDiversityPdf))-median(log(maxDiversityPdf));
    denom = denom+(denom==0);
    diversityTerm = (log(maxDiversityPdf)-median(log(maxDiversityPdf)))/...
        denom;
    denom = max(eliteObjectiveValue)-median(eliteObjectiveValue);
    denom = denom+(denom==0);
    objectiveTerm = (eliteObjectiveValue-median(eliteObjectiveValue))/...
        denom;
    if any(stringIn("goUniform",varargin))
        goPdf = ones(size(diversityTerm));
    elseif any(stringIn("pureExplore",varargin))
        goPdf = maxDiversityPdf;
    elseif any(stringIn("exploreToExploit",varargin))
        foo = evalCount/evalBudget;
        goPdf = exp((1-foo)*diversityTerm(:)-foo*objectiveTerm(:));
    elseif any(stringIn("exploreToExploitAdaptive",varargin))
        % "Inverse temperature" beta goes from infinity to zero
        % This might conceivably cause problems with coupon collection,
        % e.g. a negative value for expeditionsThisEpoch below, but handle
        % that in situ, as indicated by a <SafeExpedition> tag
        foo = evalCount/evalBudget;
        beta = -log(2)/log(1-foo);
        goPdf = exp((1-foo)*diversityTerm(:)-foo*beta*objectiveTerm(:));
    elseif any(stringIn("exploitToExplore",varargin))
        foo = evalCount/evalBudget;
        goPdf = exp(foo*diversityTerm(:)-(1-foo)*objectiveTerm(:));
    elseif any(stringIn("exploitToExploreAdaptive",varargin))
        % "Inverse temperature" beta goes from infinity to zero.
        % This might conceivably cause problems with coupon collection,
        % e.g. a negative value for expeditionsThisEpoch below, but handle
        % that in situ, as indicated by a <SafeExpedition> tag
        foo = evalCount/evalBudget;
        beta = -log(2)/log(1-foo);
        goPdf = exp(foo*diversityTerm(:)-(1-foo)*beta*objectiveTerm(:));
    elseif any(stringIn("pureExploit",varargin))
        foo = evalCount/evalBudget;
        goPdf = exp(-foo*objectiveTerm(:));
    else
        % Encourage high diversity contributions and/or low objective vals
        goPdf = exp(diversityTerm(:)-objectiveTerm(:));
    end
    % <Added 20221221> after observing a Dirac distribution in practice
    % Correct for oversparsity or infinities 
    isInf = (goPdf==Inf);
    if nnz(~isInf) < ceil(numel(goPdf)/2)
        % To play nicely with partialCouponCollection invocation
        goPdf = isInf+1/nnz(~isInf);
    elseif any(isInf)
        goPdf = isInf+sqrt(eps);
    end
    % </Added 20221221>
    % <SafeExpedition; added 20230616>
    if any(~isfinite(goPdf))
        warning('any(~isfinite(goPdf)); setting goPdf to uniform');
        goPdf = ones(size(goPdf));  % looking for NaN in particular
    end
    % </SafeExpedition>
    if numel(goPdf) == 1, goPdf = 1; end	% preclude possibility of a NaN 
    goPdf = goPdf/sum(goPdf);	% get a bona fide PDF
    goCdf = cumsum(goPdf);      % CDF for easy sampling
    
    %% Determine "go effort" expeditionsThisEpoch of sampling from goPdf    
    % Use lower bound for collecting half of the coupons from goPdf. The
    % rationale here is to systematically avoid dedicating lots of effort
    % to cells that aren't promising, but to have confidence that many
    % cells will still be explored
    [~,lowerBound,~] = partialCouponCollection(goPdf,...
        ceil(nnz(goPdf)/2),nnz(goPdf)); % scaleZero needed nnz vs numel
    expeditionsThisEpoch = ceil(lowerBound);
    % <SafeExpedition; added 20230616> Have seen some choices of goPdf go
    % awry (viz., the "adaptive" ones, though this doesn't preclude other
    % possibilities), so this establishes some extra guardrails. In the
    % future it would be better to enforce more comprehensive checks on
    % goPdf than the SafeExpedition ones currently in place.
    if expeditionsThisEpoch < 1
        warning('expeditionsThisEpoch < 1; adjusting to |E|*log(|E|)');
        expeditionsThisEpoch = ceil((numel(goPdf)/2)*log(numel(goPdf)/2));
    elseif expeditionsThisEpoch > evalBudget || ~isfinite(expeditionsThisEpoch)
        warning(['expeditionsThisEpoch > evalBudget || '...
            '~isfinite(expeditionsThisEpoch); adjusting to evalBudget/10']);
        expeditionsThisEpoch = ceil(evalBudget/10);
    end
    % </SafeExpedition>
    disp(['evalCount = ',num2str(evalCount),...
        '; expeditionsThisEpoch = ',num2str(expeditionsThisEpoch),...
        '; evalBudget = ',num2str(evalBudget)]);
    
    %% </GO>
    
    %% Memorialize extrema of objective for normalization in loop below
    globalMax = max([history.objective]);
    globalMin = min([history.objective]);
    
    %% Inner loop over expeditions
    newState_X = [];
    newSource = [];
    for expedition = 1:expeditionsThisEpoch
        disp(['expedition = ',num2str(expedition),'/',...
            num2str(expeditionsThisEpoch)]);
        
        %% Sample from goCdf to get a "base elite" from which to explore
        % Recall that entries of goCdf (are presumed to) correspond to
        % inhabited cells
        baseIndex = find(rand<goCdf,1,'first');
        
        %% <EXPLORE>

        %% Determine exploration effort, AKA "power schedule"
        % If the 'entropic' argument is passed in, this is done along the
        % lines of ENTROPIC as described in
        % https://doi.org/10.1145/3368089.3409748, but using the
        % similarity-senstive generalization of (q = 1, i.e.) Shannon
        % entropy (at cellular resolution) instead of vanilla Shannon
        % entropy. Otherwise a maximum diversity Ansatz is used if phi is
        % "obviously" constant, or the approach of goExploreDissimilarity
        % if phi is not obviously constant. 
        baseCell = elite_cell{baseIndex}'; % note transpose
        if any(stringIn("simtropic",varargin))
            % Get probability of going from base cell to other cells
            fromBaseCell = ...
                [history(ismember([history.source]',baseCell,'rows')).cell];
            % Recall that
            %     target(ind,:) = fromBaseCell'
            [target,~,ind] = unique(fromBaseCell','rows');
            % Recall that
            %     num(j) = nnz(ismember(fromBaseCell',target(j,:),'rows'))
            num = histcounts(ind);
            % Ensure statistical and numerical safety by adding (initial)
            % ones to all counts
            p = ones(1,nnz(isElite));
            for j = 1:size(target,1)
                ind_j = ismember(cell2mat(elite_cell)',target(j,:),'rows');
                p(ind_j) = p(ind_j)+num(j);
            end
            p = p/sum(p);
            % Similarity-sensitive entropy (see eqn (6.3) of the Leinster
            % book)
            if any(stringIn("scaleZero",varargin))
                % Deploy a subtle bound here that we proceed to detail:
                %
                % The diversity of order $1$ is
                % \begin{equation}
                % \label{eq:diversity1}
                % D_1^{Z}(p) := \prod_{j:p_j > 0} (Zp)_j^{-p_j}
                % \end{equation}
                % and the corresponding similarity-sensitive generalization 
                % of Shannon entropy is
                % \begin{equation}
                % \label{eq:ssEntropy1}
                % \log D_1^{Z}(p) = - \sum_{j:p_j > 0} p_j \log (Zp)_j.
                % \end{equation}
                % We want to compute the maximum value of 
                % \eqref{eq:ssEntropy1} for the common case $Z = \exp[-td]$ 
                % in the limit $t \downarrow 0$ for a generic dissimilarity 
                % $d$. However, this goal turns out to be overly ambitious.
                % 
                % The first-order approximation 
                % $Z = \exp[-td] \approx 11^T - td$ yields
                % \begin{equation}
                % \label{eq:ssEntropy1Approx}
                % \log D_1^{Z}(p) \approx t p^T d p.
                % \end{equation}
                % If $Z$ is positive definite for all sufficiently small 
                % $t$ (e.g., when $d$ is a Euclidean distance matrix) then 
                % by Theorem 3.2 of \cite{izumino2006maximization} we have 
                % that 
                % $$\arg \max_{p \in \mathbb{R}^n : 1^T p = 1} p^T d p = 
                % \frac{d^{-1}1}{1^T d^{-1} 1},$$
                % though in general this extremum will have negative 
                % components, so the best we can practically do is to bound 
                % \eqref{eq:ssEntropy1Approx} using
                % \begin{equation}
                % \label{eq:bound1}
                % \max_{p \ge 0 : 1^T p = 1} p^T d p 
                % \le \max_{p \in \mathbb{R}^n : 1^T p = 1} p^T d p 
                % = \frac{1}{1^T d^{-1} 1}.
                % \end{equation}
                % It can be shown (see, e.g., Example 5.16 of 
                % \cite{devriendt2022graph}) that maximizing 
                % \eqref{eq:ssEntropy1Approx} is $\mathbf{NP}$-hard for 
                % arbitrary $d$. There are strong reasons to believe this 
                % remains true even when $Z$ is positive definite for all 
                % sufficiently small $t$. However, in the more restrictive
                % case where $d$ is ultrametric this maximization can be 
                % performed efficiently \cite{leinster2012measuring,
                % leinster2016maximizing,leinster2021entropy} and more 
                % generally when $d$ is \emph{negative type} in the sense 
                % of \cite{leinster2017magnitude}, where a convex quadratic 
                % program provides a solution: see Proposition 5.20 of 
                % \cite{devriendt2022graph}.
                % 
                % To summarize, when $Z = \exp[-td]$ is positive definite 
                % for all sufficiently small $t$, we have the bound
                % \begin{align}
                % \label{eq:bound2}
                % \lim_{t \downarrow 0} 
                % \frac{\log D_1^Z(p)}
                % {\max_{p \ge 0 : 1^T p = 1} \log D_1^Z(p)} 
                % & \ge \lim_{t \downarrow 0} 
                % \frac{\log D_1^Z(p)}
                % {\max_{p \in \mathbb{R}^n : 1^T p = 1} \log D_1^Z(p)} 
                % \nonumber \\
                % & = (p^T d p) \cdot (1^T d^{-1} 1).
                % \end{align}
                % 
                % This bound governs the power schedule of our fuzzer at 
                % scale $t=0$.
                % \footnote{
                % Proposition 5.17 of \cite{devriendt2022graph} gives the 
                % generic bounds $p^T d p \in 
                % [\frac{1}{2},\frac{n-1}{n}] \cdot \max_{jk} d_{jk}$.
                % }
                % % \cite{izumino2006maximization} refers to 
                % % Izumino, S. and Nakamura, N. "Maximization of quadratic
                % forms expressed by distance matrices." 
                % % \cite{devriendt2022graph} refers to Devriendt, K. 
                % "Graph geometry from effective resistances."
                foo = p(:)'*dissimilarityMatrix_Y*p(:);
                bar = ones(size(p(:)'))*dissimilarityMatrix_Y*ones(size(p(:)));
                ssEntropyNorm = foo*bar;
                if ssEntropyNorm > 1
                    warning(['ssEntropyNorm = ',num2str(ssEntropyNorm),...
                        ' > 1; changing to 1 (hopefully very nearby)']);
                    ssEntropyNorm = 1;
                end
            else
                ssEntropy = -p*log(Z*p(:));
                % Divide by max possible value: at p = w, it's log(sum(w))
                ssEntropyNorm = ssEntropy/log(sum(w));
            end
            exploreEffort = ceil(maxExploreEffort*ssEntropyNorm);
        elseif any(stringIn("entropic",varargin))
            % Get probability of going from base cell to other cells
            fromBaseCell = ...
                [history(ismember([history.source]',baseCell,'rows')).cell];
            % Recall that
            %     target(ind,:) = fromBaseCell'
            [target,~,ind] = unique(fromBaseCell','rows');
            % Recall that
            %     num(j) = nnz(ismember(fromBaseCell',target(j,:),'rows'))
            num = histcounts(ind);
            % Ensure statistical and numerical safety by adding (initial)
            % ones to all counts
            p = ones(1,nnz(isElite));
            for j = 1:size(target,1)
                ind_j = ismember(cell2mat(elite_cell)',target(j,:),'rows');
                p(ind_j) = p(ind_j)+num(j);
            end
            p = p/sum(p);
            % Similarity-sensitive entropy (see eqn (6.3) of the Leinster
            % book)
            entropy = -p*log(p(:));
            % Divide by max possible value: log(numel(p)) for p uniform
            entropyNorm = entropy/log(numel(p));
            exploreEffort = ceil(maxExploreEffort*entropyNorm);
        else
            if phiIsConstant
                % Exploration effort determined by diversity contribution
                divFactor = maxDiversityPdf/max(maxDiversityPdf);
                exploreEffort = ceil(maxExploreEffort*divFactor(baseIndex));
            else
                % Get history of the base cell
                baseCell = ...
                    stateCell(dis_Y_cell,landmarks_Y,K,elite_Y(baseIndex));
                inBase = find(ismember([history.cell]',baseCell,'rows'));        
                % Determine expeditions that visit the base cell (easy),
                % rather than expeditions that start from the base cell
                % (hard).
                %
                % Note that our method of recording history (viz., state_*,
                % cell, birth, reign, and objective) makes it impossible to
                % exactly reconstruct the itinerary of expeditions (i.e.,
                % which cells hosted bases during a given epoch). While we
                % could guess at this--or, more sensibly, augment our
                % records to know this--it's not clear that this would
                % actually be useful, and it would increase the complexity
                % of the code. So we avoid this
                baseEpoch = unique([history(inBase).birth]);
                % In base cell from last expedition visiting it
                oneAgo = find([history(inBase).birth]==baseEpoch(end));
                % In base cell from penultimate expedition visiting it
                if numel(oneAgo) && numel(baseEpoch)>1
                    twoAgo = find([history(inBase).birth]==baseEpoch(end-1));
                else
                    twoAgo = oneAgo;	% could be empty in principle
                end
                % Determine normalized differential in objective over last
                % two expeditions visiting the base cell
                oneAgoObjective = [history(inBase(oneAgo)).objective];
                twoAgoObjective = [history(inBase(twoAgo)).objective];
                denom = globalMax-globalMin;
                denom = denom+(denom==0);
                bestOneAgoNormed = min((oneAgoObjective-globalMin)/denom);
                bestTwoAgoNormed = min((twoAgoObjective-globalMin)/denom);
                baseDelta = bestOneAgoNormed-bestTwoAgoNormed;  % in [-1,1]
                % Determine exploration effort on the basis of prior
                % efforts
                if epoch > 2
                    priorExploreEffort = ...
                        nnz([history(inBase).birth]==baseEpoch(end));
                else
                    % Geometric mean of 1 and maxExploreEffort seems
                    % appropriate
                    priorExploreEffort = ceil(sqrt(maxExploreEffort));
                end
                foo = max(priorExploreEffort*2^-baseDelta,1);
                exploreEffort = ceil(min(foo,maxExploreEffort));
            end
        end

        %% Probe, iteratively reducing bandwidth as appropriate
        disp("probing");
        numProbes = 2*maxExploreEffort;
        probe = cell(1,numProbes);
        if any(stringIn("unitBandwidth",varargin))
            bandwidth = 1;
            for j = 1:numel(probe)
                probe{j} = localGenerator(elite_X{baseIndex},bandwidth);
            end
        else
            % Initial bandwidth is big--we will bring it down to size. Note
            % that this works just fine for Booleans and Hamming distance
            % (or its ilk) along with the localGenerator function of
            % Example 3 above--more generally, one can generally retool the
            % localGenerator function to play nicely with this
            % initialization. <BANDWIDTH>
            bandwidth = max(dissimilarityMatrix_X(baseIndex,:));
            for j = 1:numel(probe)
                % In principle we could make the localGenerator function
                % interact with available information--for example, taking
                % the in-cell and/or nearby history with objective values
                % below some quantile, getting the mean and covariance of
                % that data, and sampling from a Gaussian with those same
                % parameters. But we aren't doing that now/yet, mainly
                % because this overspecializes.
                probe{j} = localGenerator(elite_X{baseIndex},bandwidth);
            end
            % Compare probe cells with base cell IN X SPACE
            baseCell_X = stateCell(dis_X,landmarks_X,K,elite_X(baseIndex));
            probeCellArray = stateCell(dis_X,landmarks_X,K,probe); % matrix
            probeInCell = ismember(probeCellArray,baseCell_X,'rows');
            % Reduce bandwidth until a significant number of probes remain
            % in the current cell (cf. Gaussian annulus theorem);
            disp("probing: adjusting bandwidth");
            while nnz(probeInCell) < numProbes/4
                bandwidth = bandwidth/2;
                for j = 1:numel(probe)
                    probe{j} = localGenerator(elite_X{baseIndex},bandwidth);
                end
                probeCellArray = stateCell(dis_X,landmarks_X,K,probe);
                probeInCell = ismember(probeCellArray,baseCell_X,'rows');
                if isQuantized_X && bandwidth < 1
                    % Avoid indefinitely reducing bandwidth while probing
                    % if dis_X is something like Hamming distance
                    break;
                end
            end
        end
        
        %% Ensure uniqueness of probes
        % Out of an abundance of caution, avoid probe collisions (this has
        % never happened in practice if problems weren't absurdly small or
        % localGenerator functions weren't poorly constructed). 
        disp("probing: checking uniqueness");
        % Build a graph G whose vertices correspond to probes and whose
        % edges correspond to pairs for which dis_X == 0. The connected
        % components of G will represent distinct probes. This is all done
        % with a local function.
        probe = unique_X(probe,dis_X);
        % It's possible that probe might intersect with {history.state_X}
        bar = zeros(numel(history),numel(probe));
        for j1 = 1:numel(history)
            for j2 = 1:numel(probe)
                bar(j1,j2) = dis_X(history(j1).state_X,probe{j2});
            end
        end
        probeInHistory = any(bar==0);
        if all(probeInHistory)
            % Go global
            probe = cell(1,numProbes);
            for j = 1:numel(probe)
                probe{j} = globalGenerator();
            end
            % Avoid collisions
            probe = unique_X(probe,dis_X);
        end

        %% Ensure nondegeneracy of probes
        % This turns out to be necessary
        probe = probe(~cellfun(@isempty,probe));    

        %% Pick probes to evaluate
        if any(stringIn("noPareto",varargin))
            toEvaluate = probe;
        else
            %% Biobjective data: weighting on probes and diff. magnitude
            disp('biobjective data: batch & individual probe diversities');
            % [Big comment follows]
            %
            % Since we are not performing the sort of biobjective analysis
            % of GoExploreDissimilarity via an approximation of the
            % objective, we might as well perform a *different* biobjective
            % analysis. As in GoExploreDissimilarity, a local dissimilarity
            % matrix considers all the probes at once. Meanwhile, the
            % differential magnitude of each individual probe relative to
            % elites is also considered. That is, we are considering the
            % diversity contributions of probes both relative to each other
            % and relative to preexisting elites.
            %
            % To understand this, try running the following:
            %
            %     rng('default');
            %     dim = 2;
            %     N = 100;
            %     x = rand(dim,N);
            %     %%
            %     d_x = squareform(pdist(x'));
            %     t_x = posCutoff(d_x)*(1+sqrt(eps));
            %     Z_x = exp(-t_x*d_x);
            %     w_x = Z_x\ones(size(Z_x,1),1);
            %     %%
            %     n = 10;
            %     y = rand(dim,n);
            %     %%
            %     d_y = squareform(pdist(y'));
            %     t_y = posCutoff(d_y)*(1+sqrt(eps));
            %     Z_y = exp(-t_x*d_y);
            %     w_y = Z_y\ones(size(Z_y,1),1);
            %     %% Differential magnitude
            %     diffMag = nan(size(y,2),1);
            %     for j = 1:numel(diffMag)
            %         % E.g., d_x(j,:) = vecnorm(x-x(:,j))
            %         foo = vecnorm(x-y(:,j));
            %         zeta = exp(-t_x*foo(:));
            %         diffmag(j) = (1-zeta.'*w_x)^2/(1-(zeta.')*(Z\zeta));
            %     end
            %     %%
            %     figure;
            %     scatter(x(1,:),x(2,:),10,w_x/max(w_x),'filled');
            %     hold on;
            %     scatter(y(1,:),y(2,:),30,diffmag/max(diffmag),...
            %         'filled','MarkerEdgeColor','k');
            %     text(y(1,:)+.02,y(2,:)+.02,string(1:size(y,2)));
            %     colorbar;
            %     daspect([1,1,1]);
            %     axis([0,1,0,1]);
            %     %%
            %     figure; 
            %     plot(w_y/max(w_y),diffmag/max(diffmag),'ko');
            %     hold on;
            %     text(w_y/max(w_y)+.02,diffmag/max(diffmag)+.02,...
            %         string(1:size(y,2)));
            %     axis([0,1,0,1]);
            %     xlabel('normalized local weighting component');
            %     ylabel('normalized differential magnitude');
            
            %% Weighting on probes
            local = probe;
            dissimilarityMatrix_local = zeros(numel(local),numel(local));
            try
                for j1 = 1:numel(local)
                    for j2 = 1:numel(local)
                        dissimilarityMatrix_local(j1,j2) = ...
                            dis_X(local{j1},local{j2});
                    end
                end
            catch
                assignin('base','local',local);
                assignin('base','base',elite_X(baseIndex));   
                assignin('base','base_Y',elite_Y(baseIndex));   
                error('assignin');
            end
            % Compute an appropriate cutoff scale to ensure a bona fide
            % diversity-maximizing distribution on elites. Don't make any
            % assumptions about positive definiteness.
            t_local = strongCutoff(dissimilarityMatrix_local)*(1+sqrt(eps));
            % Local similarity matrix and weighting: latter is proportional
            % to diversity-maximizing distribution
            Z_local = exp(-t_local*dissimilarityMatrix_local);
            w_local = Z_local\ones(size(Z_local,1),1);
    
            %% Differential magnitude
            % See appendix C of https://arxiv.org/abs/2201.11677 for
            % details
            if any(stringIn("scaleZero",varargin))
                % At scale zero, the magnitude is always unity, so the
                % differential magnitude must be zero
                diffMag = zeros(numel(probe),1);
            else
                diffMag = nan(numel(probe),1);
                for j = 1:numel(probe)
                    foo = nan(numel(elite_X),1);
                    for k = 1:numel(elite_X)
                        foo(k) = dis_X(probe{j},elite_X{k});
                    end
                    zeta = exp(-t_local*foo);
                    diffMag(j) = (1-zeta.'*w)^2/(1-(zeta.')*(Z\zeta));
                end
            end
            
            %%
            biObjective = [-diffMag,-w_local]';
    
            %% Determine quantitative Pareto dominance
            % We will select the least-dominated points after normalizing
            % both objectives to zero mean and unit variance a la
            biObjectiveNormalized = diag(var(biObjective,0,2))\...
                (biObjective-mean(biObjective,2));
            % we want to rank points by their Pareto domination a la
            dominatedBy = nan(size(biObjectiveNormalized,2),1);
            for j = 1:numel(dominatedBy)
                % This is slicker than repmat
                obj_j = biObjectiveNormalized(:,j)...
                    *ones(1,size(biObjectiveNormalized,2));
                dominator = min(obj_j-biObjectiveNormalized);
                dominatedBy(j) = max(dominator);
            end
    
            %% Determine new states to evaluate after this loop
            cutoff = numel(local)-numel(probe);
            % 1:cutoff corresponds to stuff already in history
            [~,dominanceInd] = sort(dominatedBy((cutoff+1):end));
            % Taking a minimum here to prevent an out-of-bounds error was
            % necessary on rare occasion (precisely once in dozens of uses
            % before instituted), for reasons that aren't immediately
            % obvious: no effort has yet been made to understand this
            numToEvaluate = min(exploreEffort,numel(dominanceInd));
            toEvaluate = probe(dominanceInd(1:numToEvaluate));
        end
        % Update in a way that lets us terminate on budget
        newState_X = [newState_X,toEvaluate]; %#ok<AGROW>
        evalCount = evalCount+numel(toEvaluate);

        %%
        newSource = [newSource,...
            repmat({baseCell'},1,numel(toEvaluate))]; %#ok<AGROW> 
        
        %% Respect the evaluation bound
        if evalCount > evalBudget
            evalExcess = evalCount-evalBudget;
            newState_X = newState_X(1:(numel(newState_X)-evalExcess));
            newSource = newSource(1:(numel(newSource)-evalExcess));
            break; 
        end
        
        %% </EXPLORE>
    end
       
    %% Apply gamma
    % *** CAN USE parfor HERE ***
    state_Y = cell(size(newState_X));
    if any(stringIn("parallel",varargin))
        parfor j = 1:numel(state_Y)
            % gamma incurs broadcast overhead but no obvious way around
            % this. To determine if this is a problem in any specific
            % instance see
            % https://www.mathworks.com/matlabcentral/answers/339728
            state_Y{j} = gamma(newState_X{j}); %#ok<PFBNS> 
        end
    else
        for j = 1:numel(state_Y)
            state_Y{j} = gamma(newState_X{j});
            disp(['evalCount (gamma): ',num2str(numel(history)+j)]); 
        end
    end

    %% Assign states to cells (cf. similar code @ initialization)
    disp('assigning states to cells'); 
    cellArray = [cellArray;...
        stateCell(dis_Y_cell,landmarks_Y,K,state_Y)]; %#ok<AGROW> 
    [~,~,cellNumber] = unique(cellArray,'rows','stable');

    %% Update data_Y and dis_Y
    disp('updating data_Y and dis_Y');
    for j = 1:numel(state_Y)
        data_Y = update_Y(data_Y,state_Y{j});
    end
    dis_Y = @(y1,y2) dis_data_Y(data_Y,y1,y2);
    
    %% Prepare history updates (cf. similar code @ initialization)
    state_X = [{history.state_X},newState_X];
    birth = [[history.birth],epoch*ones(size(newState_X))];
    % objective *** CAN USE parfor HERE *** BUT there may be no need
    objective = nan(size(state_Y));
    for j = 1:numel(objective)
        objective(j) = phi(state_Y{j});
    end
    state_Y = [{history.state_Y},state_Y]; %#ok<AGROW>
    objective = [[history.objective],objective]; %#ok<AGROW>
    source = [{history.source},newSource]; 

    %% Reign
    reign = [[history.reign],zeros(size(newState_X))];
    for j = 1:max(cellNumber)
        inCell = find(cellNumber==j);
        [~,indArgmin] = min(objective(inCell));
        reign(inCell(indArgmin)) = epoch;        
    end

    %% Update history
    newState_X = state_X;
    newCell = mat2cell(cellArray',...
        size(cellArray,2),ones(1,size(cellArray,1)));   % a cell
    newBirth = num2cell(birth);
    newReign = num2cell(reign);
    newState_Y = state_Y;
    newObjective = num2cell(objective);
    newSource = source;
    newHistory = [newState_X;newCell;newBirth;newReign;newState_Y;...
        newObjective;newSource];
    history = cell2struct(newHistory,...
        ["state_X","cell","birth","reign","state_Y","objective",...
        "source"],1);

end

end

%% END MAIN FUNCTION ******************************************************

%% Local function
function probe = unique_X(probe,dis_X)

% Produces a subset of probes with pairwise nonzero dissimilarities

d_uX = zeros(numel(probe),numel(probe));
for i1 = 1:numel(probe)
    for i2 = (i1+1):numel(probe)
        d_uX(i1,i2) = dis_X(probe{i1},probe{i2});
    end
end
d_uX = max(d_uX,d_uX');
G_uX = graph(d_uX==0);
% Keep only the first representative per connected component of G
cc_uX = conncomp(G_uX,'OutputForm','cell');
probe = probe(cellfun(@(x)x(1),cc_uX));

end