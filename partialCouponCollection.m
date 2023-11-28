function [exact,lowerBound,upperBound] = partialCouponCollection(p,m,c)

% Calculation of the expected time for the general coupon collection
% problem in which the ith kind of coupon is drawn with probability p(i),
% and m is the number of coupon types to collect. Using p = ones(1,m)/m
% thus recovers the classical result for the uniform case. The exact result
% is computed when reasonably cheap (for details, see Corollary 4.2 of
% https://doi.org/10.1016/0166-218X(92)90177-C; for the case m = numel(p),
% see Theorem 4.1), and bounds computed more generally based on the same
% approach (both for rapidly decreasing and approximately uniform p).
%
% Inputs:
%     p,      probability distribution
%     m,      number of coupons to collect
%     c,      cutoff parameter for computations (will compute 2^c terms
%             provided that c <= 16: otherwise, this will be avoided--this
%             can speed things up dramatically)
%
% Outputs:
%     exact expected time (= NaN if unknown) and lower/upper bounds
%
% For evaluation purposes it is useful to avoid taking exact results as
% bounds; see code cells starting with
%         %% OMIT THIS FOR EVALUATION PURPOSES
% and comment them out as warranted.
%
% By way of (internal) documentation, see the appended LaTeX snippet.
%
% NB. Note that gammaln is used instead of nchoosek except when actually
% producing (vs counting) the combinations. This approach is MUCH faster. 
%
% NB. Subsumes couponCollectionTotal.m and generalCouponCollection.m via
% special cases--these are thus deprecated.
%
% NB. To get just the lowerBound output argument anonymously, use something
% like this (for the second output argument):
% 
%     function nthArg = outargn(fun,n)%#ok<STOUT,INUSL> 
% 
%     % Return nth output argument from function fun. If [y_1,...,y_n,...,z]
%     % = fun(inputs) then y_n = outargn(@() fun(inputs),n). I think this is
%     % a tolerable use of eval.
% 
%     eval(['[',repmat('~,',1,n-1),'nthArg] = fun();']);
%
% Last modified 20220502 by Steve Huntsman
%
% Copyright (c) 2022, Systems & Technology Research. All rights reserved.

%% Check p matrix, finite, real, nonnegative, sums to unity
if ~ismatrix(p), error('p not matrix'); end
isReal = isfinite(p)&isreal(p);
if ~all(isReal(:)), error('p not real'); end
if any(p<0), error('p not nonnegative'); end
if abs(sum(p)-1) > sqrt(eps)
    warning('p does not sum to unity: normalizing');
end
p = p(:)'/sum(p);

%% Check m
if ~isscalar(m), error('m not scalar'); end
if ~isfinite(m), error('m not finite'); end
if ~isreal(m), error('m not real'); end
if m ~= round(m)
    error('m not integral');
end
if m < 1, error('m < 1'); end

%% Check c
if ~isscalar(c), error('c not scalar'); end
if ~isfinite(c), error('c not finite'); end
if ~isreal(c), error('c not real'); end
if c ~= round(c)
    error('c not integral');
end
if c < 1, error('c < 1'); end
if c > 16
    avoid_c = 1;
    warning('c too big to quickly compute corresponding bounds: avoiding');
else
    avoid_c = 0;
end

%% Check number of coupon types is achievable (note restriction to support)
if m > nnz(p), error('numCouponTypes > nnz(p)'); end

%% Annoying corner case
if m == 1
    warning('m = 1');
    exact = 1;
    lowerBound = 1;
    upperBound = 1;
    return;
end

%%
p = sort(p,'descend');
n = numel(p);
exact = NaN;

%% [EXACT] Produce exact result for both bounds if easy, else m = n bound
if n <= 16
    exact = 0;
    for ell = 0:(m-1)
        ind = nchoosek(1:n,ell);
        pInd = reshape(p(ind),size(ind));
        sign = (-1)^(m-1-ell);
        S = sum(1./(1-sum(pInd,2)));
        exact = exact+sign*nchoosek(n-ell-1,n-m)*S;
    end
    %% OMIT THIS FOR EVALUATION PURPOSES
    lowerBound = exact;
    upperBound = exact;
    return;
end

%% [TOTAL] Get "total" result for m = n, which is also an upper bound
% Theorem 4.1 of https://doi.org/10.1016/0166-218X(92)90177-C yields
% expectation for total coupon collection as an integral that admits
% straightforward numerical computation
% disp('Computing for case m = n');
f = @(t) 1-prod(1-exp(-p(:)*t),1);
x = 1; while f(x) > eps, x = 10*x; end  % OK upper limit for integral
t = linspace(0,x,1e4);
total = trapz(t,f(t));

%% Initialize bounds
lowerBound = 0;
upperBound = total;

%% Return total result if m = n
if m == n
    exact = total;
    %% OMIT THIS FOR EVALUATION PURPOSES
    lowerBound = exact;
    upperBound = exact;
    return;
end

%% Now we'll have to compute less trivial bounds (c permitting)...
% See notes 
if ~avoid_c
    %% Form power set of 1:c
    % disp('Generating power set of 1:c');
    powerSet = false(2^c,c);
    for j = 1:c
        powerSet(:,j) = bitget((0:(2^c-1))',j);
    end
    %% Bounds assuming rapid decay of p
    % Based on exact result of Corollary 4.2 in
    % https://doi.org/10.1016/0166-218X(92)90177-C
    % disp('Computing bounds for rapidly decaying p');
    lb0 = 0;
    ub0 = 0;
    for ell = 0:(m-1)
        % disp(['    ',num2str(ell+1),'/',num2str(m)]);
        lambda = min(c,ell);
        localPowerSet = powerSet(1:2^lambda,1:lambda);
        %%
        sum_ell_lb = 0;
        sum_ell_ub = 0;
        for j = 1:size(localPowerSet,1)
            M = localPowerSet(j,:);
            P_M = sum(p(M));
            mu = nnz(M);
            ind = 1:(ell-mu);
            if ell-mu >= 0 && ell-mu <= n-lambda
                % Indices for ell-mu smallest and largest components of p
                % that aren't already reserved for P_M
                ind_lb = ind+n-(ell-mu);
                ind_ub = ind+min(lambda,n-(ell-mu));
                % numer = nchoosek(n-lambda,ell-mu);
                % Much faster to use gammaln than nchoosek here; no warnings
                numer = exp(gammaln((n-lambda)+1)-gammaln((ell-mu)+1)...
                    -gammaln((n-lambda)-(ell-mu)+1));
                sum_ell_lb = sum_ell_lb+numer/(1-P_M-sum(p(ind_lb)));
                sum_ell_ub = sum_ell_ub+numer/(1-P_M-sum(p(ind_ub)));
            end
        end
        %%
        % coeff = nchoosek(n-ell-1,n-m)*(-1)^(m-1-ell);
        % Much faster to use gammaln than nchoosek here; no warnings
        foo = exp(gammaln((n-ell-1)+1)-gammaln((n-m)+1)...
            -gammaln((n-ell-1)-(n-m)+1));
        coeff = foo*(-1)^(m-1-ell);
        lb0 = lb0+min(coeff*sum_ell_lb,coeff*sum_ell_ub);
        ub0 = ub0+max(coeff*sum_ell_lb,coeff*sum_ell_ub);
    end
    lowerBound = max(lowerBound,lb0);
    upperBound = min(upperBound,ub0);
    %% Bounds assuming near-uniformity of p
    % Based on exact result of Corollary 4.2 in
    % https://doi.org/10.1016/0166-218X(92)90177-C. Similar mechanically to
    % preceding stuff.
    % disp('Computing bounds for nearly uniform p');
    lbu = 0;
    ubu = 0;
    delta = n*p-1;  % deviation from uniformity
    for ell = 0:(m-1)
        % disp(['    ',num2str(ell+1),'/',num2str(m)]);
        lambda = min(c,ell);
        localPowerSet = powerSet(1:2^lambda,1:lambda);
        %%
        sum_ell_lbu = 0;
        sum_ell_ubu = 0;
        for j = 1:size(localPowerSet,1)
            M = localPowerSet(j,:);
            Delta_M = sum(p(M));
            mu = nnz(M);
            % Indices for ell-mu largest magnitude deviations that aren't
            % already reserved for Delta_M
            ind = (1:(ell-mu))+min(lambda,n-(ell-mu));
            if ell-mu >= 0 && ell-mu <= n-lambda
                denomL = 1-(ell+Delta_M+sum(abs(delta(ind))))/n;
                denomU = 1-(ell+Delta_M-sum(abs(delta(ind))))/n;
                denomL = max(denomL,eps);   % to be safe
                denomU = max(denomU,eps);   % to be safe
                % numer = nchoosek(n-lambda,ell-mu);
                % Much faster to use gammaln than nchoosek here; no warnings
                numer = exp(gammaln((n-lambda)+1)-gammaln((ell-mu)+1)...
                    -gammaln((n-lambda)-(ell-mu)+1));
                sum_ell_lbu = sum_ell_lbu+numer/denomL;
                sum_ell_ubu = sum_ell_ubu+numer/denomU;
            end
        end
        %%
        % coeff = nchoosek(n-ell-1,n-m)*(-1)^(m-1-ell);
        % Much faster to use gammaln than nchoosek here; no warnings
        foo = exp(gammaln((n-ell-1)+1)-gammaln((n-m)+1)...
            -gammaln((n-ell-1)-(n-m)+1));
        coeff = foo*(-1)^(m-1-ell);
        lbu = lbu+min(coeff*sum_ell_lbu,coeff*sum_ell_ubu);
        ubu = ubu+max(coeff*sum_ell_lbu,coeff*sum_ell_ubu);
    end
    lowerBound = max(lowerBound,lbu);
    upperBound = min(upperBound,ubu);
end

%% Augment lower bound using uniform case
% This is an easy calculation from the Corollary 4.2 cited above: the fact
% that the uniform case provides a lower bound is both intuitively obvious
% and proved in https://doi.org/10.1239/jap/1437658606.
%
% (NB. It is elementary to show that the gradient of the expectation w/r/t
% coupon probabilities is zero at uniformity, and similarly that the
% Hessian is diagonal there [note that this tactic is not employed by the
% reference cited here in favor of a global argument].)
uniformLowerBound = n*(sum(1./(1:n))-sum(1./(1:(n-m))));
lowerBound = max(lowerBound,uniformLowerBound);

%% LaTeX documentation
% Per Corollary 4.2 of \cite{flajolet1992birthday}, we have that the
% expected time for the event $X_m$ of collecting $m$ of $n$ coupons via
% IID draws from the distribution $(p_1,\dots,p_n)$ satisfies
% \begin{equation} \label{eq:partialCoupon} \mathbb{E}(X_m) =
% \sum_{\ell=0}^{m-1} (-1)^{m-1-\ell} \binom{n-\ell-1}{n-m} \sum_{|L| =
% \ell} \frac{1}{1-P_L} \end{equation} with $P_L := \sum_{k \in L} p_k$.
% However, the sum \eqref{eq:partialCoupon} is generally difficult or
% impossible to evaluate in practice due to its combinatorial complexity,
% and it is desirable to produce useful bounds. \footnote{ The specific
% case $m = n$ admits an integral representation that readily admits
% numerical computation, viz. $\mathbb{E}(X_n) = \int_0^\infty \left ( 1 -
% \prod_{k=1}^n [1-\exp(-p_k t)] \right ) \ dt$. While an integral
% representation of $\mathbb{E}(X_m)$ also exists for generic $m$, it is
% also combinatorial in form and \eqref{eq:partialCoupon} (which is
% actually just the result of evaluating it symbolically) appears easier to
% compute. }
% 
% Towards this end, assume w/l/o/g that $p_1 \ge \dots \ge p_n$, and let $c
% \le n$. (For clarity, it is helpful to imagine that $c < n$ and $p_c \gg
% p_{c+1}$, but we do not assume this.) To bound $\sum_{|L| = \ell}
% (1-P_L)^{-1}$, we first note that $\{L : |L| = \ell\}$ is the union of
% disjoint sets of the form $\{L : |L| = \ell \text{ and } L \cap [\lambda]
% = M\}$ for $M \in 2^{[\lambda]}$, where $\lambda := \min \{c,\ell\}$.
% Thus \begin{equation} \label{eq:bound1a} \sum_{|L| = \ell}
% \frac{1}{1-P_L} = \sum_{M \in 2^{[\lambda]}} \sum_{\substack{|L| = \ell
% \\ L \cap [\lambda] = M}} \frac{1}{1-P_L}. \end{equation}
% 
% Now $P_L = P_{L \cap [\lambda]} + P_{L \backslash [\lambda]}$. If we are
% given bounds of the form $\pi_* \le P_{L \backslash [\lambda]} \le
% \pi^*$, we get in turn that $$\frac{1}{1-P_{L \cap [\lambda]}-\pi_*} \le
% \frac{1}{1-P_L} \le \frac{1}{1-P_{L \cap [\lambda]}-\pi^*}.$$ If
% furthermore $\pi_*$ and $\pi^*$ depend on $L$ only via $L \cap
% [\lambda]$, then \begin{align} \label{eq:bound1b} \frac{|\{L : |L| = \ell
% \text{ and } L \cap [\lambda] = M\}|}{1-P_M-\pi_*} & \le
% \sum_{\substack{|L| = \ell \\ L \cap [\lambda] = M}} \frac{1}{1-P_L}
% \nonumber \\ & \le \frac{|\{L : |L| = \ell \text{ and } L \cap [\lambda]
% = M\}|}{1-P_M-\pi^*}. \end{align} Meanwhile, writing $\mu := |M|$ and
% combinatorially interpreting the Vandermonde identity $\sum_\mu
% \binom{\lambda}{\mu} \binom{n-\lambda}{\ell-\mu} = \binom{n}{\ell}$
% yields \begin{equation} \label{eq:bound1c} |\{L : |L| = \ell \text{ and }
% L \cap [\lambda] = M\}| = \binom{n-\lambda}{\ell-\mu} \end{equation} and
% in turn bounds of the form \begin{equation} \label{eq:bound1d} \sum_{M
% \in 2^{[\lambda]}} \binom{n-\lambda}{\ell-\mu} \frac{1}{1-P_M-\pi_*} \le
% \sum_{|L| = \ell} \frac{1}{1-P_L} \le \sum_{M \in 2^{[\lambda]}}
% \binom{n-\lambda}{\ell-\mu} \frac{1}{1-P_M-\pi^*}. \end{equation}
% 
% Now the best possible choice for $\pi_*$ is
% $P_{[\ell-\mu]+n-(\ell-\mu)}$; similarly, the best possible choice for
% $\pi^*$ is $P_{[\ell-\mu]+\min\{\lambda,n-(\ell-\mu)\}}$. This
% immediately yields upper and lower bounds for \eqref{eq:partialCoupon},
% though the alternating sign term leads to intricate expressions that are
% hardly worth writing down explicitly.
% 
% The resulting bounds are hardly worth using in some situations, and quite
% good in others. We augment them with the easy lower bound obtained by
% using the uniform distribution in \eqref{eq:partialCoupon}
% \cite{anceaume2015new} and the easy upper bound obtained by taking $m =
% n$; we also use the exact results when feasible (e.g., $n$ small or $m =
% n$) as both upper and lower bounds. These basic augmentations have a
% significant effect in practice.
% 
% Experiments on exactly solvable (in particular, small) cases show that
% though the bounds for $(1-P_L)^{-1}$ are good, the combinatorics involved
% basically always obliterates the overall bounds for distributions of the
% form $p_k \propto k^{-\gamma}$ with $\gamma$ a small positive integer.
% However, the situation improves dramatically for distributions that decay
% quickly enough.
% 
% We can similarly also derive bounds along the lines above based on the
% deviations $\delta_k := n p_k - 1$. The only significant difference in
% the derivation here versus the one detailed above is that we are forced
% to consider absolute values of the deviations. In this regime we also
% have the simple and tight lower bound $\mathbb{E}(X_m) \ge n(H_n -
% H_{n-m})$, where $H_n := \sum_{k=1}^n k^{-1}$ is the $n$th harmonic
% number \cite{flajolet1992birthday,anceaume2015new}. In fact this bound is
% quite good for $n$ small, to the point that replacing harmonic numbers
% with logarithms can easily produce larger deviations from the bound than
% the error itself.