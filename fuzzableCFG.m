function [cfg,S,b] = fuzzableCFG(pdf,N,n_productions)

% Produce an annotated CFG and associated arrays that dynamicExecution.m
% can consume.
%
% INPUTS: e.g., pdf = [6;1;3;0], N = 256, n_productions = 20
%
%     The entries of pdf correspond to the respective productions
%         S; S | if b; S; fi | while b; S; end | jnz b L
%     so in particular set the last entry to zero to avoid jnz (jump
%     nonzero), i.e., "goto." The code for this has been (and is intended
%     to be) used elsewhere. To be able to _run_ jnz in MATLAB it's
%     necessary to use an approach along the lines of
%     https://www.mathworks.com/matlabcentral/answers/225928-how-do-i-use-goto-in-matlab#answer_1066860
%     which is not implemented here/yet.
%
%     N is the "size of a byte."
%
%     n_productions gives the number of times the BNF grammar above is
%     applied.
%
% OUTPUTS: cfg is the actual CFG; S and b carry data about statements and
% Boolean predicates.

%% Check inputs
% pdf
if ~ismatrix(pdf), error('pdf not matrix'); end
if numel(pdf) ~= 4, error('numel(pdf) ~= 4'); end
isNonnegative = isreal(pdf)&(pdf>=0|isfinite(pdf));
if any(~isNonnegative(:)), error('~any(~isNonnegative(:))'); end
% n_productions
if ~isscalar(n_productions), error('~isscalar(n_productions)'); end
if ~isreal(n_productions), error('~isreal(n_productions)'); end
if ~isfinite(n_productions), error('isfinite(n_productions)'); end
if n_productions ~= round(n_productions)
    error('n_productions ~= round(n_productions)'); 
end

%% 
if pdf(4) > 0, error('pdf(4) > 0: will yield unrunnable code (jnz)'); end

%% Produce program (skeleton) from probabilistic context free grammar: S ->
% S; S | if b; S; fi | while b; S; end | jnz b L
productionSet = ...
    {["S";"S"];["if b";"S";"fi"];["while b";"S";"end"];["jnz b L"]}; %#ok<NBRAK>
% Distribution on productions represented as a map for semantic clarity
keySet = cellfun(@join,productionSet);
pdf = pdf/sum(pdf); % normalize
cdf = cumsum(pdf);
prob = containers.Map(keySet,pdf);
if prob("jnz b L") > 0, warning('Will not produce executable MATLAB!'), end
% Produce program (skeleton) from PCFG
program = "S";
for j = 1:n_productions
    % Choose production and location
    ind_production = find(rand>[0;cdf],1,'last');
    foo = find(program=="S");
    if numel(foo)
        ind_location = foo(randi(numel(foo)));
        % Apply production to location (note program is a row array)
        program = [program(1:(ind_location-1)),...
            productionSet{ind_production}',...
            program((ind_location+1):end)];
    else
        break;  % all terminals (all jnz vs S is technically possible)
    end
end
program = program(:);

%% Handle jnz (= jump nonzero)
% Note that usually this won't be present
jump = find(program=="jnz b L");
LZs = find(ismember(program,["if b","while b","S"]));	% landing zones
if numel(LZs) > 0
    % All jumps go to a LZ
    for j = 1:numel(jump)
        LZ = LZs(randi(numel(LZs)));
        program(jump(j)) = join(["jnz b",LZ]);
    end
else
    % No place to land, so turn jumps into statements (vs., e.g., deleting)
    if numel(jump)
        program(program=="jnz b L") = "S";
    end
end

%% Build control flow graph
% We build a control flow graph by associating lines in a program with
% vertices in the CFG; edges are produced according to the following table:
%     src @ j | tar(1)  | tar(0)
%     ---------------------------
%     if b	  | j+1	    | [fi]+1
%     while b | j+1	    | [end]+1
%     end     | [while] | *
%	  jnz b L | L       | j+1
%     fi; S   | j+1	    | *
% Here [] = line #. Matching is done ridiculously, with multiple passes (so
% fight me, or better still, show me a lightweight LL(1) parser
% implementation in MATLAB that magically lets me bypass the construction
% of an abstract syntax tree like this hack does). Constructs besides if,
% fi, while, end, and jnz are interpreted as statements S. The weight of
% edges corresponds to the truth value of a predicate.
cfg = digraph;
program = ["% START";program(:)];
for j = 1:numel(program)
    if program(j) == "if b"
        tar1 = j+1;
        tar0 = j-1+matchParenthesis(program(j:end),"fi")+1;     % + START
        cfg = addedge(cfg,j,tar0,0);                            % 0 = false
    elseif program(j) == "while b"
        tar1 = j+1;
        tar0 = j-1+matchParenthesis(program(j:end),"end")+1;	% + START
        cfg = addedge(cfg,j,tar0,0);                            % 0 = false
    elseif program(j) == "end"
        tar1 = j-matchParenthesis(program(j:-1:1),"while b")+1;
    elseif startsWith(program(j),"jnz b")
        foo = char(program(j));
        tar1 = str2double(foo(numel(char("jnz b"))+2:end))+1;	% + START
        tar0 = j+1;
        cfg = addedge(cfg,j,tar0,0);                            % 0 = false
    else
        % fi; S
        tar1 = j+1;
    end
    cfg = addedge(cfg,j,tar1,1);                                % 1 = true
end
program = [program(:);"% HALT"];
cfg.Nodes.Name = append(program,": ",string(0:numel(program)-1)');

%% Compute number of predicates from START to a given node
cfg.Nodes.PredicateDepth = nan(size(cfg.Nodes,1),1);
for j = 1:size(cfg.Nodes,1)
    cfg.Nodes.PredicateDepth(j) = ...
        nnz(contains(cfg.Nodes.Name(shortestpath(cfg,1,j)),["if","while"])); 
end

%% Build executable program suitable for fuzzing T&E
toEval = program(2:end-1);      % omit START; HALT
toEval(toEval=="fi") = "end";	% to eventually run within MATLAB
% Populate Boolean predicates and statements with actual content
ind_b = find(endsWith(toEval,"b"));
ind_S = find(toEval=="S");
b = cell(size(toEval));
S = cell(size(toEval));
for j = 1:numel(ind_b)
    b{ind_b(j)} = randi(N);
    foo = char(toEval(ind_b(j)));
    entry = ['x(',num2str(cfg.Nodes.PredicateDepth(ind_b(j)+1)),')'];
    toEval(ind_b(j)) = ...
        string([foo(1:end-1),entry,' == b{',num2str(ind_b(j)),'}']);
end
%
for j = 1:numel(ind_S)
    entry = ['x(',num2str(cfg.Nodes.PredicateDepth(ind_S(j))),')'];
    % NB. mod((1:10)-1,10)+1 = 1:10
    toEval(ind_S(j)) = string([entry,' = mod(',entry,',',num2str(N),')+1']);
end
% Further annotate CFG nodes with eval commands and values to eval
cfg.Nodes.Eval = ["% START";toEval;"% HALT"];
cfg.Nodes.Predicate = regexp(cfg.Nodes.Eval,'x\(\d+\) == .+','match','once');
cfg.Nodes.Assignment = regexp(cfg.Nodes.Eval,'x\(\d+\) = .+','match','once');