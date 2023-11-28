function path = dynamicExecution(cfg,S,b,numSteps,x) %#ok<INUSD> 

% Dynamically execute the program specified by cfg, S, and b for numSteps
% on input x. Cf. fuzzableCFG.m

%% S and b are not in this code but are asserted to be in the workspace
assert(exist('S','var'),"exist('S','var')");
assert(exist('b','var'),"exist('b','var')");

%%
eval(['x = ',mat2str(x),';']); %#ok<EVLEQ> 
%
currentNode = "% START: 0";
path = currentNode;
for j = 1:numSteps
    if startsWith(currentNode,"% HALT")
        break;
    end
    foo = cfg.Nodes.Predicate(findnode(cfg,currentNode));
    bar = cfg.Nodes.Assignment(findnode(cfg,currentNode));
    if ~ismissing(foo)
        predicate = eval(foo);
        % Go down appropriate branch and update
        [edgeId,nodeId] = outedges(cfg,currentNode);
        branch = find(predicate==cfg.Edges.Weight(edgeId));
        currentNode = nodeId(branch); %#ok<FNDSB>
    elseif ~ismissing(bar)
        eval(join([bar,";"],"")); % avoid spamming stdout
        [~,currentNode] = outedges(cfg,currentNode);
        assert(numel(currentNode)==1,"elseif");
    else
        [~,currentNode] = outedges(cfg,currentNode);
        assert(numel(currentNode)==1,"else");
    end
    path = [path,currentNode]; %#ok<AGROW>
end