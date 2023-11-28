function data_Y = temp20230508_update(D,data_Y,state_Y)

%% Update counts
if isempty(data_Y)
    if ~isempty(state_Y), error('~isempty(state_Y)'); end
    data_Y.counts = adjacency(D);  % every arc gets a +1 initially
else
    src = state_Y(1:end-1);
    tar = state_Y(2:end);
    n = size(D.Nodes,1);
    data_Y.counts = data_Y.counts+sparse(src,tar,1,n,n);
end

%% Update metric based on counts
P = diag(sum(data_Y.counts,2))\data_Y.counts;   % row-stochastic matrix
% Pass in a second argument to get_dhp that is slightly larger than the 0.5
% value that can (and in practice will) produce pseudometric degeneracies
beta = 1/2+1e-2;                                % 0.51 seems good
dhp = get_dhp(P,beta);  % no second arg means beta = 1/2
data_Y.groundMetric = dhp;