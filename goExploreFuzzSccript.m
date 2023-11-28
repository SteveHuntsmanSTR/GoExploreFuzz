%% Be VERY sure before running...
x = input(['This will clear workspace and take days if not weeks to run. ',...
    'Enter y to continue: '],"s");
if ~matches(x,"y")
    error('aborted!');
end

%% Number of control flow graphs; trials per graph; evals per trial
numCFG = 10;    % number of programs/CFGs
numTrials = 40; % number of trials for fuzzing each program/CFG
evalBudget = 1000;

%% Flags for goExploreFuzz* and IDs for objectives
flags = [...
    "default",...                   % not recognized flag, will be ignored
    "goUniform",...                 % go via
    "pureExplore",...               % go via
    "exploreToExploit",...          % go via
    "exploreToExploitAdaptive",...  % go via
    "exploitToExplore",...          % go via
    "exploitToExploreAdaptive",...  % go via
    "pureExploit",...               % go via
    "simtropic",...                 % effort
    "entropic",...                  % effort
    "unitBandwidth",...
    "noPareto",...
    ];
objectives = [...
    "hpDist",...    % hitting prob distance from entry @(y) min(-d1(y))
    "hpDistExp",... % exp of above @(y) min(-exp(d1(y)))
    "ordDist",...   % digraph distance from entry @(y) min(-d1ordinary(y))
    "ordDistExp",...    % exp of above @(y) min(-exp(d1ordinary(y))
    "depth",...     % depth (in 'layered' drawing) @(y) min(yData(y))
    "depthExp",...  % exp of above @(y) min(-exp(max(yData)-yData(y)));
    "const",...     % constant @(y) 1
    ];

%% Parameters for the control flow graphs/programs
pdf = [.6;.1;.3;0]; % last component 0 -> no gotos 
N = 8;              % maximum value of a "byte" plus 1
n_productions = 50; % context free grammar -> control flow graph

%% Produce control flow graphs/programs
rng('default'); 
cfg_m = cell(1,numCFG);
S_m = cell(1,numCFG);
b_m = cell(1,numCFG);
m = 0;
while m < numCFG
    %% Produce control flow graph
    [cfg,S,b] = fuzzableCFG(pdf,N,n_productions);
   
    %% Restrict cfg to its largest strong component D
    % Any leading statements/subroutines and/or if's will be omitted
    scc = conncomp(cfg,'Type','strong');
    hc = histcounts(scc,'BinMethod','integers');
    [~,ind] = max(hc);
    D = subgraph(cfg,scc==ind);

    %% 
    if size(D.Nodes,1) < size(cfg.Nodes,1)-2
        continue;
    else
        m = m+1;
        cfg_m{m} = cfg;
        S_m{m} = S;
        b_m{m} = b;

        %% Plot
        % figure; 
        % plot(cfg,'Layout','layered','NodeLabel',cfg.Nodes.Eval,...
        %     'EdgeCData',cfg.Edges.Weight); 
        % colormap(flipud(jet));

        %% Plot
        % figure; 
        % plot(D,'Layout','layered','NodeLabel',D.Nodes.Eval,...
        %     'EdgeCData',D.Edges.Weight); 
        % colormap(flipud(jet));
    end
end

%% Characterize goExploreFuzz and global generator coverage
meanCoverage = cell(numCFG,numel(flags),numel(objectives));
stdCoverage = cell(numCFG,numel(flags),numel(objectives));
meanCoverage_gg = cell(numCFG,1);
stdCoverage_gg = cell(numCFG,1);
timePerCFG = nan(1,numCFG);
for m = 1:numCFG
    tic;
    cfgStr = ['CFG ',num2str(m),'/',num2str(numCFG)];

    %% Redo stuff from above
    cfg = cfg_m{m};
    S = S_m{m};
    b = b_m{m};

    %% Restrict cfg to its largest strong component D
    % Any leading statements/subroutines and/or if's will be omitted
    scc = conncomp(cfg,'Type','strong');
    hc = histcounts(scc,'BinMethod','integers');
    [~,ind] = max(hc);
    D = subgraph(cfg,scc==ind);

    %% For depth objective
    figure;
    plo = plot(D,'Layout','layered');
    yData = plo.YData;
    close;

    %% Ordinary digraph distance from (neighbor in D of) start vertex
    foo = distances(D);
    d1ordinary = foo(1,:);
    
    %% Hitting probability distance from (neighbor in D of) start vertex
        % function data_Y = temp20230508_update(D,data_Y,state_Y)
        % 
        % %% Update counts
        % if isempty(data_Y)
        %     if ~isempty(state_Y), error('~isempty(state_Y)'); end
        %     data_Y.counts = adjacency(D);  % every arc gets a +1 initially
        % else
        %     src = state_Y(1:end-1);
        %     tar = state_Y(2:end);
        %     n = size(D.Nodes,1);
        %     data_Y.counts = data_Y.counts+sparse(src,tar,1,n,n);
        % end
        % 
        % %% Update metric based on counts
        % P = diag(sum(data_Y.counts,2))\data_Y.counts;   % row-stochastic 
        % % Pass in a second argument to get_dhp that is slightly larger
        % % than the 0.5 value that can (and in practice will) produce
        % % pseudometric degeneracies
        % beta = 1/2+1e-2;                                % 0.51 seems good
        % dhp = get_dhp(P,beta);  % no second arg means beta = 1/2
        % data_Y.groundMetric = dhp;    
    data_Y = temp20230508_update(D,[],[]);
    d1 = data_Y.groundMetric(1,:);
    
    %% Inputs to goExploreFuzz* EXCEPT THE OBJECTIVE phi (see below)
    % Dynamically execute a program by traversing its CFG for numCFGevals
    % steps, then toss out any vertices in the path that are not in its
    % largest strong component D. By construction, such vertices can only
    % occur at the beginning and/or end of the path.
    numCFGevals = 40;
    gamma0 = @(x) dynamicExecution(cfg,S,b,numCFGevals,x);
    nonzero = @(x) x(x~=0);
    gamma = @(x) nonzero(findnode(D,gamma0(x)));
    % Dissimilarity on X is basically Hamming distance. An edit distance
    % would be absurdly slow in practice and overly slow even for
    % evaluation purposes. One might however consider an optimized
    % implementation of relative de Bruijn entropy as an approximation: see
    % https://arxiv.org/abs/1509.02975 (note that as a practical matter
    % this would require using [e.g.] a globalGenerator that produced
    % structs with one field a sequence and another field the output of
    % radixwordquiver on that sequence, but optimizing the implementation
    % would be the focus of effort here).
    minLen = @(x1,x2) min(numel(x1),numel(x2));
    maxLen = @(x1,x2) max(numel(x1),numel(x2));
    hamming = @(x1,x2) nnz(x1(1:minLen(x1,x2))~=x2(1:minLen(x1,x2)))...
        +maxLen(x1,x2)-minLen(x1,x2);
    dis_X = hamming;
    % Use edit distance on CFG paths with insert/delete cost =
    % max(data_Y.groundMetric,[],'all') and substitution cost =
    % data_Y.groundMetric.
    % 
    % A cheaper alternative would be Hausdorff distance on traces lifted 
    % from data_Y.groundMetric a la
    %     hausMat = @(d,y1,y2) d(unique(y1(:)),unique(y2(:))); dis_data_Y =
    %     @(data_Y,y1,y2) max(...
    %         max(min(hausMat(data_Y.groundMetric,y1,y2),[],1)),...
    %         max(min(hausMat(data_Y.groundMetric,y1,y2),[],2)));
    % However, so long as dis_Y = @(y1,y2) dis_data_Y(data_Y,y1,y2) is 
    % cheap relative to gamma, we should be OK. In practice we would want 
    % (and be able) to engineer this by only considering certain 
    % vertices/basic blocks.
    dis_data_Y = @(data_Y,y1,y2) functionalEditDistance(y1,y2,...
        @(s)max(data_Y.groundMetric,[],'all'),...
        @(s,t)data_Y.groundMetric(s,t));
    % ...and the ground metric is the hitting probabilities metric obtained
    % from the empirical Markov chain on D
    update_Y = @(data_Y,state_Y) temp20230508_update(D,data_Y,state_Y);
    % To get as many elites as possible, set K = L
    L = 6;
    T = 8;
    K = L;  
    % Uniform strings on 1:N of length equal to predicate depth (which by
    % construction is sufficient and sometimes necessary)
    globalGenerator = ...
        @() randi(N,1,max(cfg.Nodes.PredicateDepth));
    % localGenerator
        % function y = temp20230511(x,theta,n)
        % 
        % % a localGenerator that overwrites symbol positions uniformly
        % % sampled (with replacement) in a suffix determined by the 
        % % realized position of an initial symbol overwrite. Cf. 
        % % temp20230418.m
        % 
        % y = x;
        % ind = randi([1,numel(y)]);
        % y(ind) = randi(n);
        % for j = 2:ceil(theta)
        %     ind_j = randi([ind,numel(y)]);
        %     y(ind_j) = randi(n);
        % end
    localGenerator = @(x,theta) temp20230511(x,theta,N);
    %
    maxExploreEffort = 2*N;   % AFL uses 512 = 2*256: here, N is like 256
    
    %%
    numEdges = size(D.Edges,1); % NOT size(cfg.Edges,1)

    %% Common initialization to be a bit more efficient
    parfor i = 1:numTrials
        triStr = ['trial ',num2str(i),'/',num2str(numTrials)];
        disp(['*** initialization: ',cfgStr,'; ',triStr]);
        %% Generate landmarks
        data_Y = temp20230508_update(D,[],[]);
        dis_Y = @(y1,y2) dis_data_Y(data_Y,y1,y2);
        seed = i*1e6;
        rng(seed);  % for baseline w/ only globalGenerator below
        [initialStates_X{i},initialStates_Y{i},...
            allInitial_X{i},allInitial_Y{i},...
            landmarkIndex{i},~,counter{i}] = ...
            generateLandmarksFuzzVerbose(gamma,globalGenerator,dis_Y,L,T);
    end

    %% Go exploring, looping over flags and objectives
    for i_flag = 1:numel(flags)
        flag = flags(i_flag);
        flagStr = ['flag ',num2str(i_flag),'/',num2str(numel(flags))];
        for i_obj = 1:numel(objectives)
            %% Objective
            objStr = ['obj ',num2str(i_obj),'/',num2str(numel(objectives))];
            % The objective is a function of a path in the CFG
            objective = objectives(i_obj);
            if matches(objective,"hpDist")
                phi = @(y) min(-d1(y));
            elseif matches(objective,"hpDistExp")
                % Normalize just to feel a little better
                phi = @(y) min(-exp(d1(y))/max(exp(d1)));
            elseif matches(objective,"ordDist")
                phi = @(y) min(-d1ordinary(y));
            elseif matches(objective,"ordDistExp")
                % Normalize just to feel a little better
                phi = @(y) min(-exp(d1ordinary(y))/max(exp(d1ordinary)));
            elseif matches(objective,"depth")
                phi = @(y) min(yData(y));
            elseif matches(objective,"depthExp")
                % Normalize just to feel a little better
                phi = @(y) min(-exp(max(yData)-yData(y))/max(exp(max(yData)-yData)));
            elseif matches(objective,"const")
                phi = @(y) 1;
            else
                error("phi = ?")
            end

            %% Collect statistics on fuzzing
            coverage = nan(numTrials,evalBudget);   % was numEdges_*
            
            %% Fuzz and track coverage
            parfor i = 1:numTrials            
                triStr = ['trial ',num2str(i),'/',num2str(numTrials)];
                disp(['*** Go-Explore: ',...
                    cfgStr,'; ',flagStr,'; ',objStr,'; ',triStr]);

                %% Go-Explore
                warning off;
                history = goExploreFuzzVersion5(initialStates_X{i},...
                    initialStates_Y{i},landmarkIndex{i},gamma,phi,dis_X,...
                    dis_data_Y,update_Y,K,globalGenerator,...
                    evalBudget-counter{i},localGenerator,maxExploreEffort,...
                    flag);
                warning on;

                %% COVERAGE ************
                
                %% Reset, just in case
                data_Y = temp20230508_update(D,[],[]);
                dis_Y_cell = @(y1,y2) dis_data_Y(data_Y,y1,y2);
                
                %% Cell info
                landmarks_Y{i} = initialStates_Y{i}(landmarkIndex{i});
                cellArray{i} = stateCell(dis_Y_cell,landmarks_Y{i},K,initialStates_Y{i});
                
                %% Augment history with initial data
                augmentedHistory = [allInitial_Y{i}(:)',{history.state_Y}];
    
                %% Total coverage (including during initialization)
                edges = [];
                for j = 1:evalBudget
                    src = augmentedHistory{j}(1:(end-1));
                    tar = augmentedHistory{j}(2:end);
                    srctar = unique([src(:),tar(:)],'rows');
                    edges = unique([edges;srctar],'rows');
                    coverage(i,j) = size(edges,1);
                end

            end
            meanCoverage{m,i_flag,i_obj} = mean(coverage/numEdges);
            stdCoverage{m,i_flag,i_obj} = std(coverage/numEdges);
        end
    end

    %% Run the global generator and track coverage
    coverage_gg = nan(numTrials,evalBudget);   % was numEdges_*
    parfor i = 1:numTrials
        triStr = ['trial ',num2str(i),'/',num2str(numTrials)];

        %%
        seed = i*1e6;
        rng(seed);

        %%
        baseline_X = cell(1,evalBudget);
        baseline_Y = cell(1,evalBudget);
        baseline_X(1:numel(allInitial_X{i})) = allInitial_X{i};
        baseline_Y(1:numel(allInitial_Y{i})) = allInitial_Y{i};
        edges = [];
        for j = 1:evalBudget
            %% Track progress
            if mod(j,10) == 0
                baseStr = ['baseline ',num2str(j),'/',num2str(evalBudget)];
                disp(['*** just global generator: ',...
                    cfgStr,'; ',triStr,'; ',baseStr]);
            end

            %% Run
            if j > numel(allInitial_X{i})
                baseline_X{j} = globalGenerator();
                baseline_Y{j} = gamma(baseline_X{j});
            end

            %% Count
            src = baseline_Y{j}(1:(end-1));
            tar = baseline_Y{j}(2:end);
            srctar = unique([src(:),tar(:)],'rows');
            edges = unique([edges;srctar],'rows');
            coverage_gg(i,j) = size(edges,1);
        end
    end
    meanCoverage_gg{m} = mean(coverage_gg/numEdges);
    stdCoverage_gg{m} = std(coverage_gg/numEdges);

    %%
    timePerCFG(m) = toc;

    %% 
    save('20230616workspaceX'); % actually ran save('20230616workspace')
end

%%
disp('timePerCFG = ');
disp(timePerCFG);

%%
fileDir = 'Fuzzing/localGraphics/';

%% Plot average over CFGs
figure;
markers = [".","o","x","+","<",">","^","v","square","diamond",...
    "pentagram","hexagram","|","-","*"];
markers = markers(1:numel(flags));
ind = 1:10*numel(markers):evalBudget;
colors = turbo(numel(objectives));
flagStr = char(64+(1:numel(flags))); % char(64+(1:26)) = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
objStr = char(96+(1:numel(objectives))); % char(96+(1:26)) = 'abcdefghijklmnopqrstuvwxyz'
%
legstr = string.empty;
%
hold on;
for i_flag = 1:numel(flags)
    for i_obj = 1:numel(objectives)
        foo = mean(cell2mat(meanCoverage(:,i_flag,i_obj)));
        plo = plot(1:evalBudget,foo,...
            'Marker',markers(i_flag),'Color',colors(i_obj,:),...
            'MarkerIndices',ind+10*randi(numel(flags)),'LineWidth',1);
        legstr(end+1) = append(flagStr(i_flag),objStr(i_obj)); %#ok<SAGROW> 
    end
end
leg = legend(legstr,'Location','southeast','NumColumns',numel(flags)/2);
xlabel('evals','Interpreter','latex');
ylabel('normalized mean edge coverage','Interpreter','latex');
title('performance of goExploreFuzz configurations','Interpreter','latex');
ylim([0,1]);
box on;
text_X = 10*ones(1,numel(flags));
text_Y = linspace(.975,.6,numel(flags));
foo = ["";repmat("goPdf = ",[7,1]);"schedule = ";"schedule = ";"";""];
text_STR = append(string(flagStr(:)),": ",foo(:),flags(:));
text(text_X,text_Y,text_STR);
text_x = 510*ones(1,numel(objectives));
text_y = text_Y(1:numel(objectives));
text_str = append(string(objStr(:)),": \phi = ",objectives(:));
text(text_x,text_y,text_str);
%
fileDir = 'Fuzzing/localGraphics';
fileName = ['averageOverCFGs',num2str(yyyymmdd(datetime))];
print([fileDir,filesep,fileName],'-dpng');
print([fileDir,filesep,fileName],'-dpdf');

%% Plot individual CFGs
% figure;
% n12 = [2,5];
% ha = tight_subplot(n12(1),n12(2),.01,.01,.05);
% %
% markers = [".","o","x","+","<",">","^","v","square","diamond",...
%     "pentagram","hexagram","|","-","*"];
% markers = markers(1:numel(flags));
% ind = 1:numel(markers):evalBudget;
% colors = turbo(numel(objectives));
% %
% for m = 1:numCFG
%     axes(ha(m)); %#ok<LAXES> 
%     hold on;
%     for i_flag = 1:numel(flags)
%         for i_obj = 1:numel(objectives)
%             plo = plot(1:evalBudget,meanCoverage{m,i_flag,i_obj},...
%                 'Marker',markers(i_flag),'Color',colors(i_obj,:),...
%                 'MarkerIndices',ind+i_flag);
%         end
%     end
%     axis([500,1000,.25,1]);
%     box on;
% end

%% Plot mean coverage at evalBudget 
foo = cellfun(@(x)x(end),meanCoverage,'UniformOutput',false);
finalMeanCoverage = squeeze(mean(cell2mat(foo),1));
foo = cellfun(@(x)x(end),stdCoverage,'UniformOutput',false);
finalStdCoverage = squeeze(mean(cell2mat(foo),1));
[bar,ind] = sort(finalMeanCoverage(:),'descend');
[J,K] = ind2sub(size(finalMeanCoverage),ind);
str = join([string(flagStr(J)'),string(objStr(K)')],"");
% Plot
figure('Position',[0,0,560,210]);
hold on;
for i_flag = 1:numel(flags)
    for i_obj = 1:numel(objectives)
        j = find(matches(str,append(flagStr(i_flag),objStr(i_obj))));
        plot(j,bar(j),'Marker',markers(i_flag),'Color',colors(i_obj,:),...
            'LineWidth',1);
        plot([j,j],bar(j)+finalStdCoverage(ind(j))*[-1,1],...
            'Color',colors(i_obj,:),'LineWidth',1);
    end
end
% Reformat str and set xticks
staggered = str;
suf = "      ";
for j = 1:numel(staggered)
    if mod(j,2)==0
        staggered(j) = append(staggered(j),suf); 
    end
end
xticks(1:numel(staggered));
xticklabels(staggered);
% Finalize plot
xlim([0,numel(staggered)+1]);
grid on;
box on;
ylabel({'normalized mean edge',['coverage at ',num2str(evalBudget),' evals']},...
    'Interpreter','latex');
%
fileDir = 'Fuzzing/localGraphics';
fileName = ['averageOverCFGsAtEvalBudget',num2str(yyyymmdd(datetime))];
print([fileDir,filesep,fileName],'-dpng');
print([fileDir,filesep,fileName],'-dpdf');

%% Plot individual coverage at evalBudget 
figure;%('Position',[0,0,560,210]);
for m = 1:numCFG
    finalMeanCoverage = cellfun(@(x)x(end),squeeze(meanCoverage(m,:,:)));
    finalStdCoverage = cellfun(@(x)x(end),squeeze(stdCoverage(m,:,:)));
    [bar,ind] = sort(finalMeanCoverage(:),'descend');
    [J,K] = ind2sub(size(finalMeanCoverage),ind);
    str = join([string(flagStr(J)'),string(objStr(K)')],"");
    %
    subplot(5,2,m);
    hold on;
    for i_flag = 1:numel(flags)
        for i_obj = 1:numel(objectives)
            j = find(matches(str,append(flagStr(i_flag),objStr(i_obj))));
            plot(j,bar(j),'Marker',markers(i_flag),'Color',colors(i_obj,:),...
                'LineWidth',1);
            plot([j,j],bar(j)+finalStdCoverage(ind(j))*[-1,1],...
                'Color',colors(i_obj,:),'LineWidth',1);
        end
    end
%     % Reformat str and set xticks
%     suf = "      ";
%     for j = 1:numel(str)
%         if mod(j,2)==0
%             str(j) = append(str(j),suf); 
%         end
%     end
    xticks(1:numel(str));
    xticklabels(str);
    xtickangle(0);
    % Finalize plot
    xlim([0,10-1e-2+1]);
    grid on;
    box on;
%     ylabel({'normalized mean edge',['coverage at ',num2str(evalBudget),' evals']},...
%         'Interpreter','latex');
end
%
fileDir = 'Fuzzing/localGraphics';
fileName = ['individualCFGsAtEvalBudget',num2str(yyyymmdd(datetime))];
print([fileDir,filesep,fileName],'-dpng');
print([fileDir,filesep,fileName],'-dpdf');
    
%%
figure('Position',[0,0,420,280]); 
pcolor([[finalMeanCoverage,nan(size(finalMeanCoverage,1),1)];nan(1,size(finalMeanCoverage,2)+1)]); 
daspect([1,1,1]);   % need to do this before colorbar to avoid overlap
view([0,0,-90]);    % for matrix view
xticks((1:numel(text_str))+.5);
yticks((1:numel(text_STR))+.5);
xticklabels(text_str);
yticklabels(text_STR);
title(['normalized mean edge coverage at ',num2str(evalBudget),' evals'],...
    'Interpreter','latex');
cbar = colorbar;
cbar.Position = [.7,cbar.Position(2:4)];
fileDir = 'Fuzzing/localGraphics';
fileName = ['averageOverCFGsAtEvalBudget2D',num2str(yyyymmdd(datetime))];
print([fileDir,filesep,fileName],'-dpng');
print([fileDir,filesep,fileName],'-dpdf');

%%
figure('Position',[0,0,420,280]); 
% Performance relative to default
foo = finalMeanCoverage./(ones(size(finalMeanCoverage,1),1)*finalMeanCoverage(1,:));
pcolor([[foo,nan(size(foo,1),1)];nan(1,size(foo,2)+1)]); 
daspect([1,1,1]);   % need to do this before colorbar to avoid overlap
view([0,0,-90]);    % for matrix view
xticks((1:numel(text_str))+.5);
yticks((1:numel(text_STR))+.5);
xticklabels(text_str);
yticklabels(text_STR);
title({['normalized mean edge coverage at ',num2str(evalBudget),' evals'],...
    'as relative proportion to default (top row)'},...
    'Interpreter','latex');
cbar = colorbar;
cbar.Position = [.7,cbar.Position(2:4)];
