%% Demonstrate the performance of SimpleLocal in finding the target ventricle
% from a seed set

% Illustrate the process quickly on a small subgraph of the 18 million node
% Brain graph, which contains the target ventricle
%
clear
load BrainSubgraph

fprintf('You have loaded a graph with %d nodes, volume = %f \n',size(SmallBrain,1),full(sum(sum(SmallBrain))/2));

%% Display information for the graph and the target set
size(SmallBrain)
Rnew = Ventricle;
[~,~,~,condRset] = set_stats(SmallBrain,Rnew)

SmallBrainUnweighted = spones(SmallBrain);


%% Get seed set
nseeds = 75;
rng(1); % reset the random numbers
p = randperm(numel(Rnew));
Rseedrand = Rnew(p(1:nseeds));

% Grow the randomly chosen nodes by their neighborhood, to get seed set R
Rneighbs = SmallBrain(Rseedrand,:);
[~,R] = find(Rneighbs);
R = unique([R; Rseedrand]);
volA = sum(nonzeros(SmallBrain));

% Display information for the seed set
[cut,vol,edges,cond] = set_stats(SmallBrain,R,volA);

fprintf('We will run SimpleLocal on a seed set with the following stats: \n')
fprintf('cut(R) = %d\nvol(R) = %f \nConductance(R) = %f \n',cut,vol,cond);


%% Run SimpleLocal on the seed set to try to recover the target set

%for delta = .1, will take a couple minutes and return a very good conductance set
delta = .1; 

tic
[SLcut,SLcond] = SimpleLocal(SmallBrain, R, delta);
timeSL = toc;
fprintf('SimpleLocal took %f seconds when delta = %f \n',timeSL,delta);


