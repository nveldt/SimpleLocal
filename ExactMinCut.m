function xC = ExactMinCut(A,eR,alpha,delta)
% xC = ExactMinCut(A,eR,alpha,delta): finds the exact minimum cut of the
% modified augmented graph for parameters alpha and delta.
% This is an updated version of the 3-Stage-Flow procedure outlined in
%
% "A Simple and Strongly-Local Flow-Based Method for Cut Improvement"
% 
% Nate Veldt, David Gleich, Michael Mahoney 
% Proceedings of The 33rd International Conference on Machine Learning, pp. 1938?1947, 2016
%
% We do not explicitely calculate flows, but instead compute exact minimum
% cuts on the current "local graph", a subgraph of the modified augmented
% graph, and update the local graph based on which nodes to the sink are
% cut at each iteration.
%
%
% Inputs:
%       A - n x n sparse adjacency matrix of the input graph
%       eR - list of the incides of seed nodes
%       delta - locality parameter
%       alpha - parameter for augmented graph
%
% Outputs:
%       xC - minimum cut set of the modified augmented graph
%
%   Note: if no improved cut is found, we return the zero vector

n = size(A,2);          % n will always be the original number of nodes
volA = sum(nonzeros(A));

% build all the degrees of the nodes
d = sum(A,1); 

if numel(eR) == n
    % this is an indicator
    R = find(eR);
else
    % this is a set of indices
    R = eR;
end
R = unique(R); % remove any duplicates

totalNodesFullyVisited = numel(R);

% s - source node
% t - sink node
% R - seed set
% N(R) - neighbors of the seed set

% Construct the first local graph, with nodes s,t,R,N(R)
%       edges from s to R, t to N(R), and edges in R, and R to N(R)

volR = sum(d(R));
fR = volR/(volA - volR);
W = R;
ARR = A(R,R);
fullyvisited = R;

[WG,inds] = assemble_local_graph(A,R,W,d,alpha,delta,fR,ARR);

[~,xS,~] = GurobiMaxFlow(WG);
S = inds(xS(2:end-1) > 0.5);     % Finding the global indices of node set

% Now we must find the (original) indices of the nodes in WG that need to
% be expanded on.

% If a node is in S (i.e. xS(i) = 1) and it has not been previously
% expanded, then we expand on it in this iteration.

E = setdiff(S,fullyvisited);     % Expand on every node in the cut that has
                                 % not previously been expanded on

if numel(S) > 0
    [~,~,~,bestCond]  = set_stats(A,R,volA);                           
    [~,~,~,cond] = set_stats(A,S,volA);
end

while numel(E) > 0 && numel(S) > 0 %&& bestCond >= alpha  
    
    % we continue if there are nodes to be expanded on
    % and if we haven't already found an improvement cut
    
    totalNodesFullyVisited = totalNodesFullyVisited + numel(E);
    %fprintf(' %d new nodes to expand on \n',numel(E));
    
    % Update which nodes have been fully visited, i.e. expanded on
    fullyvisited = [fullyvisited;E(:)]; % add to fullyvisited
    assert(numel(intersect(W,E)) == 0);
    
    W = [W; E(:)];
    %fprintf('%d nodes we expand on \n',numel(E));
    
    
    % By assempbling a new working graph, we are doing the "expansion"
    [WG,inds] = assemble_local_graph(A,R,W,d,alpha,delta,fR,ARR);
    [~,xS,c] = GurobiMaxFlow(WG);
    S = inds(xS(2:end-1) > 0.5);    % Finding the global inices of node set
    %fprintf('\t the cut has cost %f \n',c);
    
    
    if numel(S) == 0
        % no improvement found
     %   fprintf('numel(S) = 0\n');
        continue
    end
    
    [~,~,~,condS] = set_stats(A,S,volA);
    
    E = setdiff(S,fullyvisited);
    
end
xC = zeros(n,1);
xC(S) = ones(numel(S),1);

%fprintf('\tIn one call to ExactMinCut, total nodes fully visited = %d \n',totalNodesFullyVisited);
end

function [WG,inds] = assemble_local_graph(A,R,W,d,alph,delt,fR,ARR)
    WnR = setdiff(W,R); 
    % so the indices

    % find the two boundary sets
    BR = A(:,R);
    BWnR = A(:,WnR);

    [ei,~] = find(BR);
    Rbound = setdiff(unique(ei),W); % evertyhing coming out of R that isn't in W

    [ei,~] = find(BWnR);
    WnRbound = setdiff(unique(ei),W); % everything ... (same as above)

    B = union(Rbound,WnRbound);

    nR = numel(R);
    nWnR = numel(WnR);
    nB = numel(B);

    % The graph we want to construct has three major regions
    % R - the internal structure of the reference set
    % Wn - the remainder of the working graph
    % B - the boundary of the working graph

    % ARR = A(R,R); % now included as an input
    AWR = A(WnR,R); 
    ABR = A(B,R);
    AW = A(WnR,WnR);
    ABW = BWnR(B,:);

    sR = d(R)*alph;
    sWnR = sparse(1,nWnR);
    sB = sparse(1,nB);

    tR = sparse(1,nR);
    tWnR = d(WnR)*alph*(fR+delt);
    tB = d(B)*alph*(fR+delt);

    WG = [ 0     sR   sWnR sB            0
          sR'   ARR  AWR' ABR'          tR'
          sWnR' AWR  AW   ABW'          tWnR'
          sB'   ABR  ABW  sparse(nB,nB) tB'
          0     tR   tWnR tB            0];

    inds = [R(:); WnR(:); B(:)];
end