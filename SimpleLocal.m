function [eSL,cond] = SimpleLocal(A,eR,delta)
% [eBest,cond] = SimpleLocalCuts(A,eR,delt,T): Implementation of the
% SimpleLocal Algorithm from:
% "A Simple and Strongly-Local Flow-Based Method for Cut Improvement"
% 
% Nate Veldt, David Gleich, Michael Mahoney 
% Proceedings of The 33rd International Conference on Machine Learning, pp. 1938?1947, 2016
%
% Inputs:
%       A - n x n sparse adjacency matrix of the input graph
%       eR - list of the incides of seed nodes
%       delt - locality parameter
%
% Outputs:
%       eSL - the indices of the set found by SimpleLocal, which minimizes
%             the modified relative conductance score
%       cond - the conductance score for the returned set

n = size(A,1);

if numel(eR) == n
    fprintf('Please give a set of indices of the seed set, not an indicator vector\n')
end

[~,~,~,alpha] = set_stats(A,eR);

cond = alpha;
fprintf('\n\nNew Call to SimpleLocal. Delta = %f \n', delta);
fprintf('Beginning with a cut that has conductance %f \n',full(alpha));
eSL = eR;

% Run the 3-Stage-Flow procedure to find the min-cut/max-flow of the modified
% augmented graph for the given values of alpha and delta
eS = ExactMinCut(A,eR,alpha,delta); 

alph0 = alpha;                   % conductance of R


% If we found a better set, then it has better conductance.
% If there was no better set, the zero vector was returned, and we can't
% find the conductance of the zero vector

if nnz(eS) > 0
    % conductance of S
    eSL = eS;
    [~,~,~,alpha] = set_stats(A,eS);
else
    alpha = alph0;                   % new alpha = old alpha
end

%fprintf('New alpha = %f \n', alpha);
fprintf('After first step, found a cut with conductance = %f \n',alpha);

while alpha < alph0
    eSL = eS;
    cond = full(alpha);
    alph0 = alpha;                   % new conductance := old conductance

    eS = ExactMinCut(A,eR,alpha,delta);

    if nnz(eS) > 0
        [~,~,~,alpha] = set_stats(A,eS);
    else
        alpha = alph0;                   % new alpha = old alpha
    end
    fprintf('Improvement found, new alpha parameter = %f \n',alpha);
end

fprintf('Found a cut with %d vertices and a conductance of %f \n',nnz(eSL),cond);
end