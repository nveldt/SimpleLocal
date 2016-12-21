function [F,eS,flow] = GurobiMaxFlow(A)
% GurobiMaxFlow = calculate the maximum flow matrix F for graph A
%
% If desired this can be replaced with any min-cut/max-flow subroutine
% where the input is an adjacency matrix A where the first node is the
% source node and the last node is the sink node, and the output at least
% gives a minimum cut set eS.
%
% A - adjacency matrix of a graph G
% 
% F - Flow matrix
% eS - minimum cut set

n = size(A,1);
m = nnz(A);
d = sparse(sum(A,2));
p = nnz(A)/2;
s = 1;
t = n;

k = nnz(A(:,n));    %how many nonzeros in the last column of A?

es = zeros(n,1);
es(s) = 1;

[r,c,caps] = find(A);


% Time to construct the incidence matrix based on this ordering
B = sparse(n,m+1);

% create the incidence matrix
for j = 1:m
   B(r(j),j) = -1;
   B(c(j),j) = 1;
end

% the flow is another edge sort of
B(n,m+1) = -1;
B(1,m+1) = 1;

% Constraints
Iz = [speye(m), zeros(m,1)];
G = [B; Iz];

b = [zeros(n,1);caps];

clear model
model.obj = full([zeros(m,1); 1]);

model.A = sparse(G);
model.rhs = b;
model.sense = [repmat('=',1,n)'; repmat('<',1,m)'];

model.vtype = repmat('C',1,m+1);
model.modelsense = 'max';

clear params;

params.cuts = 2;
params.outputflag = 0;
result = gurobi(model, params);

flow = result.objval;
f = result.x(1:m);

F = sparse(r,c,f,n,n);
F = F - F';

eS = result.pi;   % the dual variables of the flow computation will be the
                  % min cut vector
eS = eS(1:n);

% Need to make sure we take the correct side of the cut
if eS(1) == 0
    eS = ones(n,1) - eS;
end

end