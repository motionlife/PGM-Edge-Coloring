% Create a random adjacency matrix with node n
n = 5;
% A = round(rand(n));
A = triu(A) + triu(A,1)';
A = A - diag(diag(A));
%Plot the graph
plot(graph(A));
its = 300;
w = [1 3 1 2];
tic
Z = sumprod(A,w,its);
%X = maxprod(A,w,its)
toc