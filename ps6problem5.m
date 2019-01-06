%{
Use the spectral partitioning method to find a division of the graph G into 
two subgraphs of size n1 = 20 and n2 = 30 such that the number of edges 
between the subgraphs (the cut size) is minimized. Visualize the obtained 
partition and find the corresponding cut size. 
%}
clear; clc; close all;

n1 = 20;
n2 = 30;
n = 50;

% construct adjacency matrix
A = ones(50);
A = A - diag(ones(1, 50));
A(21:50, 1:20) = 0;
A(1:20, 21:50) = 0;
A(21, 1) = 1;
A(22, 2) = 1;
A(23, 3) = 1;
A(1, 21) = 1;
A(2, 22) = 1;
A(3, 23) = 1;

G = graph(A); % construct graph

L = laplacian(G); % generate graph Laplacian
% find 2 smallest eigenvalues and corresponding (orthonormal) eigenvectors
[Q, D] = eigs(L, 2, 'smallestabs'); 

alpha1 =(n1 - n2) / sqrt(n);
alpha2 = 2 * sqrt(n1 * n2 / n);
q1 = Q(:,1); % eigenvector corresponding to smallest eigenvalue
q2 = Q(:,2); % eigenvector corresponding to second smallest eigenvalue 

[M1, I1] = maxk(q2, n1); % get indices of largest n1 components
[M2, I2] = maxk(q2, n2); % get indices of largest n2 components
s_plus = zeros(n, 1);    % preallocate
s_minus = zeros(n, 1);   % preallocate

for i=1:n1
    % replace n1 largest components of q2 with 1
    a = I1(i);
    s_plus(a) = 1;
end
s_plus(s_plus~=1) = -1; % replace n2 smallest components of q2 with -1

for i=1:n2
    % replace n2 largest components of q2 with -1
    a = I2(i);
    s_minus(a) = -1;
end
s_minus(s_minus~=-1) = 1; % replace n1 smallest components of q2 with 1

R_plus = 0.25 * s_plus' * L * s_plus;    % calculate cut size
R_minus = 0.25 * s_minus' * L * s_minus; % calculate cut size

nodes = find(s_minus == 1);           % identify partition
nodes2 = find(s_minus == -1);         % identify partition
figure;
h = plot(G, 'NodeLabel', '');         % plot original graph
highlight(h, nodes,'NodeColor','g');  % visualize partition
highlight(h, nodes2,'NodeColor','r'); % visualize partition
highlight(h, 'Edges', [1:193], 'EdgeColor', 'g');
highlight(h, 'Edges', [194:628], 'EdgeColor', 'r');
highlight(h, 'Edges', [20, 39, 57], 'EdgeColor', 'b');
hold on;
w = zeros(3, 1);
w(1) = scatter(NaN,NaN,'og', 'filled');
w(2) = scatter(NaN,NaN,'or', 'filled');
w(3) = plot(NaN,NaN,'b-');
legend(w, 'subgraph of size 20', 'subgraph of size 30', ...
    'edges corresponding to minimum cut', 'Location', 'northwest');

% The smallest cut size is R_minus = 3, and the corresponding partition is
% given by s_minus.