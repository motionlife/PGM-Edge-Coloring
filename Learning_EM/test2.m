clear;
%%%--------------For test of Problem 1-------------------------
%% Test Graph is a tree structure
A = [0,0,0,1,1;
     0,0,0,1,0;
     0,0,0,1,0;
     1,1,1,0,0;
     1,0,0,0,0];
 L = false(size(A));
 L(1,4) = 1; % set edge(1,4) as hidden variable
%  L(2,4) = 1;
%% Test Graph is a non-tree structure
% A = [0 1 1 1;
%      1 0 0 1;
%      1 0 0 1;
%      1 1 1 0];
% L = false(size(A));
% L(1,4) = 1;% set edge(1,4) as hidden variable
% L(1,2) = 1;% set edge(1,2) as hidden variable

w = [4 1 3 2];
plot(graph(A));
M = 10000;
samples = samplegen(A,w,7000,M);%generate 50000 test samples using gipps sampling
[r, c] = find(L); samples(r,c,:) = 0;samples(c,r,:) = 0;%mask out the values of hidden variable

learned_w = colorem(A, L, samples);
fprintf('Test w is: [%s]; Learned w is: [%s]\n',sprintf('%.1f, ', w),sprintf('%.2f, ', learned_w));
