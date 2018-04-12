clear;
%%%--------------For test of Problem 1-------------------------
% %% Test Graph is a tree structure
A = [0,0,0,1,1;
     0,0,0,1,0;
     0,0,0,1,0;
     1,1,1,0,0;
     1,0,0,0,0];
%% Test Graph is a non-tree structure
% A = [0 1 1 1;
%      1 0 0 1;
%      1 0 0 1;
%      1 1 1 0];
 
w = [1 5 3 4 2 6];
plot(graph(A));

samples = samplegen(A,w,10000,50000);%generate 50000 test samples from gipps sampling
learned_w = colormle(A,samples);
fprintf('Test w is: [%s]; Learned w is: [%s]\n',sprintf('%.1f, ', w),sprintf('%.2f, ', learned_w));
