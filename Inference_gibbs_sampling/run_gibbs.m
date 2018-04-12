clear;
%%%--------------For test of part1 of Problem 2-------------------------
A = [0 1 1 1;
    1 0 0 1;
    1 0 0 1;
    1 1 1 0];

w = [1 2 3 4];
marginals = gibbs(A,w,10000,10000);
disp(marginals);
w2 = [-5 -4 -3 -2];
marginals2 = gibbs(A,w2,10000,10000);
disp(marginals2);
%%%--------------For part2 of Problem 2-------------------------
%the probability that (a, d) is colored with color 4
%using weights w=[1,2,3,4] in the graph A.
A = [0 1 1 1;
    1 0 0 1;
    1 0 0 1;
    1 1 1 0];
w = [1 2 3 4];
iterations = [2^6, 2^10, 2^14, 2^18];
len = length(iterations);
result = zeros(len,len);
tic
parfor i=1:len
    [~,samples] = gibbs(A,w,iterations(i),iterations(4)); %#ok<*PFBNS>
    for j=1:len
        result(i,j) = sum(samples(1,4,1:iterations(j))==4)/iterations(j);
    end
end
its = {'its_64','its_1024','its_16384','its_262144'};
burnins = {'burnin_64','burnin_1024','burnin_16384','burnin_262144'};
array2table(result, 'VariableNames', its, 'RowNames', burnins)
toc
