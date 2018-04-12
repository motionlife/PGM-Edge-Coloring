function samples = samplegen(A, w, burnin, its)
%The output should be an nxnxk tensor of marginals whose (ith, jth, xth) entry is equal to
%the probability that (i, j) belongs to E is colored with color xij
N = size(A, 1);
K = length(w);
if ~issymmetric(A)
    error('Adjacency matrix A must be symmetric!');
end
%Initalize variable---Not necessary a valid assignment
S = randi(K,N,N).* triu(A);
S = S + S';
for i = 1:burnin
    draw();
end
%Storing the samples
samples = zeros(N,N,its);
for i = 1:its
    samples(:,:,i) = draw();
end
%Calculate marginals
% marginals = zeros(N,N,K);
% for k = 1:K
%     for i=1:N
%         for j=i+1:N
%             if A(i,j)>0
%                 margin = sum(samples(i,j,:)==k)/its;
%                 marginals(i,j,k) = margin;
%                 marginals(j,i,k) = margin;
%             end
%         end
%     end
% end
    function state = draw()
        for n =1:N
            for m =n+1:N
                if S(n,m) > 0
                    %update variable(n,m)
                    prob = zeros(1,K);
                    for xnm=1:K
                        prob(xnm) = potential(n,m,xnm) * potential(m,n,xnm) * exp(w(xnm));
                    end
                    total = sum(prob);
                    if total ~= 0
                        s = find(rand<cumsum(prob/total),1);
                        S(n,m) = s;
                        S(m,n) = s;
                    end
                end
            end
        end
        state = S;
    end
    function pt = potential(n,m,k)  %what if there's leaf structure in the graph???Ignored
        C = S(n,:);
        C(m) = [];
        pt = ~any(C==k);
    end
end


