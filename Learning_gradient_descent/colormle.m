% function [w,W] = colormle(A, samples)
% %Using gradient descent method to learn weight parameters in the edge
% %coloring problem, using belief propagation to inference probability in
% %each iteration.
% if ~issymmetric(A)
%     error('Adjacency matrix A must be symmetric!');
% end
% N = size(A, 1);
% M = size(samples,3);
% K = max(samples(:)); %get the number of colors from samples
% W = initw(A,N,K); %initialize W matrix
% bp_its = 10; %set number of belief propagation iterations (not necessary to converge)
% gd_its = 1;
% gd_max_its = 500; %set maximun number of iterations for gradient descent terminate
% tolerance = 1.0e-4;
% step = 0.1;
% 
% while 1
%     p = marginals();%get the distribution corresponding to the lateset W
%     %get the expected value of the feature maps under the marginal distribution
%     delta = zeros(N,N,K);
%     for n1 = 1:N
%         for n2 = n1+1:N
%             if A(n1,n2)>0
%                 %compute the full sufficient statistics and then FSS - Expectation of f.m. (See MLAPP.KeviN.M P677)
%                 delta(n1,n2,:) = squeeze(delta(n1,n2,:))' + sum(samples(n1,n2,:)==(1:K),3)/M - squeeze(p(n1,n2,:))';
%                 delta(n2,n1,:) = delta(n1,n2,:);
%             end
%         end
%     end
%     
%     fprintf('Gradient round #%3d, parameter gradient: \x0394(w) = [%s]\n', gd_its, sprintf('%d, ', unifyW(delta,A)-1))
%     if (gd_its > gd_max_its) || (max(abs(delta(:))) < tolerance)
%         break
%     end
%     delta = delta * step;
%     W = W + delta;
%     gd_its = gd_its + 1;
% end
% w = unifyW(W,A);% get a perfect output from each eadge's weight
% %BP function
%     function blfX = marginals()
%         %scaling factors
%         eta1 = 0.05;%eta1 = 0.1(optional)
%         eta2 = 0.1;
%         MsgC2V = ones(N,N,K);%Initialize Cliuqe to variable message table:row->clique;col->var;page->color
%         MsgV2C = ones(N,N,K);%Initialize Variable to clique message table:row->var;col->clique;page->color
%         for t = 1:bp_its
%             % Update the message passing from clique c to variable v
%             for c = 1:N
%                 nodes = find(A(c,:));
%                 degree = numel(nodes);
%                 if degree > 1 %This is an valid clique marked by c; Filter out tree-like leaves
%                     for dg = 1:degree
%                         v = nodes(dg);
%                         Msg = squeeze(MsgV2C(nodes(nodes~=v),c,:));
%                         width = degree - 1;
%                         for xi = 1:K
%                             Xc = permn([1:xi-1 xi+1:K],width);% get all the possible value combination
%                             for r = 1:size(Xc,1)
%                                 if width == numel(unique(Xc(r,:)))
%                                     MsgC2V(c,v,xi) = MsgC2V(c,v,xi) + prod(Msg(sub2ind([width K],1:width,Xc(r,:))));
%                                 end
%                             end
%                         end
%                         MsgC2V(c,v,:) = MsgC2V(c,v,:) * eta1;
%                     end
%                 end
%             end
%             % Update the message passing from variable v to clique c
%             for v = 1:N
%                 for c = 1:N
%                     if A(v,c) > 0 %This is a valid variable marked by v->c; Filter out non-edges
%                         MsgV2C(v,c,:) = exp(squeeze(W(v,c,:))) .* squeeze(MsgC2V(v,c,:)) * eta2;
%                     end
%                 end
%             end
%         end
%         %store the beliefs of variables (a.k.a the color of each edge) key=>edge'i-j'; color domain=> 1:k
%         blfX = zeros(N,N,K);
%         for i = 1:N
%             for j = i+1:N
%                 if A(i,j)>0
%                     blfx = exp(squeeze(W(i,j,:))) .* squeeze(MsgC2V(i,j,:) .* MsgC2V(j,i,:));
%                     sm = sum(blfx);
%                     if sm ~= 0
%                         blfX(i,j,:) = blfx / sm;% normalize variable beliefs
%                     end
%                 end
%             end
%         end
%     end
% end
% function W = initw(A,N,K)
% W = zeros(N,N,K);
% w = rand(1,K);
% for i = 1:N
%     for j = i+1:N
%         if A(i,j)>0
%             W(i,j,:) = w;
%             W(j,i,:) = w;
%         end
%     end
% end
% end
% function w = unifyW(W,A)
% K = size(W,3);
% w = zeros(1,K);
% for k=1:K
%     w(k) = sum(sum(W(:,:,k)));
% end
% w = w/sum(A(:)~=0);
% w = w - min(w) + 1;%align W before output it, to make it start with a minimum value of 1, not required
% end
% function M = permn(V, N)
% nV = numel(V);
% if nV==0 || N == 0
%     M = zeros(nV,N);
% elseif N == 1
%     M = V(:);
% else
%     [Y{N:-1:1}] = ndgrid(1:nV);
%     M = V(reshape(cat(N+1,Y{:}),[],N));
% end
% end


%% Below is version 1
function w = colormle(A, samples)
%Using gradient descent method to learn weight parameters in the edge
%coloring problem, using belief propagation to inference probability in
%each iteration.
if ~issymmetric(A)
    error('Adjacency matrix A must be symmetric!');
end
N = size(A, 1);
M = size(samples,3);
K = max(samples(:)); %get the number of colors from samples
w = rand(1,K); %initialize W vector
bp_its = 17; %set number of belief propagation iterations (not necessary to converge)
gd_its = 1;
gd_max_its = 2000; %set maximun number of iterations for gradient descent terminate
tolerance = 1.0e-6;
step = 0.1;

while 1
    p = marginals();%get the distribution corresponding to the lateset w
    %get the expected value of the feature maps under the marginal distribution
    delta = zeros(1, K);
    for n1 = 1:N
        for n2 = n1+1:N
            if A(n1,n2)>0
                %compute the full sufficient statistics and then FSS - Expectation of f.m. (See MLAPP.KeviN.M P677)
                delta = delta + sum(samples(n1,n2,:)==(1:K),3)/M - squeeze(p(n1,n2,:))';
            end
        end
    end
    fprintf('Gradient round #%3d, parameter gradient: \x0394(w) = [%s]\n', gd_its, sprintf('%d, ', delta))
    if (gd_its > gd_max_its) || (max(abs(delta)) < tolerance)
        break
    end
    delta = delta * step;
    w = w + delta;
    gd_its = gd_its + 1;
end
w = w - min(w) + 1;%align w before output it, to make it start with a minimum value of 1, not required
%BP function
    function blfX = marginals()
        %scaling factors
        eta1 = 0.05;%eta1 = 0.1(optional)
        eta2 = 0.1;
        MsgC2V = ones(N,N,K);%Initialize Cliuqe to variable message table:row->clique;col->var;page->color
        MsgV2C = ones(N,N,K);%Initialize Variable to clique message table:row->var;col->clique;page->color
        for t = 1:bp_its
            % Update the message passing from clique c to variable v
            for c = 1:N
                nodes = find(A(c,:));
                degree = numel(nodes);
                if degree > 1 %This is an valid clique marked by c; Filter out tree-like leaves
                    for dg = 1:degree
                        v = nodes(dg);
                        Msg = squeeze(MsgV2C(nodes(nodes~=v),c,:));
                        width = degree - 1;
                        for xi = 1:K
                            Xc = permn([1:xi-1 xi+1:K],width);% get all the possible value combination
                            for r = 1:size(Xc,1)
                                if width == numel(unique(Xc(r,:)))
                                    MsgC2V(c,v,xi) = MsgC2V(c,v,xi) + prod(Msg(sub2ind([width K],1:width,Xc(r,:))));
                                end
                            end
                        end
                        MsgC2V(c,v,:) = MsgC2V(c,v,:) * eta1;
                    end
                end
            end
            % Update the message passing from variable v to clique c
            for v = 1:N
                for c = 1:N
                    if A(v,c) > 0 %This is a valid variable marked by v->c; Filter out non-edges
                        MsgV2C(v,c,:) = exp(w) .* squeeze(MsgC2V(v,c,:))' * eta2;
                    end
                end
            end
        end
        %store the beliefs of variables (a.k.a the color of each edge) key=>edge'i-j'; color domain=> 1:k
        blfX = zeros(N,N,K);
        for i = 1:N
            for j = i+1:N
                if A(i,j)>0
                    blfx = exp(w) .* squeeze(MsgC2V(i,j,:) .* MsgC2V(j,i,:))';
                    sm = sum(blfx);
                    if sm ~= 0
                        blfX(i,j,:) = blfx / sm;% normalize variable beliefs
                    end
                end
            end
        end
    end
end

function M = permn(V, N)
nV = numel(V);
if N > nV
    error("N can't be greater than the length of V");
end
if nV==0 || N == 0
    M = zeros(nV,N);
elseif N == 1
    M = V(:);
else
    [Y{N:-1:1}] = ndgrid(1:nV);
    M = V(reshape(cat(N+1,Y{:}),[],N));
end
end
