function X = maxprod(A, w, its)
%function [Z, blfX, blfC,MsgC2V,MsgV2C] = sumprod(A, w, its)
%The function should return the approximate partition function attained after its
%iterations by plugging the appropriately normalized beliefs into the Bethe free energy.
% A: adjacency matrix of the graph G; w: vector representing the weights, its: the number of iterations
N = size(A, 1);
K = length(w);
if ~issymmetric(A)
    error('Adjacency matrix A must be symmetric!');
end

%scaling factors
eta1 = 0.1;
eta2 = 0.1;
MsgC2V = ones(N,N,K);%Initialize Cliuqe to variable message table:row->clique;col->var;page->color
MsgV2C = ones(N,N,K);%Initialize Variable to clique message table:row->var;col->clique;page->color
for t = 1:its
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
            if (A(v,c) > 0)%This is a valid variable marked by v->c; Filter out non-edges
                MsgV2C(v,c,:) = exp(w) .* squeeze(MsgC2V(v,c,:))' * eta2;
            end
        end
    end
    
end

%store the beliefs of variables (a.k.a the color of each edge) key=>edge'i-j'; color domain=> 1:k
blfX = zeros(N,N,K);
for i = 1:N
    for j = i+1:N
        if (A(i,j)>0)
            blfx = exp(w) .* squeeze(MsgC2V(i,j,:) .* MsgC2V(j,i,:))';
            sm = sum(blfx);
            if sm ~= 0
                blfX(i,j,:) = blfx / sm;% normalize variable beliefs
            end
        end
    end
end

%store the beliefs of cliques (a.k.a the assignments of color for the incidented edges of node i);
blfC = cell(1,N);
for c = 1:N
    V = find(A(c,:));
    num = numel(V);
    if num > 1
        blfc = zeros(ones(1,num)*K);
        Msg = squeeze(MsgV2C(V,c,:));
        XC = permn(1:K,num);% get all the possible value combination
        for r = 1:size(XC,1)
            Xr = XC(r,:);
            if num == numel(unique(Xr))
                temp = num2cell(Xr);
                ind1 = sub2ind(size(blfc),temp{:});
                ind2 = sub2ind([num K],1:num,Xr);
                blfc(ind1) = prod(Msg(ind2));
            end
        end
        sm = sum(blfc(:));
        if(sm ~= 0)
            blfc = blfc / sm;
        end
        blfC{c} = blfc;
    end
end

X = A;

for i =1:N
    nodes = find(A(i,:));
    num = numel(nodes);
    if num > 1
        blf_C = blfC{i};
        [~,position] = max(blf_C(:));
        clen = length(size(blf_C));
        indices = cell(1,clen);
        [indices{:}] = ind2sub(size(blf_C),position);
        for n=1:clen
            X(i,nodes(n)) = indices{n};
            X(nodes(n),i) = indices{n};
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