function Assignment = maxprod(A, w, its)
N = size(A, 1);
K = length(w);
if ~issymmetric(A)
    error('Adjacency matrix A must be symmetric!');
end
MsgC2V = zeros(N,N,K);%Initialize Cliuqe to variable message table:row->clique;col->var;page->color
MsgV2C = zeros(N,N,K);%Initialize Variable to clique message table:row->var;col->clique;page->color
backtracing = cell(N,N,K);
for t = 1:its
    % Update the message passing from clique c to variable v
    for c = 1:N
        nodes = find(A(c,:));
        degree = numel(nodes);
        if degree > 1 %This is an valid clique marked by c; Filter out tree-like leaves
            for dg = 1:degree
                v = nodes(dg);
                neibrs = nodes(nodes~=v);
                Msg = squeeze(MsgV2C(neibrs,c,:));
                width = degree - 1;
                for xi = 1:K
                    Xc = permn([1:xi-1 xi+1:K],width);% get all the possible value combination
                    max_msg = 0;
                    for r = 1:size(Xc,1)
                        xc = Xc(r,:);
                        if width == numel(unique(xc))
                            sum_msg = sum(Msg(sub2ind([width K],1:width,xc)));
                            if(sum_msg > max_msg)
                                max_msg = sum_msg;
                                if t==its %back-tracing storage
                                    backtracing{c,v,xi} = xc;
                                end
                            end
                        end
                    end
                    MsgC2V(c,v,xi) = max_msg /2;
                end
            end
        end
    end
    % Update the message passing from variable v to clique c
    for v = 1:N
        for c = 1:N
            if (A(v,c) > 0)%This is a valid variable marked by v->c; Filter out non-edges
                MsgV2C(v,c,:) = (w + squeeze(MsgC2V(v,c,:))') /2;
            end
        end
    end
    
end

Assignment = zeros(N,N);
for i = 1:N
    for j = i+1:N
        if (A(i,j)>0)
            maxmaginal = w + squeeze(MsgC2V(i,j,:) + MsgC2V(j,i,:))';
            [~,m_color] = max(maxmaginal);
            if Assignment(i,j) == 0
                Assignment(i,j) = m_color;
                %then we could get the configuration of clique i and clique j based on the m_color and back-tracing storage
                inodes = find(A(i,:));
                ineibrs = inodes(inodes~=j);
                if all(Assignment(i,ineibrs)==0)
                    Assignment(i,ineibrs) = backtracing{i,j,m_color};
                    Assignment(ineibrs,i) = backtracing{i,j,m_color};
                end
            end
            if Assignment(j,i) == 0
                Assignment(j,i) = m_color;
                jnodes = find(A(j,:));
                jneibrs = jnodes(jnodes~=i);
                if all(Assignment(j,jneibrs)==0)
                    Assignment(j,jneibrs) = backtracing{j,i,m_color};
                    Assignment(jneibrs,j) = backtracing{j,i,m_color};
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
