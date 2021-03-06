function w = colorem(A, L, samples)
%Using EM algorithm learn w from data with missing variables;
if ~issymmetric(A)
    error('Adjacency matrix A must be symmetric!');
end
N = size(A, 1);
M = size(samples,3);
K = max(samples(:)); %get the number of colors from samples
missed = sum(L(:));
if missed == 0
    error('At least set one latent variable!');
end
range = K ^ missed;
likely = zeros(1,range);
fullData = zeros(N,N, M*range); %extend to full sample space considering all possible situation
dataWeights = zeros(M*range,1);%storage the weight for each data sample from full data space
colors = permn(1:K,missed);%fliplr(permn(1:K,missed));
w = rand(1,K); %initialize W vector
bp_its = 9; %set number of belief propagation iterations (not necessary to converge)
gd_max_its = 500; %set maximun number of iterations for gradient descent terminate
em_max_its = 200;
tolerance = 1.0e-7;
step = 0.1;
em_its = 1;
while 1
    w_old = w;
    %% 1. E Step: get expected sufficient statistic (ESS)
    for m = 1:M
        smpl = samples(:,:,m);
        mark = (m - 1) * range;
        for idx = 1: range
            Xm = colors(idx,:);
            smpl(L)  = Xm;
            smpl(L') = Xm;
            likely(idx) = inference(smpl,N,w,Xm);
            fullData(:,:,mark + idx) = smpl;
        end
        dataWeights(mark + (1 : range)) = likely / sum(likely(:));% normalization for each group
    end
    %% 2. M Step: Gradient ascent get the most likely parameter based on the ESS from above;
    gd_its = 1;
    while 1
        p = marginals();%get the distribution corresponding to the lateset w
        %get the expected value of the feature maps under the marginal distribution
        delta = zeros(K, 1);
        for n1 = 1:N
            for n2 = n1+1:N
                if A(n1,n2)>0
                    %compute the expected sufficient statistics for variable(n1,n2)
                    delta = delta + squeeze(fullData(n1,n2,:)==1:K) * dataWeights/M - squeeze(p(n1,n2,:));% ESS - Expectation of f.m.
                end
            end
        end
        %fprintf('Gradient round #%3d, parameter gradient: \x0394(w) = [%s]\n', gd_its, sprintf('%d, ', delta))
        if (gd_its > gd_max_its) || (max(abs(delta)) < tolerance)
            break
        end
        delta = delta * step;
        w = w + delta';
        gd_its = gd_its + 1;
    end
    
    %% 3. check if reaches a local minima
    changes = abs(w - w_old);
    fprintf('EM round #%3d, w change: \x0394w = [%s]\n', em_its, sprintf('%d, ', changes))
    if (max(changes) < 1.0e-9) || (em_its > em_max_its)
        break
    end
    em_its = em_its + 1;
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
%Inference functon used in E step
function prob = inference(data,N,w,Xm)
prob = 0;
for c = 1:N
    clq = data(c,:);
    Xc = clq(clq~=0);
    num = numel(Xc);
    if num > 1 % Then node c is a potentially valid 'clique'
        if(num ~= length(unique(Xc)))
            return
        end
    end
end
prob = prod(exp(w(Xm)));
end

function M = permn(V, N)
nV = numel(V);
if nV==0 || N == 0
    M = zeros(nV,N);
elseif N == 1
    M = V(:);
else
    [Y{N:-1:1}] = ndgrid(1:nV);
    M = V(reshape(cat(N+1,Y{:}),[],N));
end
end

