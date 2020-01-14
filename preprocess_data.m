function region = preprocess_data(region)
[nrows,ncols,~]  = size(region);

time_start  = 1;
time_end    = size(region,3);

g = imageGraph([nrows,ncols],8);

affinity_self       = sum(region(:,:,time_start:time_end).^2,3) + eps;
affinity_up         = sum(region(1:end-1,:,time_start:time_end) .* region(2:end,:,time_start:time_end),3) ...
    ./ (sqrt(affinity_self(1:end-1,:)) .* sqrt(affinity_self(2:end,:)));
affinity_right      = sum(region(:,1:end-1,time_start:time_end) .* region(:,2:end,time_start:time_end),3) ...
    ./ (sqrt(affinity_self(:,1:end-1)) .* sqrt(affinity_self(:,2:end)));
affinity_up_right   = sum(region(1:end-1,2:end,time_start:time_end)  .* region(2:end,1:end-1,time_start:time_end),3) ...
    ./ (sqrt(affinity_self(1:end-1,2:end)) .* sqrt(affinity_self(2:end,1:end-1)));
affinity_down_right = sum(region(2:end,2:end,time_start:time_end)  .* region(1:end-1,1:end-1,time_start:time_end),3) ...
    ./ (sqrt(affinity_self(2:end,2:end)) .* sqrt(affinity_self(1:end-1,1:end-1)));

up_inds         = (g.Edges.EndNodes(:,2)-g.Edges.EndNodes(:,1)==1);
right_inds      = (g.Edges.EndNodes(:,2)-g.Edges.EndNodes(:,1)==nrows);
up_right_inds   = (g.Edges.EndNodes(:,2)-g.Edges.EndNodes(:,1)==(nrows-1));
down_right_inds = (g.Edges.EndNodes(:,2)-g.Edges.EndNodes(:,1)==(nrows+1));

g.Edges.Weight(up_inds)         = affinity_up(:);
g.Edges.Weight(right_inds)      = affinity_right(:);
g.Edges.Weight(up_right_inds)   = affinity_up_right(:);
g.Edges.Weight(down_right_inds) = affinity_down_right(:);

nn = numnodes(g);
[sOut,tOut] = findedge(g);
K = sparse(sOut,tOut,g.Edges.Weight,nn,nn);
K = K + K';
K = spdiags(ones(nrows*ncols,1),0,K);

D = sum(K,2)+eps;
one_over_D = spdiags((1./D),0,size(K,1),size(K,2));
P = one_over_D * K;
%
tmp_flat = reshape(region(:,:,time_start:time_end),[],time_end-time_start+1);
tmp_flat = P*tmp_flat;

tmp_flat = imfilter(tmp_flat,[1 1 1 1]/4,'symmetric');
region = reshape(tmp_flat,nrows,ncols,[]);

return

%%