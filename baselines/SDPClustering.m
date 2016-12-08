% Returns n x k cluster indicator matrix Zround.
% Uses the SDP for general cluster structure in Hajek Wu Xu '16

function Zround = SDPClustering(R, k, sum_cluster_sizes, sum_squared_cluster_sizes)

n = size(R, 1);
cvx_begin
    variable X(n,n)
    maximize( trace(X' * R) )
    subject to
        X == semidefinite(n)
        X >= 0
        X <= 1
        trace(eye(n) * X) == sum_cluster_sizes
        trace(ones(n) * X) == sum_squared_cluster_sizes
cvx_end

[Z, c] = eigs(X, k);
Zround = round(abs(Z * sqrt(c)));


end