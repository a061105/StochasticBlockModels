
%synthetic data
NOISE_LEVEL = 1e-1;
LAMBDA = 0.1;
n=15;
k=3;

m = floor(n/3);
z1 = [ones(1,m) zeros(1,2*m)]'; 
z2 = [zeros(1,m) ones(1,m) zeros(1,m) ]'; 
z3=  [zeros(1,2*m) ones(1,m)]';
Z_true = [z1 z2 z3];
Z=[z1 z2 z3]; 
M_true = Z*Z';

%R= M_true + NOISE_LEVEL*randn(n); 
%R=(R+R')/2

% %solve
% [c,Z] = boolLasso(R,LAMBDA); 
% 
% %find top support
% [c2,ind]=sort(c,'descend'); 
% Z2 = Z(:,ind);
% 
% Z2(:,1:3)

%% BASELINES

E = triu(rand(n) < NOISE_LEVEL);
R = mod(M_true + E + E', 2);

% Spectral Clustering
addpath SpectralClustering/files
C = SpectralClustering(R, 3, 2);

% SDP Clustering
addpath /scratch/cluster/ianyen/cvx/
cvx_setup
sum_cluster_sizes = n;
sum_squared_cluster_sizes = sum(M_true(:));
Zround = SDPClustering(R, k, sum_cluster_sizes, sum_squared_cluster_sizes);

