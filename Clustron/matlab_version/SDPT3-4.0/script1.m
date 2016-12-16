
%synthetic data
NOISE_LEVEL = 1e-1;
LAMBDA = 0.1;
n=15;

m = floor(n/3);
z1 = [ones(1,m) zeros(1,2*m)]'; 
z2 = [zeros(1,m) ones(1,m) zeros(1,m) ]'; 
z3=  [zeros(1,2*m) ones(1,m)]';
Z_true = [z1 z2 z3];

M_true = Z_true*Z_true';
%R= M_true + NOISE_LEVEL*randn(n); 
%R=(R+R')/2
E = triu(rand(n) < NOISE_LEVEL);
R = mod(M_true + E + E', 2);
R(R==1)=0.9;
R(R==0)=0.1;
%solve
[c,Z] = boolLasso(R,LAMBDA, M_true); 

%find top support
[c2,ind]=sort(c,'descend'); 
Z2 = Z(:,ind);

Z2(:,1:3)
