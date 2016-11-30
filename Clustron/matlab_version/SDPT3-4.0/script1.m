
%synthetic data
NOISE_LEVEL = 1e-1;
LAMBDA = 0.1;
n=15;

z1 = [1 1 1 1 0 0 0 0 0 0 0 0 0 0 0]'; 
z2 = [0 0 0 0 1 1 1 1 1 1 1 1 0 0 0]'; 
z3=  [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1]'; 
Z=[z1 z2 z3]; 

M_true = Z*Z';
R= M_true + NOISE_LEVEL*randn(n); 
R=(R+R')/2

%solve
[c,Z] = boolLasso(R,LAMBDA); 

%find top support
[c2,ind]=sort(c,'descend'); 
Z2 = Z(:,ind);

Z2(:,1:3)
