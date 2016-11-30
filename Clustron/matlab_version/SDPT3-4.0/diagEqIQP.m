%% Solving the problem:
%%
%%     max  x'Cx 
%%     s.t. x_i^2 = 1, i=1...n
%%
%% by solving SDP that maximizes tr(C'X) and then do max-cut-style randomized rounding.

function [x, obj] = diagEqIQP(C)

n = length(C);

[X, obj] = diagEqSDP(C);

V = chol(X);

r = randn(n,1);
r = r ./ norm(r);

x = sign(V'*r);
obj = x'*C*x;
