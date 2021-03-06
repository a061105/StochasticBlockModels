%% solving the following problem:
%%
%%     min \frac{1}{2}\|R-M\|^2_F  + lambda*\|M\|_S
%% 
%% where \|.\|_S is atomic norm based on the atomic set:
%%
%%     S = { zz' | z_i\in {0,1} }.
%%
%% Our solver reterns c and Z=[z_1,...,z_k] s.t.  M=\sum_j c_j* z_jz_j'.

function  [c, Z] = boolLasso(R, lambda, M_true)

TOL = 1e-2;
T = 100;
T2 = 100;

n = length(R);
Z = [];
c = [];

err_ref = sum(sum(abs(M_true-R)));
M = zeros(n);
for t = 1:T
	%compute gradient
	grad_M = gradient_M_logi( R, M );
	
	%find greedy direction & add to active set
	[z, obj] = boolIQP(-grad_M);
	Z = [Z z];
	c = [c;0];
	%Z = [1 1;1 1;1 1;1 1;1 1;1 0;1 0];
	%c = [1;1];
	
	%fully corrective by prox-GD
	k = length(c);
	h = diag_hessian(Z);
	eta = 2.5e-1/(max(h)*k); %step size
	for t2 = 1:T2*k
		grad_c = gradient_c_logi(R,Z,c);
		c = prox( c - eta*grad_c, eta*lambda );
	end
	
	%shrink c and Z for j:cj=0
	Z = Z(:,c'>TOL);
	c = c(c>TOL);
	
	%dump info
	M = compute_M(c,Z);
	obj = boolLassoObj(R,lambda,M,c);
	%err = sum(sum(abs(M-R)));
	%M_avg = mean(mean(M));
	%l = -sum(sum(  R.*log( M*p+(1-M)*q ) +  (1-R).*log( M*(1-p)+(1-M)*(1-q) ) ));
	%['loss=' num2str(l) ', reg=' num2str(lambda*sum(c)) ', err(M,R)=' num2str(err) ', avg_M=' num2str(M_avg)]
	['t=' num2str(t) ', obj=' num2str(obj)]
end

end

function M = compute_M(c,Z)
	
	[n,k] = size(Z);
	M = zeros(n,n);
	for j = 1:k
		M = M + c(j)*Z(:,j)*Z(:,j)';
	end
end

function grad_M = gradient_M_LS(R,M)
	
	grad_M = -(R-M);
end


function grad_M = gradient_M_logi(R,M)
	p = 0.9;
	q = 0.1;
	
	temp1 = (p-q)*M+q;
	temp2 = (1-q)-(p-q)*M;
	
	grad_M = -(p-q)*R ./ temp1 + (p-q)*(1-R) ./ temp2;
end

function grad_c = gradient_c_LS(R,Z,c)
	
	k = length(c);
	grad_c = zeros(k,1);
	%for j=1:k
	%	ZTzj = Z'*Z(:,j); %k by 1
	%	grad_c(j) = -Z(:,j)'*R*Z(:,j) + 0.5*c'*(ZTzj.^2) + 0.5*c(j)*ZTzj(j).^2;
	%end
	grad_M = gradient_M(R,compute_M(c,Z));
	for j=1:k
		grad_c(j) = Z(:,j)'*grad_M*Z(:,j);
	end
end

function grad_c = gradient_c_logi(R,Z,c)
	
	k = length(c);
	grad_c = zeros(k,1);
	%for j=1:k
	%	ZTzj = Z'*Z(:,j); %k by 1
	%	grad_c(j) = -Z(:,j)'*R*Z(:,j) + 0.5*c'*(ZTzj.^2) + 0.5*c(j)*ZTzj(j).^2;
	%end
	grad_M = gradient_M_logi(R,compute_M(c,Z));
	for j=1:k
		grad_c(j) = Z(:,j)'*grad_M*Z(:,j);
	end
end
function h = diag_hessian(Z)
	k = size(Z,2);
	h = zeros(k,1);
	for i = 1:k
		h(i) = (Z(:,i)'*Z(:,i)).^2;
	end
end

function c2 = prox( c, lambda )
	
	c2 = c;
	c2(c<=lambda)=0;
	c2(c>lambda) = c(c>lambda)-lambda;
end
