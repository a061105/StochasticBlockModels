function obj = boolLassoObj(R,lambda,M,c)
	
	obj = loss_logi(R,M) + lambda*sum(c);
end

function l = loss_sq(R,M)

	D = R-M;
 	l = sum(sum(D.*D))/2 
end

function l = loss_logi(R,M)

p = 0.9;
q = 0.1;

l = -sum(sum(  R.*log( M*p+(1-M)*q ) +  (1-R).*log( M*(1-p)+(1-M)*(1-q) ) ));

end
