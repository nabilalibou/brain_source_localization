function [Sout] =sissy(X, A, T, lambda, alpha)
rho=1;
Niter=60;


P = sparse(rho*(T.'*T+speye(size(T,2))));
APi=A/P;
L=chol(eye(size(A,1))+APi*A.','lower');
s1=A.'*X;



u =zeros(length(T),1);
v = zeros(length(A),1); 
z = zeros(length(T),1);
s = zeros(length(A),1);
y = zeros(length(A),1);


for  i = 1:Niter
    b = s1+rho*(T.'*(z+u/rho)+y+v/rho);
    s = P\b-APi'*(L'\(L\(APi*b)));
    z = prox_op(T*s-(1/rho)*u,'L1', lambda/rho);
    y = prox_op(s-(1/rho)*v, 'L1', lambda*alpha/rho);
    u = u +rho*(z-T*s);
    v = v + rho*(y-s);
   
end
Sout =s;