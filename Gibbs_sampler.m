function [SOut,LambdaOut]=Gibbs_sampler(X,A, sigma_s2, sigma_n2, lambda);


%number of iterations
Niter=200;

%get number of sensors and number of dipoles
[N,D]=size(A);

% constants 
nA = sum(A.^2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


S = zeros(D,1);
%si=0;
Q = zeros(D,Niter);

% Iterate and converge
for j=1:Niter
    e = X - A*S; 
    L = 0;
    for i=1:D
        ei=e + A(:,i) * S(i,:);
        sigma_i2 = sigma_n2*sigma_s2/(sigma_n2+sigma_s2*nA(i));
        mui = (sigma_i2/sigma_n2)*A(:,i)'*ei;
        nui=lambda*sqrt((sigma_i2/sigma_s2))*exp(mui^2/(2*sigma_i2));
        lambda_i = nui/(nui+1-lambda);
        

        %bernoulli
        test = rand(1);
        if test<lambda_i
            Q(i,j)=1;
        end


        if Q(i,j)==1
            S(i) = mui + sqrt(sigma_i2) * randn(1);
            L = L+1;
        else
            S(i) = 0;
        end
        e = ei - A(:,i)*S(i);
        
    end

end

 
%   SOut = S;
%   LambdaOut=0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make the estimation from Q and S using the MAP criterium
% Q(:,Niter/2:Niter)

for j = 1:Niter/2
    q = Q(:,Niter/2+j); 
    idx = find(q==1); 
    R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
    S_opt = (R\(A(:,idx)'*X))/sigma_n2; 
    cout(j) = - S_opt'*R*S_opt/sigma_n2; 
end

[Val, OptIdx] = min(cout); 
q = Q(:,Niter/2+OptIdx);
idx = find(q==1); 
R = (A(:,idx)'*A(:,idx))/sigma_n2 + eye(length(idx))/sigma_s2; 
SOut = zeros(D,1); 
SOut(idx) = (R\(A(:,idx)'*X))/sigma_n2;
LambdaOut = length(idx)/D; 
