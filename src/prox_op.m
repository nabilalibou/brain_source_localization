function x=prox_op(y,type,lambda)
%proximal operator
%
%prox_op(y,type,lambda,mu)
%
%computes the proximal operator used in FISTA for L1-norm regularization 
%and L12-norm regularization terms.
%
%INPUT: y - input vector for the proximal operator
%       type - specifies the type of the regularization term:
%               'L1' - L1-norm regularization
%               'L12' - L12-norm (mixed norm) regularization: L1-norm over
%                       first dimension, L2-norm over second dimension
%       lambda - regularization parameter
%
%OUTPUT: x - solution to the proximal problem
%
% Hanna Becker, October 2014

switch type
    case 'L1'
        %weight vector
        w=ones(size(y));
        
        %determine elements that are approximately 0
        idx=find(abs(y)>1e-15);
        
        %initialization of x
        x=zeros(size(y));
        
        %update non-zero elements of x
        x(idx)=(y(idx)./abs(y(idx))).*max(abs(y(idx))-lambda*w(idx),zeros(size(y(idx))));
    
    case 'L12'
        [N,K]=size(y);
        
        %weight vector
        w=ones(N,1);
        
        %determine rows whose L2-norm is approximately 0
        ny=normm(y);
        idx=find(ny>1e-15);
        
        %initialization of x
        x=zeros(size(y));
        
        %update non-zero elements of x
        x(idx,:)=y(idx,:).*(max(ones(length(idx),1)-lambda*sqrt(w(idx))./ny(idx),zeros(length(idx),1))*ones(1,K));

end