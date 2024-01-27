function [SOut]=MNE(X,A, lambda);


SOut = A'*inv(A*A'+lambda*eye(length(X)))*X;