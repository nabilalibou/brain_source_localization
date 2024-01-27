function y=betadist(x,a,b)
% probability density function of the beta distribution with parameters a
% and b

y=(x.^(a-1).*(1-x).^(b-1))/beta(a,b);