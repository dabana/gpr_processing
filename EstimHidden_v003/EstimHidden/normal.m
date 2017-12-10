function X = normal(sigma,n,p)

if nargin<3, p=1; end;

X = sigma*randn(n,p);