function X = symexp(sigma,n,p)

if nargin<3, p=1; end;

X = exprnd(1/sqrt(2),n,p);

X = sign(rand(size(X))-0.5).*X;

X = sigma*X;