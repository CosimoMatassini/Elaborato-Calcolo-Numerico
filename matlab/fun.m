function [f, jacobian] = fun(x)
%
%	[f, jacobian] = fun(x);
%
%	Calcola il gradiente e la matrice Hessiana di f(x)
%
%	Input:
%		x - vettore delle ascisse
%
%	Output:
%		f - gradiente della funzione f(x)
%		jacobian - matrice Hessiana di f(x)

x = x(:);
n = length(x);
Q = 4 * eye(n) + diag(ones(n-1, 1), 1) + diag(ones(n-1, 1), -1);
e = ones(n, 1);
alfa = 2;
beta = -1.1;
grad = @(x) Q * x - alfa * e .* sin(alfa * x) - beta * e .* exp(-x);
Jac = @(x) Q - alfa^2 * diag(e .* cos(alfa * x)) + beta * diag(e .* exp(-x));
f = grad(x);
jacobian = Jac(x);
return