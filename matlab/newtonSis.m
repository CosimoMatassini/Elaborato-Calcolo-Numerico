function [x, nit] = newtonSis(fun, x0, tol, maxit)
%
%	[x, nit] = newtonSis(fun, x0, tol, maxit)
%
%	Metodo di newton per la risoluzione di sistemi di equazioni non lineari
%
%	Input:
%		fun: [f, jacobian] = fun(x) se il sistema da risolvere e f(x)=0
%			f: gradiente di una funzione f(x) di cui vogliamo approssimare una radice
%			jacobian: matrice Hessiana di f(x);
%		x0: vettore valori iniziali
%		tol: tolleranza
%		maxit: numero massimo di iterazioni
%
%	Output:
%		x: soluzione del sistema
%		nit: numero di iterazioni eseguite

if nargin < 2
	error('Numero di argomenti insufficiente');
elseif nargin == 2
	tol = 1e-6;
	maxit = 1000;
elseif nargin == 3
	maxit = 1000;
elseif maxit <= 0 || tol <= 0
	error('Dati in ingresso errati');
end

x0 = x0(:);
nit = maxit;
for i = 1:maxit
	[f, jacobian] = fun(x0);
	delta = mialum(jacobian, -f);
	x = x0 + delta;
	if norm(delta ./ (1 + abs(x0)), inf) <= tol
		nit = i;
		break
	end
	x0 = x;
end
return

function x = mialum(A, b)
%	x = mialum(A, b);
%
%	Risolve il sistema lineare Ax = b con fattorizzazione LU senza pivoting parziale.
%
% 	Input:
%		A - matrice dei coefficienti
%		b - vettore dei termini noti
%
%	Output:
%		x - vettore soluzione

[m, n] = size(A);
if m ~= n
	error('La matrice non e quadrata!');
end
[m, o] = size(b);
if m ~= n || o ~= 1
	error('Dimensione del vettore e della matrice non compatibili')
end

for i = 1:n-1
	if A(i, i) == 0
		error('Matrice singolare');
	end
	for j = i+1:n
		A(j, i) = A(j, i) / A(i, i);
		A(j, i+1:n) = A(j, i+1:n) - A(j, i) * A(i, i+1:n);
	end
end

x = b(:);
for i = 2:n
	x(i:n) = x(i:n) - A(i:n, i-1) * x(i-1);
end

for i = n:-1:1
	x(i) = x(i) / A(i, i);
	x(1:i-1) = x(1:i-1) - A(1:i-1, i) * x(i);
end
return
