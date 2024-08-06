function [x, iterazioni] = secanti(f, x0, x1, tol, itmax)
%
%	[x, iterazioni] = secanti(f, xo, x1, tol, itmax)
%
%	Calcola una approssimazione della radice di f(x)
%	con tolleranza tol
%
%	Input:
%		f - funzione da cui ricavare la radice;
%		x0, x1 - punti iniziali;
%		tol - accuratezza richiesta
%		itmax - numero massimo di iterazioni
%
%	Output:
%		x - approssimazione della soluzione
%		iterazioni - numero delle iterazioni eseguite

if nargin < 3
	error('numero di argomenti in ingresso errato')
elseif nargin == 3
	tol = 10e-16;
	itmax = 1000;
elseif nargin == 4
	itmax = 1000;
end

if tol <= 0
	error('tolleranza errata')
end
if itmax <= 0
	error('numero di iterazioni massimo errato')
end

f0 = f(x0);
f1 = f(x1);
for iterazioni = 1:itmax
	if f1 == f0
		error('impossibile eseguire il metodo');
	end
	x = (f1 * x0 - f0 * x1) / (f1 - f0);
	if abs(x - x1) <= tol
		break
	elseif iterazioni < itmax
	x0 = x1;
	f0 = f1;
	x1 = x;
	f1 = f(x1);
	end
end
if abs(x - x1) > tol
	error('soluzione non trovata nelle iterazioni massime disponibili');
end
return
