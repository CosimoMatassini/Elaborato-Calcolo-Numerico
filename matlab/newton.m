function [x, iterazioni] = newton(f, df, x0, tol, imax)
%
%	[x, iterazioni] = newton(f, df, x0, tol, imax)
%
%
%	Input:
%		f - funzione da cui ricavare la radice
%		df - derivata della funzione f
%		x0 - punto di partenza
%		tol - tolleranza richiesta
%		imax - numero di iterazioni max
%
%	Output:
%		x - approssimazione della radice
%		iterazioni - numero delle iterazioni eseguite

if nargin < 3
	error('numero di argomenti insufficiente');
elseif nargin == 3
	tol = 10e-16;
	imax = 1000;
elseif nargin == 4
	imax = 1000;
end

x = x0;

for iterazioni = 1:imax
	fx = f(x);
	fxderivata = df(x);
	if fxderivata == 0
		break;
	end
	x = x - fx/fxderivata;
	if abs(x - x0) <= tol * (1 + abs(x0))
		return;
	else
		x0 = x;
	end
end

warning('soluzione non trovata nelle iterazioni massime disponibili');
return