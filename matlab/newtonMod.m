function [x, iterazioni] = newtonMod(f, m, df, x0, tol, imax)
%
%	[x, iterazioni] = newton(f, m, df, x0, tol, imax)
%
%
%	Input:
%		f - funzione da cui ricavare la radice
%		m - molteplicità della radice
%		df - derivata della funzione f
%		x0 - punto di partenza
%		tol - tolleranza richiesta
%		imax - numero di iterazioni max
%
%	Output:
%		x - approssimazione della radice
%		iterazioni - numero delle iterazioni eseguite
%

if nargin < 4
	error('numero di argomenti insufficiente');
elseif nargin == 4
	tol = 10e-16;
	imax = 1000;
elseif nargin == 5
	imax = 1000;
end

if m < 1
	error('molteplicità della radice non corretta')
end

x = x0;

for iterazioni = 1:imax
	fx = f(x);
	fxderivata = df(x);
	if fxderivata == 0
		error('derivata uguale a zero');
	end
	x = x - m*fx/fxderivata;
	if abs(x - x0) <= tol * (1 + abs(x0))
		return;
	else
		x0 = x;
	end
end

warning('soluzione non trovata nelle iterazioni massime disponibili');
return