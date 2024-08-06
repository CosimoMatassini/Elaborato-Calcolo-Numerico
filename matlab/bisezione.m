function [x, iterazioni] = bisezione(f, a, b, tol, itmax)
%
%	[x, iterazioni] = bisezione(f, a, b, tol, itmax)
%
%	calcola approssimazione della radice di
%	f(x) con tolleranza tol
%
%	Input:
%		f - funzione da cui calcolare la radice
%		a, b - punti iniziali
%		tol - tolleranza richiesta
%		itmax - numero iterazioni max
%
%	Output:
%		x - approssimazione della soluzione
%		iterazioni - numero delle iterazioni eseguite

if nargin < 3
	error('numero di argomenti in ingresso errato');
elseif nargin == 3
	tol = 1e-6;
	itmax = ceil(log2(b - a) - log2(tol));
elseif nargin == 4
	itmax = ceil(log2(b - a) - log2(tol));
end

if(b <= a)
	error('intervallo iniziale errato');
end

if tol <= 0
	error('tolleranza non valida');
end

iterazioni = 0;

fa = f(a);
fb = f(b);

if fa == 0
	x = a;
	return;
end
if fb == 0
	x = b;
	return;
end

if fa * fb > 0
	error('intervallo non accettabile');
end

for iterazioni = 1:itmax
	x = (a + b) / 2;
	fx = f(x);
	fxderivata = abs(fb - fa) / (b - a);

	if abs(fx) / abs(fxderivata) <= tol
		break;
	elseif fa * fx > 0
		a = x;
		fa = fx;
	else
		b = x;
		fb = fx;
	end
end
return