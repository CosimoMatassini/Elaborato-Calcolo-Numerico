function x = chebyshev(n, a, b)
%
%	x = chebyshev(n, a, b)
%
%	calcola le n+1 ascisse di Chebyshev sull'intervallo [a, b]
%
%	Input:
%		n: numero di ascisse che vogliamo calcolare
%		a, b: intervallo in cui vengono calcolate le ascisse di Chebyshev
%
%	Output:
%		x: ascisse di Chebyshev calcolate sull'intervallo [a, b]

if a >= b || n <= 0
	error('Dati errati');
end

x = (a+b)/2 + ((b-a)/2) * cos((2*[n:-1:0] + 1)/((2*(n+1)))*pi);
return