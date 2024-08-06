function dP = hornerDerivata(ascisse, coefficienti, xi)
%
%	dP = hornerDerivata(ascisse, coefficienti, xi)
%
%	Calcola la derivata di un polinomio in forma di Newton in un punto specifico
%
%	Input:
%		ascisse - Vettore di ascisse [x0, x1, ..., xn]
%		coefficienti - Vettore dei coefficienti [a0, a1, ..., an]
%		xi - Ascissa su cui valutare la derivata
%	Output:
%		dP - Valore della derivata del polinomio in xi
if nargin < 3
    error("Numero di parametri insufficienti");
end
n = length(coefficienti);
if n ~= length(ascisse)
    error("Dimensione degli input errata");
end
P = coefficienti(n);
dP = 0;
for k = n-1:-1:1
	dP = dP .* (xi - ascisse(k)) + P;
	P = P .* (xi - ascisse(k)) + coefficienti(k);
end
return