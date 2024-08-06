function YQ = lagrange(X, Y, XQ)
%
%	YQ = lagrange(X, Y, XQ)
%
%	Calcola il polinomio interpolante in forma di Lagrange definito dalle
%	coppie (Xi, Yi) nei punti del vettore XQ
%
%	Input:
%		(X,Y): dati del problema
%		XQ: vettore in cui calcolare il polinomio
%
%	Output:
%		YQ: polinomio interpolante in forma di Lagrange

n = length(X);
if length(Y) ~= n || n <= 0
	error('Dati inconsistenti');
end
if length(unique(X)) ~= n
	error('Ascisse non distinte');
end
YQ = zeros(size(XQ));
for i=1:n
	YQ = YQ + Y(i) * lin(XQ, X, i);
end
return