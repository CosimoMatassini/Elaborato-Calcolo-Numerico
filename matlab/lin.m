function L = lin(x, xi, i)
%
%	L = lin(x, xi, i)
%
%	Calcola il polinomio di base di Lagrange in funzione degli argomenti
%	passati
%
%	Input:
%		x: vettore in cui calcolare il polinomio
%		xi: vettore ascisse
%
%	Output:
%		L: polinomio di base di Lagrange

L = ones(size(x));
n = length(xi) - 1;
xii = xi(i);
xi = xi([1:i-1, i+1:n+1]);
for k=1:n
	L = L.*(x - xi(k))/(xii - xi(k));
end
return

