function [x, nr] = miaqr(A, b)
%
%	[x, nr] = miaqr(A,b)
%
%	Esegue la fattorizzazione QR di A
%	restituendo la soluzione ai minimi quadrati del sistema
%	e la norma del corrispondente vettore residuo
%
%	Input:
%		A: matrice m x n
%		b: vettore dei termini noti
%
%	Output:
%		x: soluzione del sistema Ax = b
%		nr: norma vettore residuo

[m, n] = size(A);
[k, o] = size(b);
if m ~= k || o ~= 1
	error('Dimensione del vettore e della matrice non compatibili')
end

for i = 1:n
	alfa = norm(A(i:m, i));
	if alfa == 0
		error('Matrice non a rango massimo');
	end
	if A(i,i) >= 0
		alfa = -alfa;
	end
	v1 = A(i,i) - alfa;
	A(i,i) = alfa;
	A(i+1:m, i) = A(i+1:m, i) / v1;
	beta = -v1 / alfa;
	A(i:m, i+1:n) = A(i:m, i+1:n) - (beta * [1; A(i+1:m, i)]) * ([1 A(i+1:m, i)'] * A(i:m, i+1:n));
	b(i:m) = b(i:m) - (beta * [1 A(i+1:m, i)'] * b(i:m)) * [1; A(i+1:m, i)];
end

x = b(:);
for i = n:-1:1
	x(i) = x(i) / A(i,i);
	x(1:i-1) = x(1:i-1) - A(1:i-1, i) * x(i);
end
nr = norm(x(n+1:m));
x = x(1:n);
return