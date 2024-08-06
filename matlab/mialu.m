function x = mialu(A,b)
%
%	x = mialu(A,b)
%
%	Metodo di fattorizzazione LU con pivoting parziale
%
%	Input:
%		A: matrice n x n
%		b: vettore dei termini noti
%
%	Output:
%		x: soluzione del sistema Ax = b
%

[m, n] = size(A);
if m ~= n
	error('La matrice non e quadrata!');
end
[m, o] = size(b);
if m ~= n || o ~= 1
	error('Dimensione del vettore e della matrice non compatibili')
end
P = (1:n)';

for i = 1:n-1
	[mi, ki] = max(abs(A(i:n, i)));
	if mi == 0
		error('Matrice singolare');
	end
	ki = ki+i-1;
	if ki > i
		A([i ki], :) = A([ki i], :);
		P([i ki]) = P([ki i]);
	end
	A(i+1:n, i) = A(i+1:n, i) / A(i,i);
	A(i+1:n, i+1:n) = A(i+1:n, i+1:n) - A(i+1:n, i) * A(i, i+1:n);
end

x = b(P);
for i = 2:n
	x(i:n) = x(i:n) - A(i:n, i-1) * x(i-1);
end
for i = n:-1:1
	x(i) = x(i) / A(i,i);
	x(1:i-1) = x(1:i-1) - A(1:i-1,i) * x(i);
end
return