function x = mialdl(A,b)
%
%	x = mialdl(A,b)
%
%	Risolve Ax = b con fattorizzazione LDLt
%
%	Input:
%		A: matrice n x n
%		b: vettore dei termini noti
%
%	Output:
%		x: soluzione del sistema Ax = b

[m, n] = size(A);
if m ~= n
	error('La matrice non e quadrata!');
end
[m, o] = size(b);
if m ~= n || o ~= 1
	error('Dimensione del vettore e della matrice non compatibili')
end

if A(1,1) <= 0
 error('Matrice non sdp');
end

A(2:n, 1) = A(2:n, 1) / A(1,1);
for j = 2:n
	v = (A(j, 1:j-1).') .* diag(A(1:j-1, 1:j-1));
	A(j, j) = A(j, j) - A(j, 1:j-1) * v;
	if A(j, j) <= 0
		error('Matrice non sdp');
	end
	A(j+1:n, j) = (A(j+1:n, j) - A(j+1:n, 1:j-1) * v) / A(j, j);
end

d = diag(A);
if ~all(d > 0)
	error('Matrice non sdp');
end
x = b(:);
for i = 2:n
	x(i:n) = x(i:n) - A(i:n, i-1) * x(i-1);
end
x = x./d;
for i = n:-1:2
	x(1:i-1) = x(1:i-1) - A(i, 1:i-1)' * x(i);
end
return