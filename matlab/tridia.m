function x = tridia(a, b, c, x)
%
%	x = tridia(a, b, c, x)
%
%	risolve il sistema tridiagonale
%	b(i)*x(i-1) + a(i)*x(i) + c(i)*x(i+1) = d(i), i = 1...n
%	con x(0)=x(n+1)=0

n = length(a);
for i = 1:n-1
	b(i) = b(i)/a(i);
	a(i+1) = a(i+1) - b(i)*c(i);
	x(i+1) = x(i+1) - b(i)*x(i);
end
	x(n) = x(n)/a(n);
for i = n-1:-1:1
	x(i) = (x(i) - c(i)*x(i+1))/a(i);
end
return