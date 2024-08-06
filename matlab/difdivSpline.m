function df = difdivSpline(X, Y, it)
%
%	df = difdivSpline(X, Y, it)
%
%	Calcola le differenze divise sulle coppie (xi, fi)
%	fermandosi alla it-esima iterazione
%
%	Input:
%		x: vettore delle ascisse
%		f: vettore delle ordinate
%		it: numero di iterazioni
%	Output:
%		df: vettore delle differenze divise

n = length(X);
if length(Y) ~= n
	error('Dati errati');
end
n = n-1;
df = Y;
for j=1:it-1
	for i = n+1:-1:j+1
		df(i) = (df(i) - df(i-1))/(X(i) - X(i-j));
	end
end
df = df(1, it:n+1);
return