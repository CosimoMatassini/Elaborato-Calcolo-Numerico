function yy = hermite(xi, fi, f1i, xx)
%
%	yy = hermite(xi, fi, f1i, xx)
%
%	Calcola il polinomio interpolante di Hermite definito dalle
%	coppie (xi, yi) nei punti del vettore xx
%
%	Input:
%		(xi, fi, f1i): dati del problema
%		xx: vettore in cui calcolare il polinomio
%
% Output:
%		yy: polinomio interpolante di Hermite

if length(fi) ~= length(xi) || length(xi) <= 0 || length(xi) ~= length(f1i)
	error('Dati inconsistenti');
end

if length(unique(xi)) ~= length(xi)
	error('Le ascisse non sono distinte');
end

fi = repelem(fi, 2);
for i = 1:length(f1i)
	fi(i*2) = f1i(i);
end
df = difdivHermite(xi, fi);
n = length(df)-1;
yy = df(n+1) * ones(size(xx));
for i = n:-1:1
	yy = yy.*(xx - xi(round(i/2))) + df(i);
end
return