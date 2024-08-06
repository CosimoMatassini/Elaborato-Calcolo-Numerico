function ll = lebesgue(a,b, nn, type)
%
%	ll = lebesgue(a, b, nn, type)
%
%	restituisce le approssimazioni della costante di Lebesgue
%	sull'intervallo [a, b] per i polinomi di grado specificato da nn
%
%	Input:
%		a, b: intervallo
%		nn: vettore contenente i gradi dei polinomi per l'approssimazione
%		type: ascisse equidistanti (= 0), chebyshev (= 1)
%
%	Output:
%		ll: approssimazioni della costante di Lebesgue

max=10001;
x=linspace(a, b, max);
ll = nn;
for i=1:length(nn)
	if type == 0
		xi = linspace(a, b, nn(i));
	elseif type == 1
		xi = chebyshev(nn(i), a, b);
	end
	leb = zeros(1, max);
	for j=1:nn(i)
		leb = leb + abs(lin(x, xi, j));
	end
	ll(i) = norm(leb);
end
return