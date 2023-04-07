function [xk, rho] = BSOR(A, b, w, tol, max_iter, xo)
% Ta funkcja jest numeryczna metoda rozwiazywania ukladow rownan liniowych.
% Jest to implementacja metody BSOR (Backward Successive Over Relaxation).
%
% Wejscia:
% A: macierz nxn, reprezentujaca wspolczynniki ukladow rownan liniowych
% b: wektor nx1, reprezentujacy prawa strone ukladow rownan liniowych
% w: parametr relaksacji, wartosc skalarna pomiedzy 0 a 2
% tol: tolerancja, wartosc skalarna wskazujaca pozadany poziom
% dokladnosci rozwiazania
% max_iter: maksymalna liczba iteracji do wykonania przed zatrzymaniem
% xo: poczatkowe oszacowanie rozwiazania, wektor nx1
%
% Wyjscia:
% xk: przyblizone rozwiazanie ukladow rownan liniowych, wektor nx1
% rho: promien spektralny macierzy iteracji, wartosc skalarna
n = size(A,2);

D = diag(diag(A));
rho = max(abs(eig(inv(D)*(D-A))))
assert(rho < 1, "Metoda nie jest zbiezna")

xk = xo;
for i = 1:n
    xk(i) = xo(i) + 1;
end

k = 0;
while norm(xk-xo) > tol && k <= max_iter
    x = xk;
    for i = n:-1:1
        sum1 = 0;
        sum2 = 0;
        for j = 1:n
            if j < i
                sum1 = sum1 + A(i,j) * xk(j);
            elseif j > i
                sum2 = sum2 + A(i,j) * xo(j);
            end
        end
        xk(i) = (1-w) * xo(i) + w * (b(i) - sum1 - sum2) / A(i,i);
    end
    xo = x;
    k = k+1;
end

D = diag(diag(A));
rho = max(abs(eig(inv(D)*(D-A))));
k
cond = norm(inv(A)) * norm(A)

end