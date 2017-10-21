function [] = main()
addpath('./util')
%% Description
% interpolate f(x)=1/(1+x^2), x \in [-5,5] by order 10 Lagrange polynomial
% with two different set of sample points
% 1. equally spaced points
% 2. zero points of order 11 Chebyshev polynomial

%% ground truth
x = linspace(-5,5);
y = f(x);
plot(x,y)
hold on

%% equally spaced points
xs1 = -5:1:5;
ys1 = f(xs1);
y1 = lagrange(x, xs1, ys1);
plot(x, y1, '--')
hold on

%% zero points of order 11 Chebyshev polynomial
xs2 = 5*zero_point_of_chebyshev(11, 1:11);
ys2 = f(xs2);
y2 = lagrange(x, xs2, ys2);
plot(x, y2, '-.m')
hold off
legend('ground truth', 'equally spaced', 'chebyshev')
title('Lagrange Interpolation of $f(x)=\frac{1}{1+x^2}, x \in [-5, 5]$','Interpreter','LaTex')
end

function [y] = f(x)
y = [];
for n=1:length(x)
    y = [y, 1./(1.+x(n)^2)];
end
end

function [res] = zero_point_of_chebyshev(n, i)
assert(all(i<=n), ['order n chebyshev polynomial has only n zero point'])
res = [];
for k=1:length(i)
    res = [res, cos(((2*i(k)-1)*pi)/(2*n))];
end
end