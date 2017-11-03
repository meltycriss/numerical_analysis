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
plot(x,y,'linewidth',2)
hold on

%% equally spaced points
xs1 = -5:1:5;
ys1 = f(xs1);
y_lagrange = lagrange(x, xs1, ys1);
plot(x, y_lagrange, '--','linewidth',2)
hold on

%% hermite
[fs, xs] = hermite_preprocess(xs1, @f, @df);
y_hermite = hermite(fs,xs,x);
plot(x,y_hermite,'linewidth',2)
hold on

%% zero points of order 11 Chebyshev polynomial
xs2 = 5*zero_point_of_chebyshev(11, 1:11);
ys2 = f(xs2);
y_chebyshev = lagrange(x, xs2, ys2);
plot(x, y_chebyshev, '-.m','linewidth',2)
hold off

l = legend('ground truth', 'equally spaced', 'hermite' ,'chebyshev');
%title('Lagrange Interpolation of $f(x)=\frac{1}{1+x^2}, x \in [-5, 5]$','Interpreter','LaTex', 'fontsize', 30);
set(l,'Fontsize',30, 'box', 'off');
end

function [y] = f(x)
y = zeros(1, length(x));
for n=1:length(x)
    curr = x(n);
    y(n) = 1./(1.+curr^2);
end
end

function [y] = df(x)
y = zeros(1, length(x));
for n=1:length(x)
    curr = x(n);
    y(n) = (-1.)*(2.*curr)/(1.+curr^2)^2;
end
end

function [res] = zero_point_of_chebyshev(n, i)
assert(all(i<=n), ['order n chebyshev polynomial has only n zero point'])
res = zeros(1, length(i));
for k=1:length(i)
    curr = i(k);
    res(k) = cos(((2*curr-1)*pi)/(2*n));
end
end