function [] = main()
addpath('./util')
%% Description
% compute factors for different interpolation method
% error analysis
% plot

%% interpolation: compute factors
% lagrange
xs_lag = -5:1:5;
ys_lag = f(xs_lag);

% hermite
[fs, xs] = hermite_preprocess(xs_lag, @f, @df);

% zero points of order 11 Chebyshev polynomial
xs_cheb = 5*zero_point_of_chebyshev(11, 1:11);
ys_cheb = f(xs_cheb);

% piecewise lagrange
v_xs_pw_lag = cell(1,length(xs_lag)-1);
v_ys_pw_lag = cell(1,length(xs_lag)-1);
for i=1:length(v_xs_pw_lag)
    v_xs_pw_lag{i} = xs_lag(i:i+1);
    v_ys_pw_lag{i} = ys_lag(i:i+1);
end

% piecewise hermite
v_fs_pw_hermite = cell(1,length(xs_lag)-1);
v_xs_pw_hermite = cell(1,length(xs_lag)-1);
for i=1:length(v_fs_pw_hermite)
    [v_fs_pw_hermite{i}, v_xs_pw_hermite{i}] = hermite_preprocess(xs_lag(i:i+1), @f, @df);
end

%% error analysis
err_x = linspace(-5,5,101);
err_lagrange = mean(abs(f(err_x)-lagrange(xs_lag, ys_lag, err_x)));
err_hermite = mean(abs(f(err_x)-hermite(fs,xs,err_x)));
err_chebyshev = mean(abs(f(err_x)-lagrange(xs_cheb, ys_cheb, err_x)));
err_pw_lagrange = mean(abs(f(err_x)-piecewise_lagrange(v_xs_pw_lag, v_ys_pw_lag, err_x)));
err_pw_hermite = mean(abs(f(err_x)-piecewise_hermite(v_fs_pw_hermite, v_xs_pw_hermite, err_x)));
fprintf([
        'mean absolute error of different interpolation methods:\n'...
        'lagrange           = %.4f\n'...
        'hermite            = %.4f\n'...
        'chebyshev          = %.4f\n'...
        'piecewise lagrange = %.4f\n'...
        'piecewise hermite  = %.4f\n'...
        ], err_lagrange, err_hermite, err_chebyshev, err_pw_lagrange, err_pw_hermite)

%% plot
need_plot = true; %set this var to true if you want to plot
if need_plot
    y_lim = [-.5, 4.];
    
    % sample points
    x = linspace(-5,5,1000);
    y = f(x);
    y_lagrange = lagrange(xs_lag, ys_lag, x);
    y_hermite = hermite(fs,xs,x);
    y_chebyshev = lagrange(xs_cheb, ys_cheb, x);
    y_pw_lagrange = piecewise_lagrange(v_xs_pw_lag, v_ys_pw_lag, x);
    y_pw_hermite = piecewise_hermite(v_fs_pw_hermite, v_xs_pw_hermite, x);
    
    % all
    figure(1)
    plot(x,y,'linewidth',2, 'DisplayName', 'ground truth');
    hold on;
    plot(x, y_lagrange, '--','linewidth',2, 'DisplayName', 'lagrange');
    hold on;
    plot(x,y_hermite, ':', 'linewidth',2, 'DisplayName', 'hermite');
    hold on;
    plot(x, y_chebyshev, '-.','linewidth',2, 'DisplayName', 'chebyshev');
    hold on;
    plot(x, y_pw_lagrange, '--', 'linewidth', 2, 'DisplayName', 'p.w. lagrange');
    hold on;
    plot(x, y_pw_hermite, '-.', 'linewidth', 2, 'DisplayName', 'p.w. hermite');
    hold off;
    l = legend('show');
    set(l,'Fontsize',30, 'box', 'off');

    % lagrange
    figure(2)
    plot(x,y,'linewidth',2, 'DisplayName', 'ground truth');
    hold on;
    plot(x, y_lagrange, '--','linewidth',2, 'DisplayName', 'lagrange');
    hold off;
    l = legend('show');
    ylim(y_lim);
    set(l,'Fontsize',30, 'box', 'off');

    % hermite
    figure(3)
    plot(x,y,'linewidth',2, 'DisplayName', 'ground truth');
    hold on;
    plot(x,y_hermite, '--', 'linewidth',2, 'DisplayName', 'hermite');
    hold off;
    l = legend('show');
    ylim(y_lim);
    set(l,'Fontsize',30, 'box', 'off');

    % chebyshev
    figure(4)
    plot(x,y,'linewidth',2, 'DisplayName', 'ground truth');
    hold on;
    plot(x, y_chebyshev, '--','linewidth',2, 'DisplayName', 'chebyshev');
    hold off;
    l = legend('show');
    ylim(y_lim);
    set(l,'Fontsize',30, 'box', 'off');

    % p.w. lagrange
    figure(5)
    plot(x,y,'linewidth',2, 'DisplayName', 'ground truth');
    hold on;
    plot(x, y_pw_lagrange, '--', 'linewidth', 2, 'DisplayName', 'p.w. lagrange');
    hold off;
    l = legend('show');
    ylim(y_lim);
    set(l,'Fontsize',30, 'box', 'off');

    % p.w. hermite
    figure(6)
    plot(x,y,'linewidth',2, 'DisplayName', 'ground truth');
    hold on;
    plot(x, y_pw_hermite, '--', 'linewidth', 2, 'DisplayName', 'p.w. hermite');
    hold off;
    l = legend('show');
    ylim(y_lim);
    set(l,'Fontsize',30, 'box', 'off');
end
    
    
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