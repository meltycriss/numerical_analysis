function [y] = lagrange(x, xs, ys)
%% L(x) = \sum{y_i * l_i(x)}
% @param x - input
% @param xs - x_1, x_2, ...
% @param ys - y_1, y_2, ...

assert(all(size(xs)==size(ys)), ['size of xs and ys should be identical']);
y = [];
for n=1:length(x)
    ls = [];
    for i=1:length(xs)
        ls = [ls,lagrange_base(x(n), xs, i)];
    end
    y = [y, sum(ys.*ls)];
end
