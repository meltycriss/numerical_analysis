function [y] = lagrange(x, xs, ys)
%% L(x) = \sum{y_i * l_i(x)}
% @param x - input
% @param xs - x_1, x_2, ...
% @param ys - y_1, y_2, ...
% @retval y - L(x)

assert(all(size(xs)==size(ys)), ['size of xs and ys should be identical']);
y = zeros(1, length(x));
for n=1:length(x)
    ls = zeros(1, length(xs));
    for i=1:length(xs)
        ls(i) = lagrange_base(x(n), xs, i);
    end
    y(n) = sum(ys.*ls);
end
