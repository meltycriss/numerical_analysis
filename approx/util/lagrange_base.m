function [y] = lagrange_base(x, xs, i)
%% l_i(x) = \sum{(x-x_1)...(x-x_{i-1})(x-x{i+1}...(x-x_{n+1}))} / \sum{(x_i-x_1)...(x_i-x_{i-1})(x_i-x{i+1}...(x_i-x_{n+1}))}
% @param x - input x
% @param xs - x_1, x_2, ...
% @param i - 1-based ith lagrange base
% @retval y - value of l_i(x)

numerator = x-xs;
numerator(i) = [];
numerator = prod(numerator);

denominator = xs(i)-xs;
denominator(i) = [];
denominator = prod(denominator);

y = numerator/denominator;
end