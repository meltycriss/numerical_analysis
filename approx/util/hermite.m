function [y] = hermite(fs, xs, x)
%% hermite interpolation: H(x) = f(x_0) + f[x_0, x_0]*(x-x_0) + f[x_0, x_0, x_1]*(x-x_0)^2 + ...
% @param fs - f(x_0), f[x_0, x_0], f[x_0, x_0, x_1], ..., f[x_0, x_0, ...,
% x_n, x_n]
% @param xs - x_0, x_0, ... , x_{n-1}, x_{n-1}, x_n
% @param x - x
% @retval y - H(x)

assert(length(fs)==length(xs)+1);

y = zeros(1, length(x));
for n=1:length(x)
    curr_x = x(n);
    y(n) = fs(1);
    for i=1:length(xs)
        curr_term = 1.;
        for j=1:i
            curr_term = curr_term * (curr_x-xs(j));
        end
        y(n) = y(n) + fs(i+1)*curr_term;
    end
end

end