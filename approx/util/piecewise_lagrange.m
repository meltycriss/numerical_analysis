function [y] = piecewise_lagrange(v_xs, v_ys, x)
%% f(x) = L_i(x)
% @param v_xs - {[x_{00}, x_{01}, ...], [x_{10}, x_{11}, ..], ...}
% @param v_ys - {[y_{00}, y_{01}, ...], [y_{10}, y_{11}, ..], ...}
% @param x - input
% @retval y - f(x)

y = zeros(1,length(x));
for k=1:length(x)
    curr_x = x(k);
    n = 1;
    for i=1:length(v_xs)
        curr_xs = v_xs{i};
        l = curr_xs(1);
        r = curr_xs(end);
        if curr_x>=l && curr_x<=r
            n = i;
            break;
        end
    end
    xs = v_xs{n};
    ys = v_ys{n};
    y(k) = lagrange(curr_x, xs, ys);
end

end