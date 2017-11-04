function [y] = piecewise_hermite(v_fs, v_xs, x)
%% f(x) = H_i(x)
% @param v_fs - vector of fs
% @param v_xs - vector of xs
% @param x - input
% @param y - f(x)

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
    fs = v_fs{n};
    xs = v_xs{n};
    y(k) = hermite(fs, xs, curr_x);
end

end