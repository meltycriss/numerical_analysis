function [fs, xs] = hermite_preprocess(v, f, df)
%% compute constants for hermite
% @param v - sample point
% @param f - function to be interpolated
% @param df - f'
% @retval fs - f(x_0), f[x_0, x_0], f[x_0, x_0, x_1], ..., f[x_0, x_0, ...,
% x_n, x_n]
% @retval xs - x_0, x_0, ... , x_{n-1}, x_{n-1}, x_n

assert(length(v)>1, 'number of sample points should > 1')
vv = zeros(1, 2*length(v));
for i=1:length(v)
    vv(2*i-1) = v(i);
    vv(2*i) = v(i);
end
xs = vv(1:end-1);

fs = zeros(1, 2*length(v));
fs(1) = f(v(1));
for i=2:length(vv)
    fs(i) = avg_diff(vv,1,i,f,df);
end

end