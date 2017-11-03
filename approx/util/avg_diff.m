function [y] = avg_diff(v, i, j, f, df)
%% f[x_0, x_1, .., x_{n-1}, x_n] = f[x_1, ..., x_n]-f[x_0, ..., x_{n-1}]/(x_n-x_0)
% @param v - x_0, x_1, ..., x_n
% @param i - begin index
% @param j - end index
% @param f - f
% @param df - f'
% @retval y - f[x_0, x_1, .., x_{n-1}, x_n]

assert(i<j, 'i should less than j');
a = v(i);
b = v(j);
if j==i+1    
    if a==b
        y = df(a);
    else
        y = (f(b)-f(a))/(b-a);
    end
else
    y = (avg_diff(v, i+1, j, f, df)-avg_diff(v, i, j-1, f, df))/(b-a);
end
end