function [y] = avg_diff(v, i, j, f, df)
%% f[x_0, x_1, .., x_{n-1}, x_n] = f[x_1, ..., x_n]-f[x_0, ..., x_{n-1}]/(x_n-x_0)
% @param v - x_0, x_1, ..., x_n
% @param i - begin index
% @param j - end index
% @param f - f
% @param df - f'
% @retval y - {INFI, f[v_i,v_{i+1}], f[v_i, v_{i+1}, v_{i+2}], ..., f[v_i,
% ..., v_j]}

assert(i<j, 'i should less than j');
global INFI;
INFI = 1e10;
global aux_mat;
aux_mat = INFI*ones(length(v));

aux_avg_diff(v, i, j, f, df);
y = aux_mat(i,i:j);
aux_mat = [];
end

function [] = aux_avg_diff(v, i, j, f, df)
% record computed result in aux_mat so as to speed up
global INFI;
global aux_mat;
a = v(i);
b = v(j);
if j==i+1    
    if a==b
        aux_mat(i,j) = df(a);
    else
        aux_mat(i,j) = (f(b)-f(a))/(b-a);
    end
else
    if(aux_mat(i+1, j)==INFI)
        aux_avg_diff(v, i+1, j, f, df);
    end
    if(aux_mat(i, j-1)==INFI)
        aux_avg_diff(v, i, j-1, f, df);
    end
    aux_mat(i,j) = (aux_mat(i+1, j)-aux_mat(i, j-1))/(b-a);
end
end