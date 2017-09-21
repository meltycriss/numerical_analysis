function [H] = mat_householder(v,beta)
%% generate householder matrix from v and beta
% @param v - v[1] = 1
% @param beta - 2/v.'v
% @retval trans - H = I - beta * vv.'

n = length(v);
H = eye(n) - beta*v*v';

end