function [v,beta] = householder(x,y)
%% H = I - beta * vv.'
% @param x - original vector
% @param y - projected vector
% @retval v - v[1] = 1
% @retval beta - 2/dot(v,v)

% norm check
eps = 1e-10;
assert(abs(norm(x)-norm(y))<eps, ['norm(x) != norm(y)']);  

v = x-y;
if any(v(:))
    beta = 2/dot(v,v);
else
    beta = 0
end

end
