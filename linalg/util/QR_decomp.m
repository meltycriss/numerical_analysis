function [Q, R] = QR_decomp(A)
%% A = Q[R;0]
% @param A - m>=n
% @retval Q - orthogonal matrix
% @retval R - upper triangular matrix

[m,n] = size(A);

Q_tmp = eye(m);

% uncommenting matrix size check depends on requirements
%assert(m>=n, ['m<n']);
%for i=1:n
for i=1:min(m, n)
    x = A(:,i);
    y = x;
    y(i) = norm(y(i:m));
    if i+1<=m
        y(i+1:m) = 0;
    end
    [v, beta] = householder(x, y);
    H = mat_householder(v, beta);
    Q_tmp = H*Q_tmp;
    A = H*A;
end

Q = Q_tmp';
R = A;

end