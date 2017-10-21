addpath('./util')
%% Descrpition
% 1. verify correctness of QR decomposition by checking
%   a) R upper triangular
%   b) Q orthogonal
%   c) QR = A
% 2. use QR decomposition to solve least square problem

%% verify QR decomposition
A = [1,2;3,4;5,6];
[Q, R] = QR_decomp(A);
% verify whether R is a upper triangular matrix
disp('R =') 
disp(R)
% verify whether Q is orthogonal matrix
disp('Q''*Q =')
disp(Q'*Q)
% verify whether QA = A
disp('QA-A =')
disp(Q*R-A)

%% QR decomposition based solution to least square problem
[m, n] = size(A);
assert(rank(A)==n, ['A should be full rank when solving least square problem'])
b = [1;1;1];
% least square problem(QR decomposition solution)
R_sub = R(1:n,1:n);
Q_t = Q';
Q_t_b_sub = Q_t*b;
Q_t_b_sub = Q_t_b_sub(1:n);
x = inv(R_sub)*Q_t_b_sub;
% least square problem ground truth(normal equation solution)
gt = inv(A'*A)*A'*b;
% verify whether x = gt
disp('x is from QR decomposition, gt is from normal equation')
disp('x-gt =')
disp(x-gt)

