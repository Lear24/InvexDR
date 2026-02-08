function [e_cell, e_all] = generate_noise(m,n,sigma,varargin)
% [e_cell, xi_cell] = generate_noise(m,n,sigma,q)
% Assign noise on each client according to sample size and noise variance
% Input:
% m - number of clients
% n - m-varibles (m×1)
% sigma - the variance of noise in each client
% varargin = 'normal', 't' or 'uniform'
% Output:
% e_cell: m-by-1 cell with cell(m_i) is an n-by-1 matrix

e_cell = cell(m,1);
e_all = zeros(sum(n),1);
ni=1;
for i =1:m
    nn = n(i);
    sig = sigma(i);
    e = randn(nn,1)*sig;
    e_cell(i,1)=mat2cell(e,nn,1);
    e_all(ni:ni+nn-1,:) = e;
    ni = ni+ nn;
end




end