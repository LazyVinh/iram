function [ H, V, f, lucky ] = arnoldi_v1( A, k, v1)
%ARNOLDI_V1 Summary of this function goes here
%   Detailed explanation goes here

lucky = 0;

[m,n] = size(A);

if m~=n
    error('A is not square')
end

if nargin == 2
    v1 = rand(n,1);
end

v1 = v1/norm(v1);

V = zeros(n,k);

V(:,1) = v1;

H = zeros(k,k);

for j = 1:k
    for i = 1:j
        H(i,j) = V(:,i)'*(A*V(:,j));
    end
    f = A*V(:,j) - V(:,1:j)*H(1:j,j);
    if j == k
        return
    else
        H(j+1,j) = norm(f);
    end
    if H(j+1,j) ~= 0
        V(:,j+1) = f/H(j+1,j);
%        V(:,1:j+1) = orth(V(:,1:j+1));
    else
        disp('LUCKY BREAK DOWN!!!')
        H(j+1,j)
        lucky = 1;
        H = H(1:k+j-1,1:k+j-1);
        V = V(:,1:k+j-1);
        return
    end
end

