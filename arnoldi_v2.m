function [ H,V,f, lucky ] = arnoldi_v2(A,Hk,Vk,fk,p)

lucky = 0;

[m,n] = size(A);

k = length(Hk);

V = zeros(n, k+p);
V(:,1:k) = Vk;

H = zeros(k+p, k+p);
H(1:k,1:k) = Hk;

f = fk;

for j = 1:p
    beta = norm(f);
    if beta == 0
        disp('LUCKY BREAKDOWN!!!')
        beta
        lucky = 1;
        H = H(1:k+j-1,1:k+j-1);
        V = V(:,1:k+j-1);
        return;
    end
    V(:,k+j) = f/beta;
    H(k+j,k+j-1) = beta;
    
    w = A*V(:,k+j);
    
    H(1:k+j, k+j) = V(:,1:k+j)'*w;
    
    f = w - V(:,1:k+j)*H(1:k+j, k+j);
end

end

