function [theta,x] = iram(AFUN,n_eig,v1)
%IRAM  Implicitly restarted Arnoldi algorithm, computes the n_eig largest 
%   eigenvalues in modulus of the matrix A given by A*x = AFUN(x).
%
%   v1 is a column vector of length n, where [m,n] = size(A), m = n.
%
%   TOLERANCES: You can specify tolerances for the lucky breakdown and the
%   deflation by defining the global variables TOL_BRKDWN and TOL_DEFL, 
%   respectively.

if n_eig < 1
    theta = [];
    x = [];
    return
end

global TOL_BRKDWN TOL_DEFL
if isempty(TOL_BRKDWN)
    TOL_BRKDWN = 1e-12;
    fprintf('No tolerance for the lucky breakdown specified. Setting TOL_BRKDWN := %0.5g.\n',TOL_BRKDWN)
end
if isempty(TOL_DEFL)
    TOL_DEFL = 1e-12;
    fprintf('No tolerance for the deflation specified. Setting TOL_DEFL := %0.5g.\n',TOL_DEFL)
end

dim = numel(v1);

k = n_eig;
p = max(6,n_eig);
m = k+p;

% compute k-step Arnoldi decomposition
[Vk,Hk,fk] = arnoldi(AFUN,k,v1);
if norm(fk) <= TOL_BRKDWN
    disp('TERMINATING IRAM BECAUSE OF LUCKY BREAKDOWN.')
    [y,theta] = eig(Hk);
    theta = diag(theta);
    x = Vk*y;
    return
end
for iter = 1:min(dim,1000)

    % expand to 2*n_eig-step decomposition
    [Vm,Hm,fm] = arnoldi(AFUN,p,Vk,Hk,fk);
    if norm(fm) <= TOL_BRKDWN
        fprintf('Termination after %d iterations because of lucky breakdown.\n',iter)
        [y,theta] = eig(Hm);
        theta = diag(theta);
        x = Vm*y;
        return
    end
    % check for deflation
    subdiag = diag(Hm,-1);
    [subdiag_minval, subdiag_min] = min(abs(subdiag));
    
    if subdiag_minval <= TOL_DEFL
        fprintf('Deflation after %d iterations in the %dth subdiagonal entry.\n',iter,subdiag_min)
        V1 = Vm(:,1:subdiag_min);
        V2 = Vm(:,subdiag_min:end);
        H1 = Hm(1:subdiag_min,1:subdiag_min);
        H2 = Hm(subdiag_min+1:m,subdiag_min+1:m);
        
        [y,theta] = eig(H1);
        theta = diag(theta);
        x = V1*y;
        
        %A_factors2 = {V1,V1',A_factors{:}};
        %A_ops2 = [dot,dot,A_ops];
        %[theta_rec,x_rec] = iram(A_factors,A_ops,n_eig-numel(theta), Vm(:,subdiag_min+1));
        
        AFUN_new = @(x) newFunc(V1,AFUN,x);
        
        [theta_rec,x_rec] = iram(AFUN_new,n_eig-numel(theta), Vm(:,subdiag_min+1));
        
        
        theta = [theta;theta_rec];
        x = [x, x_rec];
        
        return
        
    end
    
    theta = eig(Hm);
    
    %arnoldires = norm(fm)*abs(y(m,:));
    
    % apply shifts to weaken the influence of unwanted eigenvalues
    [~,ind] = sort(abs(theta),'descend');
    ind = ind(k+1:m);
    
    et = zeros(m,1);
    et(m) = 1;
    for i = 1:numel(ind)
        [Q,R] = qr(Hm-theta(ind(i))*eye(m));
        Hm = R*Q+theta(ind(i))*eye(m);
        Vm = Vm*Q;
        et = Q'*et;
    end
    
    beta = Hm(k+1,k);
    sigma = et(k);
    fk = beta*Vm(:,k+1) + sigma*fm;
    Hk = Hm(1:k,1:k);
    Vk = Vm(:,1:k);
    
    iter
    Hk
    eig(Hk)
end

[y,theta] = eig(Hk);
theta = diag(theta);
x = Vk*y;

% Helper function for computing (I-V*V')*AFUN(X). This is necessary because
% the product V*V' would be a full matrix with the same size as A.
function ret = newFunc(V,AFUN,X)
AFUNX = AFUN(X);
ret = AFUNX - V*(V'*AFUNX);