function [ theta, x ] = iram( A, k, p, v1)


% erzeuge k-Schritt Zerlegung
[ Hk, Vk, fk, lucky ] = arnoldi_v1( A, k, v1);

if lucky
    [y,theta] = eig(Hk);
    theta = diag(theta);
    x = Vk*y;
    return
end

m = k+p;

for iter = 1:1000
    
    % expandiere auf m-Schritt Zerlegung
    [ Hm, Vm, fm, lucky ] = arnoldi_v2(A,Hk,Vk,fk,p);
    if lucky
        [y,theta] = eig(Hm);
        theta = diag(theta);
        x = Vm*y;
        return
    end
    
    % checke auf Deflation
    subdiag = diag(Hm,-1);
    
    [subdiag_minval, subdiag_min] = min(abs(subdiag));
    
    rtol = 0.1;
    
    if subdiag_minval <= 0.001
        disp('deflation')
        V1 = Vm(:,1:subdiag_min);
        V2 = Vm(:,subdiag_min:end);
        H1 = Hm(1:subdiag_min,1:subdiag_min);
        H2 = Hm(subdiag_min+1:m,subdiag_min+1:m);
        
        if norm(H1)<0.01
            theta = [];
            return
        end
        
        [y,theta] = eig(H1);
        
        theta = diag(theta);
        
        x = V1*y;
        
        disp(sprintf('Found addtional %d eigenvalues',numel(theta)))
        
        cont = input('continue?')
        
        if cont == 0
            return
        else
            A_rec  = (eye(length(A))-V1*V1')*A;
            [ theta_rec, x_rec ] = iram( A_rec, k, p, Vm(:,subdiag_min+1));
        end
    
        theta = [theta;theta_rec];
        
        x = [x, x_rec];
        
        return
        
    end
    
    % berechne Eigenvektoren und Eigenwerte von Hm
    [y,theta] = eig(Hm);
    theta = diag(theta);
    
    arnoldires = norm(fm)*abs(y(m,:));
    
    % sortiere in gewünschte und ungewünschte Eigenwerte
    %[~,ind] = sort(abs(theta),'ascend');
    
    [~,ind] = sort(theta,'descend');
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
end

theta = theta(ind(1:k));
x = Vm*y(:,ind(1:k));