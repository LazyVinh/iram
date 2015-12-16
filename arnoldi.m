function [V,H,f] = arnoldi(varargin)
%ARNOLDI  Computes the Arnoldi decomposition AV = VH + f*e_k, where A is
%   square, V is isometric, H a k-by-k upper Hessenberg matrix and e_k the
%   k-th unit vector. The decomposition is uniquely determined by V(:,1).
%
%   If A is a (sparse) matrix...
%
%   [V,H,f] = ARNOLDI(A) returns a 6-step Arnoldi decomposition of A.
%
%   [V,H,f] = ARNOLDI(A,k) returns a k-step Arnoldi decomposition of A.
%
%   [V,H,f] = ARNOLDI(A,k,v1) returns the unique k-step Arnoldi
%   decomposition of A with V(:,1) = v1.
%
%   [V,H,f] = ARNOLDI(A,p,V,H,f) returns the expansion of the k-step
%   Arnoldi decomposition to an m-step decomposition with m=k+p.
%
%   If the matrix A is only given by a function handle Ax = AFUN(x)...
%
%   [V,H,f] = ARNOLDI(AFUN,k,v1) returns the unique k-step Arnoldi
%   decomposition of A with V(:,1) = v1. Essentially, the size of v1 tells 
%   us the dimension of A.
%
%   [V,H,f] = ARNOLDI(AFUN,p,V,H,f) returns the expansion of the k-step
%   Arnoldi decomposition to an m-step decomposition with m=k+p.

[A,Amatrix,k,p,V,H,f] = check_input(varargin{:});

global TOL_BRKDWN
if isempty(TOL_BRKDWN)
    TOL_BRKDWN = 1e-12;
    fprintf('No tolerance for the lucky breakdown specified. Setting TOL_BRKDWN := %0.5g.\n',TOL_BRKDWN)
end

for j = 1:p
    beta = norm(f);
    if beta <= TOL_BRKDWN
        disp('LUCKY BREAKDOWN!!!')
        V = V(:,1:k+j-1);
        H = H(1:k+j-1,1:k+j-1);
        return;
    end
    V(:,k+j) = f/beta;
    H(k+j,k+j-1) = beta;
    if Amatrix
        w = A*V(:,k+j);
    else
        w = A(V(:,k+j));
    end
    H(1:k+j, k+j) = V(:,1:k+j)'*w;
    f = w - V(:,1:k+j)*H(1:k+j, k+j);
    % correction steps to conserve orthogonality
    c = V(:,1:k+j)'*f;
    f = f - V(:,1:k+j)*c;
    H(1:k+j, k+j) = H(1:k+j, k+j) + c;
end

end

function [A,Amatrix,k,p,V,H,f] = check_input(varargin)

% Process the input A or AFUN

if isfloat(varargin{1})
    A = varargin{1};
    Amatrix = true;
else
    % By checking the function A with fcnchk, we can now use direct
    % function evaluation on the result, without resorting to feval
    A = fcnchk(varargin{1});
    Amatrix = false;
end
if Amatrix
    [m,n] = size(A);
    if (m ~= n)
        error('MATLAB:arnoldi:NonSquareMatrixOrFunction',...
            'A must be a square matrix or a function.')
    end
    k = 1;
    if nargin == 1
        p = 5;
        
        V = zeros(n,k+p);
        H = zeros(k+p);
        V(:,1) = rand(n,1);
        V(:,1) = V(:,1)/norm(V(:,1));
        H = V(:,1)'*A*V(:,1);
        f = A*V(:,1) - V(:,1)*H;
        return
    end
    if nargin == 2
        p = varargin{2}-1;
        
        V = zeros(n,k+p);
        H = zeros(k+p);
        V(:,1) = rand(n,1);
        V(:,1) = V(:,1)/norm(V(:,1));
        H = V(:,1)'*A*V(:,1);
        f = A*V(:,1) - V(:,1)*H;
        return
    end
    if nargin == 3
        p = varargin{2}-1;
        
        V = zeros(n,k+p);
        H = zeros(k+p);
        v1 = varargin{3};
        if abs(norm(v1)-1)>eps
            warning('MATLAB:arnoldi:StartVectorNotNormalized', 'Normalizing start vector, because norm is not 1.')
            v1 = v1/norm(v1);
        end
        H(1,1) = v1'*A*v1;
        f = A*v1 - v1*H(1,1);
        return
    end
    if nargin > 3
        p = varargin{2};
        V_ = varargin{3};
        H_ = varargin{4};
        [k,k2] = size(H_);
        if k~=k2
            error('MATLAB:arnoldi:HessenbergNotSquare',...
                'H is not square.')
        end
        V = zeros(n,k+p);
        H = zeros(k+p);
        
        V(:,1:k) = V_;
        H(1:k,1:k) = H_;
        f = varargin{5};
        return
    end
else
    if nargin <= 2
        error('MATLAB:arnoldi:NotEnoughInputs', 'Not enough input arguments.')
    end
    if nargin == 3
        k=1;
        p = varargin{2}-1;
        v1 = varargin{3};
        if abs(norm(v1)-1)>eps
            warning('MATLAB:arnoldi:StartVectorNotNormalized', 'Normalizing start vector, because norm is not 1.')
            v1 = v1/norm(v1);
        end
        n = numel(v1);
        V = zeros(n,k+p);
        V(:,1) = v1;
        H = zeros(k+p);
        H(1,1) = v1'*A(v1);
        f = A(v1) - v1*H(1,1);
        return
    end
    if nargin > 3
        p = varargin{2};
        f = varargin{5};
        n = numel(f);
        V_ = varargin{3};
        H_ = varargin{4};
        [k,k2] = size(H_);
        if k~=k2
            error('MATLAB:arnoldi:HessenbergNotSquare',...
                'H is not square.')
        end
        V = zeros(n,k+p);
        H = zeros(k+p);
        V(:,1:k) = V_;
        H(1:k,1:k) = H_;
        return
    end
end
end