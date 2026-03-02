function [H,B,JB] = GegSob_generate(N,mu,lambda,method,balance)
%GegSob_generate generates the recurrence matrix for the Gegenbauer-Sobolev
%orthogonal polynomials as well as the matrix pencil (J*B,B)
%   
%INPUT:
%   N = number of Gegenbauer SOPs in the sequence, i.e., compute {S0,S1,...,S_{N-1}}
%   mu = parameter of the Gegenbauer measure (1-x^2)^mu
%   lambda = weight of the integral involving the derivative
%   The inner product: int(f(x)g(x)(1-x^2)^mu dx + lambda int(f'(x)g'(x)(1-x^2)^mu dx
%   method = select which method to use to compute the parameter alpha:
%          (default) 1 = expression that relates alpha to a sequence of OPs
%                    2 = the analytical recurrence relation for alpha
%   balance = boolean indicating to use balancing for eigenvalue problem:
%                    false = no balancing, then recurrence for monic poly
%         (default)  true = balance the matrix and pencil
%OUTPUT:
%   H = Lower-Hessenberg matrix of size NxN
%   B = The 3-banded matrix B of size NxN (first term in matrix pencil)
%   JB = product of J (Jacobi matrix) with B (3-banded matrix) of size
%           NxN (second term in matrix pencil, after truncation)

if nargin<4
    method = 1;
    balance = true;
elseif nargin <5
    balance = true;
end

% Compute the parameters alpha and gamma
alphas = computeAlpha(N-2,mu,lambda,method);
nn = 1:N-1;
if mu==1/2
    gamma = [1/2,1/4*ones(1,N-2)];
else
    gamma = nn.*(nn+2*(mu-1))./(4*(nn+(mu-1)+1/2).*(nn+(mu-1)-1/2));
end

if balance == false
    % Construct the 3-banded matrix B
    B = diag(ones(N,1));
    B = B+ diag(alphas(1:N-2),-2);

    % Construct the product J*B directly
    JB = diag(ones(N,1),1);  JB = JB(1:N,1:N+1);
    JB = JB + [diag(gamma+alphas,-1),zeros(N,1)];
    JB = JB + [diag(gamma(3:end).*alphas(1:end-2),-3),zeros(N,1)];


    % Construct the Hessenberg recurrence matrix H
    Binv = diag(ones(N,1));
    J = floor((N-1)/2);
    for jj = 1:J
        k = 1:jj;
        for i=1:N-2*jj
            Binv(i+2*jj,i) = (-1)^jj*prod(alphas(i+2*k-2));
        end
    end
    H = Binv*JB;
    H = H(1:N,1:N);
    JB = JB(1:N,1:N);

else % balance

    % Construct the balanced 3-banded matrix B
    B = eye(N,N);
    for ll=1:N-2
        if ll==1
            B(ll+2,ll) = 1/sqrt((gamma(1)+alphas(1))*(gamma(2)+alphas(2)-alphas(1)))*alphas(1);
        else
            B(ll+2,ll) = 1/sqrt((gamma(ll)+alphas(ll)-alphas(ll-1))*(gamma(ll+1)+alphas(ll+1)-alphas(ll)))*alphas(ll);
        end
    end

    % Construct the balanced product J*B directly
    JB = zeros(N,N);
    for ii = 1:N-1
        if ii == 1
            JB(ii,ii+1)  = sqrt(gamma(1)+alphas(1)); JB(ii+1,ii) = sqrt(gamma(1)+alphas(1));
        else
            JB(ii,ii+1)  = sqrt(gamma(ii)+alphas(ii)-alphas(ii-1));JB(ii+1,ii) = 1/sqrt(gamma(ii)+alphas(ii)-alphas(ii-1))*(gamma(ii)+alphas(ii));
        end
    end
    for jj=1:N-3
        if jj==1
            JB(jj+3,jj) = (gamma(3)*alphas(1))/sqrt((gamma(1)+alphas(1))*(gamma(2)+alphas(2)-alphas(1))*(gamma(3)+alphas(3)-alphas(2)));
        else
            JB(jj+3,jj) = (gamma(jj+2)*alphas(jj))/sqrt((gamma(jj)+alphas(jj)-alphas(jj-1))*(gamma(jj+1)+alphas(jj+1)-alphas(jj))*(gamma(jj+2)+alphas(jj+2)-alphas(jj+1)));
        end
    end

    % Construct the balanced Hessenberg recurrence matrix H
    H = B\JB;
end

