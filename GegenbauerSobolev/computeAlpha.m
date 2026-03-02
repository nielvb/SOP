function alphas = computeAlpha(N,mu,lam,method)
%COMPUTEALPHA computes the parameter alpha in the 3-banded matrix of the
%matrix decomposition of the Hessenberg recurrence matrix of
%Gegenbauer-Sobolev polynomials
%INPUT:
%   N = largest integer for which to compute the alphas, i.e., 0,1,2,...,N
%   mu = Gegenbauer weight (must be >0) in integrals of Sobolev inner product
%   lam = parameter lambda in Sobolev inner product, positive real number
%   method = select which method to use to compute the parameter alpha:
%          (default) 1 = expression that relates alpha to a sequence of OPs
%                    2 = the analytical recurrence relation for alpha
%OUTPUT
%   alphas = alphas for 0,1,2,...,N


c = @(n,mu) n*(n-1)/(4*(n+mu-1/2)*(n+mu-3/2));
gam = @(n,mu) n*(n+2*mu)/(4*(n+mu+1/2)*(n+mu-1/2));

if mu<=0
    error('Invalid value for mu, must satisfy mu>0')
end

if nargin<4
    method = 1;
end
alphas = [];
if method == 1 % Expression in terms of orthogonal polynomials
    a = @(n) (n+1/2)/n*(n+mu)/(n+mu+1/2)*1/(256*(n+mu/2+3/4)*(n+mu/2+1/4)^2*(n+mu/2-1/4));
    b = @(n) 1/16*1/(n+mu/2+3/4)*((n+1/2)/((n+mu/2+1/4)*(n+mu+1/2))+(n+mu+1)/((n+1)*(n+mu/2+5/4)));
    atil = @(n) n/(n-1/2)*(n+mu-1/2)/(256*(n+mu)*(n+mu/2+1/4)*(n+mu/2-1/4)^2*(n+mu/2-3/4));
    btil = @(n) 1/(16*(n+mu/2+1/4))*(n/((n+mu)*(n+mu/2-1/4))+(n+mu+1/2)/((n+1/2)*(n+mu/2+3/4)));
    for ii=0:N
        if ii==0
            alpha = -c(2,mu);
        else
            if mod(ii,2)==0 % even indices
                n = ii/2;
                Phat_prev = 1;
                Phat_n = lam+(mu+1)/((2*mu+5)*(2*mu+3));
                if n>1
                    for ii=1:n-1
                        Phat_next = (lam+b(ii))*Phat_n-a(ii)*Phat_prev;
                        Phat_prev = Phat_n;
                        Phat_n = Phat_next;
                    end
                end
                alpha =-c(2*n+2,mu)*gam(2*n,mu)/(4*n^2)*Phat_prev/Phat_n;
            else % odd indices
                Qhat_prev = 1;
                Qhat_n = lam+1/(2*mu+3);
                n=(ii+1)/2;
                if n>1
                    for ii=1:n-1
                        Qhat_next = (lam+btil(ii))*Qhat_n-atil(ii)*Qhat_prev;
                        Qhat_prev = Qhat_n;
                        Qhat_n = Qhat_next;
                    end
                end
                alpha =-c(2*n+1,mu)*gam(2*n-1,mu)/((2*n-1)^2)*Qhat_prev/Qhat_n;
            end
        end
        alphas = [alphas,alpha];
    end

elseif method == 2 % Analytical 3-term recurrence relation
    for ii = 0:N
        if ii==0
            alpha =  -2/((2*mu+1)*(2*mu+3));
        elseif ii==1
            alpha = -6/((2*mu+3)^2*(2*mu+5)*(lam+1/(2*mu+3)));
        else
            alpprev = alphas(end-1);
            numerator = -c(ii+2,mu)*gam(ii,mu);
            denominator = lam*ii^2 + c(ii,mu)*ii/(ii+2*mu-1)  + alpprev*ii/(ii+2*mu-1)  + gam(ii,mu);
            alpha = numerator/denominator;
        end
        alphas = [alphas,alpha];
    end

end



