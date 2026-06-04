function[iterj1,x1,conv1]=Comput_zeros_LagSobF(n,alpha,lambda,tol)
% Computation of the zeros of
% Laguerre-Sobolev polynomials by the method described in
% T. Laudadio, N. Mastronardi, F. Marcell\'an, N. Van Buggenhout,  P. Van Dooren,
% Computation of the zeros of Laguerre--Sobolev polynomials by  the Ehrlich--Aberth method 
%
%INPUT 
%       n = the degree of the polynomial
%       x = vector of the initial approximation of the zeros
%       alpha = the Laguerre parameter
%       lambda = the parameter of the Laguerre-Sobolev inner product
%       tol = the required tolerance
%
%OUTPUT
%       iterj = vector with the required number of iteration per each zero
%       x1 = vector of zeros
%       conv1 = vector with entries equal to 0 if convergence occurs, 1
%          if it does not converge in 30 iterations
%
%Author: Nicola Mastronardi (2026)

%%%%%%%%%%%%%
 
[g2F,z,a1,b1,mu1]=gaussq(6,n,alpha,alpha,0,0);

for i=1:n-1
    z(i+1)=(g2F(i)+g2F(i+1))/2;
end
z(1)=g2F(1)/2;
[iterj1,x1,conv1] = EA_method_F(n,z,alpha,lambda,tol); % Linear 



function [iterj,x,conv1] = EA_method_F(n,x,alpha,lambda,tol)
% Ehrlich–Aberth Method for computing the 
% zeros  of Laguerre-Sobolev polynomials
%INPUT 
%       n = the degree of the polynomial
%       x = vector of the initial approximation of the zeros
%       alpha = the Laguerre parameter
%       lambda = the parameter of the Lagierre-Sobolev inner product
%       tol = the required tolerance
%
%OUTPUT
%       iterj = vector with the required number of iteration per each zero
%       x = vector of zeros
%       conv1 = vector with entries equal to 0 if convergence occurs, 1
%               if it not converges in 30 iterations

x0=x;
X=[];
conv=zeros(n,1);
conv1=zeros(n,3);
iter=0;
iterj=zeros(n,1);
while sum(conv)<n
    iter=iter+1;
    %  conv,pause
    for j=1:n
        if conv(j)==0
            x0=x(j);
            [QN] = poly_LagueSob_ratioN(n,alpha,lambda,x0);
            
            p01=QN(n);
            sum1=0;
            for k=1:n
                if k~=j
                    sum1=sum1+1/(x(j)-x(k));
                end
            end
            %  x1(j)=x(j)-p01/(1-p01*sum1); Jacobi
            xt=x(j)-p01/(1-p01*sum1); % Gauss-Seidel
            errj=abs(x(j)-xt);
            x(j)=xt;
            if abs(p01)< tol | errj<tol | iter >30
                conv(j)=1;
                iterj(j)=iter;
            end
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [QN] = poly_LagueSob_ratioN(n,alpha,lambda,x)
% Computation of the ratio Q_n(x)/Q'_n(x) in x
% Q_n(x) Laguerre-Sobolev polynomials
%INPUT
%       n = the degree of the polynomial
%       alpha = the Laguerre parameter
%       lambda = the parameter of the Lagierre-Sobolev inner product
%       x = point in which the ratio is computed
%
%OUTPUT  
%       QN = the value of the  ratio Q_n(x)/Q'_n(x)   
        for i =2:n+1
            beta1(i)=2*(i-1)+alpha+1;
            gamma1(i)=(i-1)*((i-1)+alpha);
        end
        
        d(1)=1;
        for i=1:n
            d(i+1)=((i+1)*(i+alpha))/(i*(2+lambda)+alpha-d(i));
        end

        M(1)=1/(x-(alpha+1));
        for i=2:n
            M(i)=1/(x-beta1(i)-gamma1(i)*M(i-1));
        end


        N(1)=x-(alpha+1);
        D(1)=0;
        K(1)=M(1);
        S(1)=0;
        for i=2:n
            N(i)=N(i-1)/(M(i)*(x-beta1(i)-gamma1(i)*D(i-1)+N(i-1)));
            D(i)=1/(x-beta1(i)-gamma1(i)*D(i-1)+N(i-1));
        end
        for i=2:n
            KK=(1+i*M(i))/M(i)*((1+d(i-1)*K(i-1))/(1+(i-1)*M(i-1)))-d(i);
            K(i)=1/KK;
            SS=(1+i*D(i))/D(i)*((1+d(i-1)*S(i-1))/(1+(i-1)*D(i-1)))-d(i);
            S(i)=1/SS;
        end

        QN(1)=x-(alpha+1);
        for i=2:n
            QN(i)=N(i)*((1+i*M(i))/(1+i*D(i)))*((1+d(i)*S(i))/(1+d(i)*K(i)));
        end