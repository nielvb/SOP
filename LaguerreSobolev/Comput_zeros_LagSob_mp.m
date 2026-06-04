function[iterj1,x1,conv1]=Comput_zeros_LagSob_mp(n,alpha,lambda,tol,ndigit)
% Computation of the zeros of
% Laguerre-Sobolev polynomials by   the Ehrlich--Aberth method 
% by using the four term recurrence relation 
%in multiple precision with Advanpix LLC. Multiprecision Computing Toolbox for Matlab.
%
% input: n: the degree of the polynomial
%        alpha: the Laguerre parameter
%        lambda: the parameter of the Lagierre-Sobolev inner product
%        tol: the required tolerance
%        ndigit: number of digit in multiple precision
%
% output:  iterj: vector with the required number of iteration per each zero
%              x1: vector of zeros
%          conv1: vector with entries equal to 0 if convergence occurs, 1
%          if it not converges in 30 iterations^

path_to_advanpix = '';
if isempty(path_to_advanpix)
    error("Provide path to advanpix (extended precision) or use code that does not require advanpix (double precision code)")
end
addpath(path_to_advanpix);
mp.Digits(ndigit); 
[g2F,z,a1,b1,mu1]=gaussq(6,n,alpha,alpha,0,0);

for i=1:n-1
    z(i+1)=(g2F(i)+g2F(i+1))/2;
end
z(1)=g2F(1)/2;
z=mp(z);
[iterj1,x1,conv1] = EA_method_F_mp(n,z,alpha,lambda,tol,ndigit); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [iterj,x,conv1] = EA_method_F_mp(n,x,alpha,lambda,tol,ndigit)
% Ehrlich–Aberth Method for computing the 
% zeros  of Laguerre-Sobolev polynomials by   the Ehrlich--Aberth method 
% 9/1/2026
% 

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
            
            [p,p1,d,q,q1] = poly_LagueSob_mp(n,alpha,lambda,x0);
            
             p01=q(n+1)/q1(n+1);
          

            sum1=0;
            for k=1:n
                if k~=j
                    sum1=sum1+mp('1')/(x(j)-x(k));
                end
            end
            %  x1(j)=x(j)-p01/(1-p01*sum1); Jacobi
            xt=x(j)-p01/(mp('1')-p01*sum1); % Gauss-Seidel
            errj=abs(x(j)-xt);
            x(j)=xt;
            
            if abs(p01)< tol | errj<tol | iter >30
                conv(j)=1;
               
                iterj(j)=iter;
            end
           
        end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [p,p1,d,q,q1] = poly_LagueSob_mp(n,alpha,lambda,xc,ndigit)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%  addpath('D:\Nicola\matlab\Multiprecision Computing Toolbox\');
%  mp.Digits(ndigit);
 %alpha=mp('-1/2');
 % lambda=mp('10^8');
    x=mp(xc);
    beta1=mp(zeros(n+1,1));
    gamma1=mp(zeros(n+1,1));
    p=mp(zeros(n+1,1));
    p1=mp(zeros(n+1,1));
    p=mp(zeros(n+1,1));
    p(1)=mp('1');
    p1(1)=mp('0');
    p(2)=x-(alpha+mp('1'));
    p1(2)=mp('1');
    
    for i =2:n+1
        beta1(i)=mp('2')*(i-1)+alpha+mp('1');
        gamma1(i)=(i-1)*((i-1)+alpha);
        p(i+1)=x*p(i)-beta1(i)*p(i)-gamma1(i)*p(i-1);
       p1(i+1)=x*p1(i)+p(i)-beta1(i)*p1(i)-gamma1(i)*p1(i-1);
 
    end
 d1(1)=mp('1'); 
 for i=1:n
     d1(i+1)=((i+1)*(i+alpha))/(i*(2+lambda)+alpha-d1(i));
 end
d(1)=mp('1'); 
for i=1:n
    d(i+1)=((i+1)*(i+alpha))/(i*(2+lambda)+alpha-d(i));
end
%d,d1,alpha,pause,
q(1)=p(1);
q(2)=p(2);
q1(1)=mp('0');
q1(2)=mp('1');

for i=2:n
    q(i+1)=p(i+1)+(i)*p(i)-d(i)*q(i);
     q1(i+1)=p1(i+1)+(i)*p1(i)-d(i)*q1(i);
end

