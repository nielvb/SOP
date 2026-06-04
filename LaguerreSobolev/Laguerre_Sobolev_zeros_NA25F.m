function [x2,iterj,x4,TT] = Laguerre_Sobolev_zeros_NA25F(n,alpha,lambda,tol)
% Computation of the zeros of
% Laguerre-Sobolev polynomials by the two methods described  in
% T. Laudadio, N. Mastronardi, F. Marcell\'an, N. Van Buggenhout,  P. Van Dooren,
% On computing the zeros of Laguerre-Sobolev polynomials, Numer. Algorithms 100 (2025), 1507--1526.
%
% 
%INPUT
%       n = the degree of the polynomial
%       alpha = the Laguerre parameter
%       lambda = the parameter of the Lagierre-Sobolev inner product
%       tol = the required tolerance
%
%OUTPUT  
%       iterj = vector with the required number of iteration per each zero
%       x2 = vector of zeros computed by the Elrich-Aberth method 
%       conv1 = vector with entries equal to 0 if convergence occurs, 1
%                  if it does not converge in 30 iterations
%       x4 = zeros computed as the generlized eigenvlued =s of the
%            balanced problem
%       TT = execution time of both methods
%
%Author: Nicola Mastronardi (2026)

med11=[];
med00=[];
l6=zeros(n,1);


  tic

        [d,d1,beta,gamma,delta,B,H] = recurrence_Laguerre_Sobolev(n+1,alpha,lambda);
tic
        b1=diag(H);
        c1=diag(H,-1);c1=[0;c1];
        d1=diag(H,-2);d1=[0;0;d1];
        e1=diag(B,-1);e1=[0;e1];

        v(1)=0;
        for i=2:n+1,
            v(i)=c1(i)-e1(i)*(b1(i-1)-e1(i-1));
        end

        delta(1)=1; for i=2:n+1, delta(i)=delta(i-1)*sqrt(v(i));end
        delta1=zeros(n+1,1);
        for i=1:n+1, delta1(i)=sqrt(v(i));end

        Hf=H;Bf=B;
        
        for i=2:n+1
            c1(i)=c1(i)/delta1(i);
            e1(i)=e1(i)/delta1(i);
            a1(i-1)=delta1(i);
        end
        for i=3:n+1
            d1(i)=d1(i)/(delta1(i-1)*delta1(i));
        end
        Hf1=diag(a1(1:n),1)+diag(b1(1:n+1))+diag(c1(2:n+1),-1)+diag(d1(3:n+1),-2);
        Bf1=diag(ones(n+1,1))+diag(e1(2:n+1),-1);
        for i=2:n+1
            Hf(i,i-1)=Hf(i,i-1)/delta1(i);
            Bf(i,i-1)=Bf(i,i-1)/delta1(i);
            Hf(i-1,i)=Hf(i-1,i)*delta1(i); end
        for i=3:n+1
            Hf(i,i-2)=H(i,i-2)/(delta1(i-1)*delta1(i));
        end
        t2=toc;
tic
x4=sort(eig(Hf1(1:n,1:n),Bf1(1:n,1:n)));
 TT(1)=toc+t2;
 tic
[llag,z,a1,b1,mu1]=gaussq(6,n,alpha,0,0,0);
        llag=llag';
        if alpha==0
            l6(1)=0;
        else
            l6(1)=llag(1)/2;
        end
        for i=2:n
            l6(i)=(llag(i-1)+llag(i))/2;
        end

        [iterj,x2] = EA_method_HB(Hf1(1:n,1:n+1),Bf1(1:n,1:n+1),l6,tol);
        TT(2)=toc+t2;
end


function [iterj,x] = EA_method_HB(H,B,x,tol)
% Ehrlich–Aberth Method for computing the 
% zeros  and the left and right eigenvector of the
% Hessenberg matrix aassociated to
% multiple orthogonal polynomials 
%  
x0=x;
[m,n]=size(H);
conv=zeros(m,1);
% V=zeros(m);
% U=zeros(m);
iter=0;
iterj=zeros(m,1);
while sum(conv)<m,
    iter=iter+1;
  %  conv,pause
    for j=1:m
        if conv(j)==0
            x0=x(j);
            [p0,p1] = one_step_Newton_nullspace_HB(H,B,x0);
           
            p01=p0(n)/p1(m);
         %   j,p0(n),p01,disp('vedi'),pause
            sum1=0;
            for k=1:m,
                if k~=j
                    sum1=sum1+1/(x(j)-x(k));
                end
            end
          %  x1(j)=x(j)-p01/(1-p01*sum1); Jacobi
            xt=x(j)-p01/(1-p01*sum1); % Gauss-Seidel
            errj=abs(x(j)-xt);
            x(j)=xt;
            if abs(p01)< tol | abs(p0(n))<tol | errj<tol | iter >30
                conv(j)=1;
               
                iterj(j)=iter;
            end
          
        end
    end

    
end
end
%%%%%%%%%%
% function [H,Q,p,p1a] = one_step_Newton_nullspace(H,lambda)
function [p,p1] = one_step_Newton_nullspace_HB(H,B,lambda)


[m,n]=size(H);

Q=eye(n);
Blam=lambda*B;
H=Blam-H;
H0=H;
CS=zeros(m,2);
p=zeros(n,1);

for i=1:m,
    G=givens(H(i,i),H(i,i+1));
    H(:,i:i+1)= H(:,i:i+1)*G';
  %  Q(i:i+1,:)=G*Q(i:i+1,:);
   CS(i,1)=G(1,1); CS(i,2)=G(1,2);
end

p(n)=CS(n-1,1);
t=1;
for i=n-1:-1:2,
    t=-t*CS(i,2);
    p(i)=CS(i-1,1)*t;
end
p(1)=-t*CS(1,2);

p1=-H0(1:m,2:m+1)\(B(1:m,1:m)*p(1:m));



end

function [d,d1,beta,gamma,delta,B,H] = recurrence_Laguerre_Sobolev(n,alpha,lambda)
% [d,d1,beta,gamma,delta,B,H] = recurrence_Laguerre_Sobolev(n,alpha,lambda)
%
% computation of the four term recurrence relation coefficients
% associated to Laguerre-Sobolev polynomials
%
%INPUT
%   n = degree of the Laguerre-Sobolev polynomial
%   alpha = parameter characterizing Laguerre-Sobolev polynomials: $e^{-x} x^{\alpha}$
%   lambda = parameter charachterizing the Sobolev inner product

beta=zeros(n,1);
gamma=zeros(n-1,1);
delta=zeros(n-2,1);
d(1)=2*(alpha+1)/(lambda+alpha+1);
for i=2:n
    d(i)=((i+1)*(i+alpha))/(i*(2+lambda)+alpha-d(i-1));
end
d=d';

d1(1)=1;
for i=1:n-1
    d1(i+1)=((i+1)*(i+alpha))/((i)*(2+lambda)+alpha-d1(i));
end
d1=d1';
for i=1:n, 
    beta(i)=2*(i-1)+alpha+d1(i);
end
beta=beta';

for i=1:n-1, 
    gamma(i)=(i)*(i+alpha-1)+(2*(i)+alpha)*d1(i);
end
gamma=gamma';

for i=1:n-2,
    delta(i)=(i+1)*(i+1+alpha-1)*d1(i);
end
delta=delta';

B=eye(n);
B=B+diag(d1(1:n-1),-1);
H=diag(ones(n-1,1),1)+diag(beta)+diag(gamma,-1)+diag(delta,-2);
end