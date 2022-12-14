% ALTHAMMER Althammer polynomials generated by the Stieltjes algorithm.
% This code is used to generate Figure 3 and Figure 4
close all
clearvars
m=60; g=1/10; s=1; nd=[m m]; a0=0; same=1;
ab=r_jacobi(m); 
zw=gauss(m,ab);
xw=[zw(:,1) zw(:,1) zw(:,2) g*zw(:,2)];


n = 60;
%% Methods by Gautschi and Zhang [Gautschi,Zhang 1995]
% Discretized Stieltjes procedure
B=stieltjes_sob(m,s,nd,xw,a0,same);

% Modified Chebyshev
gamma = g;
mom=zeros(2,2*n);
mom(1,1) = 2; mom(2,1)=2*gamma;
abm = r_jacobi(2*n-1);
B_cheb = chebyshev_sob(n,mom,abm);

%% New proposed methods
% These are based on matrix manipulation, so first the starting vector and
% Jordan matrix are formed. These generate the Krylov subspace

% Starting vector
w = zeros(2*m,1); w(2:2:2*m) = sqrt(zw(:,2));
% Jordan matrix
Z = zeros(2*m);
Z(1:2:end,1:2:end) = diag(zw(:,1));
Z(2:2:end,2:2:end) = diag(zw(:,1));
for k=2*m:-1:1%1:2*m
    if mod(k,2)==1
       Z(k,k+1) = sqrt(g);
    end
end

% Arnoldi iteration
[V,H] = Arnoldi(Z,w,m+1);
% Updating procedure
Hup = updating(Z,w,'PR');

%% plot computed zeros
figure
subplot(2,2,1)
plot(real(sobzeros(n,m,B)),imag(sobzeros(n,m,B)),'b*','MarkerSize',10)
subplot(2,2,2)
plot(real(eig(H(1:n,1:n))),imag(eig(H(1:n,1:n))),'ro','MarkerSize',10)
subplot(2,2,3)
plot(real(sobzeros(n,m,B_cheb)),imag(sobzeros(n,m,B_cheb)),'g+','MarkerSize',10)
subplot(2,2,4)
plot(real(eig(Hup(1:n,1:n))),imag(eig(Hup(1:n,1:n))),'k*','MarkerSize',10)



