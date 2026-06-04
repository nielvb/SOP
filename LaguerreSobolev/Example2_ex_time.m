% Code that generates Figure 1 in the paper (reduced largest n to 400 for faster execution time, paper goes up to 2000)
%   This figure compares the execution time of three methods.
% Paper: T. Laudadio, N. Mastronardi, F. Marcell\'an, N. Van Buggenhout,  P. Van Dooren,
% On computing the zeros of Laguerre-Sobolev polynomials, Numer. Algorithms 100 (2025), 1507--1526.

clear all
close all
ind=0;
alpha=-.5;
lambda=.5;
tol=1e-13;
    for n = [100 200 300 400] %n=[100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000]
        n
        tol=tol/sqrt(n);
ind=ind+1;
        
      tic
        [iterj1,x1,conv1]=Comput_zeros_LagSobF(n,alpha,lambda,tol);
        Tim(ind,1)=toc;

       
        [x3,iterj3,x4,TT] = Laguerre_Sobolev_zeros_NA25F(n,alpha,lambda,tol);
        Tim(ind,2)=TT(2);
        Tim(ind,3)=TT(1);
    end


    % 
figure(1),semilogy([100:100:n],Tim(:,1),'*',[100:100:n],Tim(:,2),'+',[100:100:n],Tim(:,3),'o')
title('Execution time','Interpreter','latex','FontSize',14)
xlabel('$n$','Interpreter','latex','FontSize',14)
ylabel('$\log(sec)$','Interpreter','latex','FontSize',14)
legend('Method $A$','Method $B$','Method $C$','Interpreter','latex','FontSize',14)
