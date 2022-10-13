% LAGUERRE SOBOLEV Computes and compares computed roots to values available
% in literature.
% This code generates Table 2 and associated errors
clearvars
close all
%% Define the problem and construct Gauss quadrature rule
m=10; g=1; s=1; nd=[m m]; same=1;
alpha = -1/2;  a0=alpha+1;
ab=r_laguerre(m,alpha); 
zw=gauss(m,ab);
xw=[zw(:,1) zw(:,1) zw(:,2) g*zw(:,2)];



%% Discretized Stieltjes procedure [Gautschi, Zhang 1995]
B=stieltjes_sob(m,s,nd,xw,a0,same);


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




%% plot zeros of degree 10 Sobolev orthogonal polynomial
n = 10;
figure
subplot(2,2,1)
plot(real(sobzeros(n,m,B)),imag(sobzeros(n,m,B)),'b*','MarkerSize',10)
title('Stieltjes procedure')
subplot(2,2,2)
plot(real(eig(H(1:n,1:n))),imag(eig(H(1:n,1:n))),'ro','MarkerSize',10)
title('Arnoldi iteration')
subplot(2,2,3)
plot(real(eig(Hup(1:n,1:n))),imag(eig(Hup(1:n,1:n))),'k*','MarkerSize',10)
title('Updating procedure')



%% compute zeros and compare to https://reader.elsevier.com/reader/sd/pii/0377042795002340?token=43EF1A9F427F4E1FFB2DE834B914EB6A1F44ABA829EA2355CE66BA2496F009E44D2F859554B22CFD26AD03EFA72AABC8&originRegion=eu-west-1&originCreation=20220911193847

% Table 1
nulpunten = [];
for degree = 1:10
    nulStieltjes = min(sobzeros(degree,m,B));
    nulArnoldi = min(eig(H(1:degree,1:degree)));
    nulUp = min(eig(Hup(1:degree,1:degree)));
    nulpunten = [nulpunten; [nulStieltjes,nulArnoldi,nulUp]];
end
nulpunten

% Make table
methods = {'Stieltjes','Arnoldi','Updating'};
% Display errors
Tab = array2table(nulpunten,'VariableNames',methods,'RowName',{'1','2','3','4','5','6','7','8','9','10'});
disp(Tab)

%%  Compute error on zeros compared to paper
% Table 1
err = [];
exact = [0.5; 0.051597; -0.070946; -0.087491; -0.07999; -0.068983; ...
            -0.059147; -0.051200; -0.044918; -0.039929];
for degree = 1:10
    nulStieltjes = min(sobzeros(degree,m,B));
    nulArnoldi = min(eig(H(1:degree,1:degree)));
    nulUp = min(eig(Hup(1:degree,1:degree)));
    err = [err; [nulStieltjes,nulArnoldi,nulUp]-exact(degree)];
end


% Make table
methods = {'Stieltjes','Arnoldi','Updating'};
% Display errors
Tab = array2table(err,'VariableNames',methods,'RowName',{'1','2','3','4','5','6','7','8','9','10'});
disp('-----Errors on computed zeros (up to 8 digits)-----')
disp(Tab)