% Example of usage of the code to generate Gegenbauer-Sobolev polynomials

%% Set parameters
mu = 1/2; % Gegenbauer parameter
lambda = 1; % Parameter weighing the contribution of the integral over derivatives
N = 20; % Degree of Gegenbauer-Sobolev polynomial

%% Compute recurrence relation
[H,B,JB] = GegSob_generate(N,mu,lambda); % Computes the recurrence matrix H and pencil (JB,B)

%% Compute zeros
[X,GegSobzeros,Y] = eig(JB,B);
GegSobzeros = diag(GegSobzeros);
disp('--------The zeros of the Gegenbauer-Sobolev polynomials are:------------')
GegSobzeros(:).'


%% Compute metrics, the conditioning of the generalized eigenvalue problem and the defect from normality
conds = [];
for ii = 1:N
    xi_ii = Y(:,ii)'*JB*X(:,ii)/(norm(X(:,ii))*norm(Y(:,ii)));
    phi_ii = Y(:,ii)'*B*X(:,ii)/(norm(X(:,ii))*norm(Y(:,ii)));
    conds = [conds, sqrt(norm(JB)^2+norm(B)^2)/sqrt(abs(xi_ii)^2+abs(phi_ii)^2)];
end
defect = sum(svd(B\JB).^2) + sum(svd(JB/B).^2) - 2*sum(abs(eig(JB,B)).^2);



