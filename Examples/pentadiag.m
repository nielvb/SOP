% PENTADIAG computes the higher order recurrence relation, represented as a
% pentadiagonal matrix, as described in Section 5.1.2
clearvars
close all
%% Define the problem and construct Gauss quadrature rule
m=6; 
alpha = 0; 
ab=r_laguerre(m,alpha); 
zw=gauss(m,ab);

M = 1;
N = 1;

c = -1;

%% Construct the matrix and vector correspoding to the discretized inner product
% Construct starting vector
w = sqrt(zw(:,2));
w = [w;0;sqrt(M)];
% Construct Jordan matrix
Z = diag(zw(:,1));
Z = blkdiag(Z,[c,sqrt(N/M);0,c]);

m = m+2;
% shift
Z = Z-c*eye(m);

%% Method 3: Arnoldi iteration
[V,H] = Arnoldi(Z,w,m-2);
H = H(1:end-1,1:end);
%H = H(1:m,1:m);
%% Method 4: Updating procedure 
Hup = updating(Z,w,'PR');

%% pentadiag


HH = H^2; HH = HH(1:5,1:5);
HHup = Hup^2; HHup = HHup(1:5,1:5);

penta = [5/2,11/(2*sqrt(2)), 1/2*sqrt(89/2),0,0;...
        11/(2*sqrt(2)), 19/2, 129/sqrt(89), 1/2*sqrt(35705/178),0;...
        1/2*sqrt(89/2), 129/sqrt(89), 5331/178, 1503493/(178*sqrt(71410)), 4*sqrt(26690173/3177745);...
        0, 1/2*sqrt(35705/178),1503493/(178*sqrt(71410)), 415128273/6355490, 72140663342/35705*sqrt(2/2375425397);...
        0, 0, 4*sqrt(26690173/3177745), 72140663342/35705*sqrt(2/2375425397),108116532681297/952972626965];

norm(HH-penta,"fro")/norm(penta,'fro')
norm(HHup-penta,"fro")/norm(penta,"fro")