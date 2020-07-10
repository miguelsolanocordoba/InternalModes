function [C,Cg,L,Weig,Ueig] = compute_eigen(rho,zf,f,om)
%[C,Cg,L,Weig,Ueig] = COMPUTE_EIGEN computes eigenvalues using Ashok & Bhaduria (2009)
%
% S=COMPUTE_EIGEN(rho,zf,f,om) solves the omega-constant eigenvalue problem 
% using the finite difference scheme from Ashok & Bhaduria, which solves for 
% the vertical eigenfunction (Weig) at the cell faces of a non-uniform vertical 
% grid and the horizontal eigenfunctions (Ueig) are computed as the dWeig/dz. 
% 
% INPUT: 
% rho  := mean density at grid centers [1xN] 
% zf   := grid points at grid faces [1xN+1]
% f    := Coriolis frequency (constant)
% om   := anglular frequency (constant)
% 
% OUTPUT
% L    := wavelength 
% C    := phase speed (celerity)
% Cg   := group speed
% Weig := vertical eigenvalues (at cell faces 'zf')
% Ueig := horizontal eigenvalues (at cell centers 'zc')
%
% Created: Oladeji, S. June 2020
% Modified: Solano, M. July 2020

%om = 2*pi/(T*3600); % angular freq
%f = 2*7.29*(10^-5)*sin(lat*pi/180); % coriolis freq

dz = -diff(zf); % layer thickness vector (top to bottom)
zc = zf(1:end-1)/2 + zf(2:end)/2; % layer centers' depth vector
drho = rho(2:end) - rho(1:end-1);
dzc = zc(2:end) - zc(1:end-1);
N2 = -(9.81/1025)*(drho./dzc);  % brunt vaisala frequency
N2 = [N2(1);N2;N2(end)];

% Compute parameters
H = nansum(dz);
N = numel(dz);
B = diag (-N2(2:end-1))/(om^2-f^2);

% Main solver
A = [];
A(1,1) = -2*1/(dz(1)*dz(2));
A(1,2) = 2*1/(dz(2)*(dz(1) + dz(2)));
for i = 2:N-2
    A(i,-1+i) =  2*1/(dz(i)*(dz(i) + dz(i+1)));
    A(i, 0+i) = -2*1/(dz(i)*dz(i+1));
    A(i, 1+i) =  2*1/(dz(i+1)*(dz(i) + dz(i+1)));
end
i = i+1;
A(i,-1+i) =  2*1/(dz(i)*(dz(i) + dz(i+1)));
A(i, 0+i) = -2*1/(dz(i)*dz(i+1));

ll = size(A,1);

% Solve the EVP
[W1,k2] = eig(A,B);

k = abs(sqrt(diag(k2)));
[k,Is] = sort(k,'ascend');
W1 = W1(:,Is);
C  = om./k;
L  = 2*pi./k;
Cg = sqrt((om^2-f^2)./(k.^2));

%% Get W structure functions
W2 = [zeros(1,ll); W1; zeros(1,ll)];    % W at cell faces

%% Get U structure
dW2 = W2(2:end,:)-W2(1:end-1,:);
dzu = repmat(dz,[1,size(dz,1)-1]);
Ueig1 = dW2./dzu;
AA = repmat(sum(Ueig1.^2.*dzu,1)./H,[N 1]).^(1/2);
AA(AA==0) = Inf;
Ueig2 = Ueig1./AA;
Ueig2(:,Ueig2(N,:)<0) = -Ueig2(:,Ueig2(N,:)<0);

% Prepare for Weig & Ueig
Ueig = Ueig2;
Weig = W2;

