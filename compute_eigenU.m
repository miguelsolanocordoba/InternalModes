function S = compute_eigenU(rho,rho0,zf,f,om)
%S = COMPUTE_EIGENU computes eigenvalues using Ashok & Bhaduria (2009)
%
% S=COMPUTE_EIGENU(rho,zf,f,om) solves the omega-constant eigenvalue problem
% using the finite difference scheme from Ashok & Bhaduria, which solves for
% the horizontal eigenfunction (Ueig) at the cell centers of a non-uniform vertical
% grid.
%
% INPUT:
% rho  := mean density at grid centers [1xN]
% zf   := grid points at grid faces [1xN+1]
% f    := Coriolis frequency (constant)
% om   := anglular frequency (constant)
%
% OUTPUT
% S.L    := wavelength
% S.C    := phase speed (celerity)
% S.Cg   := group speed
% S.Ueig := horizontal eigenvalues (at cell centers 'zc')
%
% Created: Oladeji, S. June 2020
% Modified: Solano, M. July 2020

%om = 2*pi/(T*3600); % angular freq
%f = 2*7.29*(10^-5)*sin(lat*pi/180); % coriolis freq

dz = -diff(zf); % layer thickness vector (top to bottom)
zc = zf(1:end-1)/2 + zf(2:end)/2; % layer centers' depth vector
drho = rho(2:end) - rho(1:end-1);
dzc = zc(2:end) - zc(1:end-1);
N2 = -(9.81/rho0)*(drho./dzc);  % brunt vaisala frequency
N2 = [N2(1);N2;N2(end)];
N2 =(N2(1:end-1)/2 + N2(2:end)/2);
N2(N2<1e-10)=1e-10; % Mask negative and very small values**

% Compute parameters
H = nansum(dz);
N = numel(dz);
B = diag(repmat(-1,[1,N-1]));
B = padarray(B,[1 1],0);

% Main solver
A=[];
A(1,1) = 1;
A(1,2) = -1;
for i=2:N
    A(i,-1+i) =  2*1/(N2(i-1)*dz(i-1)*(dz(i-1)+dz(i)));
    A(i, 0+i) = -2*(N2(i-1)*dz(i-1) + N2(i)*dz(i))/(N2(i-1)*dz(i-1)*N2(i)*dz(i)*(dz(i-1)+dz(i)));
    A(i, 1+i) =  2*1/(N2(i)*dz(i)*(dz(i-1)+dz(i)));
end
i=i+1;
A(i,-1+i) = -1;
A(i, 0+i) = 1;

% Solve the EVP
[psi1,invCe2] = eig(A,B);
Ce2 = 1./(diag(invCe2));
Ce = sqrt(Ce2);

% k = sqrt(diag(k2));
[Ce,Is]=sort(Ce,'descend'); %sort
Ce(1) = []; Ce(end) = [];
psi = psi1(:,Is);
psi(:,1) = [];
psi(:,end) = [];

k = abs(sqrt((om^2-f^2))./Ce);
L  = 2*pi./k;
C  = om./k;
Cg = (om^2-f^2)./(om*k);

%% Get U structure
Ueig1 = psi(1:end-1,:)/2 + psi(2:end,:)/2;

%% Normalize U structure
dzu = repmat(dz,[1,size(dz,1)-1]);
AA = repmat(sum(Ueig1.^2.*dzu,1)./H,[N 1]).^(1/2);
AA(AA==0) = Inf;
Ueig2 = Ueig1./AA;
Ueig2(:,Ueig2(N,:)<0) = -Ueig2(:,Ueig2(N,:)<0);
Ueig = Ueig2;

%% Save output in a structure
S = struct('N2',N2,'k',k,'L',L,'C',C,'Cg',Cg,'Ueig',Ueig);
