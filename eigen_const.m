clc; clear; close all; warning off 
%EIGEN_CONST solves the EVP for constant stratification.  
% EIGEN_CONST is solves the the omega-constant vertical mode
% eigenvalue problem: 
%
% PDE: W'' + N^2/c_n^2 W = 0
% BCs: W(0) = 0; W(H) = 0
% 
% Where the prime (') denotes derivative with respect to coordinate 
% (z), which is z=0 at the surface and z=H at the bottom. N is the  
% constant stratification frequency and c_n is the eigenspeed. 
% 
% c_n = sqrt(omega^2 - f^2)/k_n 
% N = sqrt(-g/rho_0 rho'); 
% 
% where (rho) is the density, g is gravity and k_n is the wavenumber. 
% Wave properties are computed as follows: 
% 
% L = 2*pi/k_n         % wave-length
% C = omega/k_n        % phase speed
% Cg= c_n^2*k_n/omega  % group speed
%
% Output: 
% const_rho.png           % plots density/stratification 
% comp_eigen_ola.png      % compares W from Oladeji-analytical 
% comp_eigen_maar.png     % compares W from Maarten-analytical 
% comp_eigen_jeffrey.png  % compares W from Jeffrey-analytical
% k.png                   % plots O-M-J-analytical wavenumber
% L.png                   % plots O-M-J-analytical wavelength
% C.png                   % plots O-M-J-analytical phase speed
% Cg.png                  % plots O-M-J-analytical group speed
%
% Created: Miguel Solano, May 14, 2020.

%% Defaults
set(groot,'defaultLineLineWidth',2);

% Paths 
addpath /data/msolano/Matlab/ocean_physics
addpath /data/msolano/matfiles
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/
figpath = '/data/msolano/figures/modes_analytical/'; % figure bin

%% Constants
g = 9.81;                       % gravity 
rho0 = 1025;                    % constant sea-water density  
omega = 12.1408331767/24/3600;  % M2 frequency
lat = 45;                       % latitude 


% Initialize profiles using pre-defined function
% NOTE: rho and N2 are FUNCTION HANDLES (rho@(z))
[rhoFunc,N2Func,zIn] = ...
     InternalModes.StratificationProfileWithName('constant');

% Variables
N0 = N2Func(0);   % Constant stratification NOTE: N0 = N^2 !!!
n = 2*64;     % Number of modes and vertical grid points
zf = linspace(zIn(1),zIn(2),n); % Output grid 
dz = diff(zf); 
zc = zf(1:end-1) + dz/2;
nz = numel(zc); 
H = nansum(dz);

rho = rhoFunc(zf); 

% Compute modes using Jeffrey's toolkit
imf = InternalModes(rho,zf,zf,lat,'method','finiteDifference'); %,...
imc = InternalModes(rho,zf,zc,lat,'method','finiteDifference'); %,...

f =  imf.f0;  %  Coriolis frequency 

% Plot density 
figure
plot(rhoFunc(zIn),zIn); 
xlabel('\rho'); ylabel('Depth[m]')
text(1035,-1000,['N_0 = ' num2str(N0)])
print('const_rho.png','-r300','-dpng') 

%% EVP: Analytical solution to non-hydrostatic with uniform stratification
% Eigenfunctions
% nmodes := # of modes solved
nmodes = nz + 1; % number of modes 
Weig = zeros(nz+1,nmodes); % Weig(1)=0 [bottom BC]; Weig(end)=0 [surface BC]
Ueig = zeros(nz,nz); 

% Weig = sin(pi*z*nn/H) for nn=1,2,3,..,nmodes
% Ueig = d(Weig)/dz (Horizontal eigenfunction at cell centers)
nn = [1:nz+1]; % modes nn = 1,2,3,...,nmodes 
Weig(2:end-1,:) = sin(pi*zf(2:end-1)'*nn/H); 
Ueig = compute_ueig(Weig,dz'); % see compute_ueig.m for details

% Eigenvalues
alpha = (omega^2-f^2)/(N0-omega^2);
lambda = pi*nn/H; 

% Solutions from Gerkema & Zimmerman lecture notes
kn = lambda*sqrt(alpha);                % wave-number 
L = 2*pi./kn;                           % wave-length 
C = (H*omega./(nn*pi)) * sqrt(1/alpha); % phase-speed
%Cg = (H./(nn*pi)) * (1/(omega*(N0-f^2))) * sqrt(omega^2-f^2) * (N0-omega^2)^(3/2);
Cg = (H./(nn*pi)) * (1/(omega*(N0-f^2))) * sqrt(omega^2-f^2) * (N0)^(3/2);


%% Maarten (sturm_liouville_hyd_normalize.m) ; 
N2M = repmat(sqrt(N0),1,nz+1); 
[CM,CgM,LM,WeigM,UeigM] = sturm_liouville_hyd_normalize(omega,N2M,dz(1),f);
WeigM = WeigM/(max(max(WeigM)));
kM = 2*pi./LM; 

%% Oladeji
N2 = N2Func(zc); 
rho = rhoFunc(zc); 
S = compute_eigen(rho',zf',f,omega); 
kO = S.k; LO = S.L; CO = S.C; CgO = S.Cg; WeigO = S.Weig; UeigO = S.Ueig;
%[kO,CO,CgO,LO,WeigO,UeigO] = compute_eigen(rho',zf',f,omega); 

%% Jeffrey (InternalModes.m)
%imc.normalization = 'uMax'; 
%[UeigJ,~,~,~] = imc.ModesAtFrequency(omega); 

imf.normalization = 'wMax'; 
[~,WeigJ,hJ,kJ] = imf.ModesAtFrequency(omega); 
UeigJ = compute_ueig(WeigJ,dz'); 
umax = max(UeigJ(2,:)); 

LJ = 2*pi./kJ; 
CJ = omega./kn; 
CgJ = (omega^2-f^2)./(omega*kJ); 
%CeJ = sqrt(g*hJ); % eigen-speed 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmodes = 5;    % number of modes to plot
n = [1:nmodes]; 

%% EIGENFUNCTIONS
% Oladeji
figure 
subplot(121)
plot(UeigO(:,n),zc); hold on
plot(Ueig(:,n),zc,'+k','MarkerSize',2);
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]'); xlim([-umax umax])

subplot(122)
plot(WeigO(:,n),zf); hold on
plot(Weig(:,n),zf,'+k','MarkerSize',2);
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
print('const_eigen_ola.png','-r300','-dpng') 


% Maarten 
figure 
subplot(121)
plot(UeigM(:,n+1),fliplr(zc)); hold on
plot(Ueig(:,n),zc,'+k','MarkerSize',2);
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]'); xlim([-umax umax]) 
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
   
subplot(122)
plot(-WeigM(:,n+1),fliplr(zc)); hold on
plot(Weig(:,n),zf,'+k','MarkerSize',2);
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
print('const_eigen_maar.png','-r300','-dpng') 


% Jeffrey Early  
figure 
subplot(121)
plot(UeigJ(:,n),zc); hold on
plot(Ueig(:,n),zc,'+k','MarkerSize',2);
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]'); xlim([-umax umax])
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
   
subplot(122)
plot(WeigJ(:,n),zf); hold on
plot(Weig(:,n),zf,'+k','MarkerSize',2);
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
print('const_eigen_jeff.png','-r300','-dpng') 


%% EIGENVALUES
% Compare wave-lengths 
figure
plot(n,kO(n),'b'); hold on
plot(n,kM(n),'r'); 
plot(n,kJ(n),'g'); 
plot(n,kn(n),'*k'); 
title('Wavenumber (k)') 
ylabel('[rad/m]'); xlabel('Mode'); 
xticklabels({'1',' ','2',' ','3',' ','4',' ','5'})
legend('Oladeji','Maarten','Jeffrey','Analytical','Location','NorthWest')
print('k.png','-r300','-dpng')

figure
plot(n,LO(n)/1000,'b'); hold on
plot(n,LM(n)/1000,'r'); 
plot(n,LJ(n)/1000,'g'); 
plot(n,L(n)/1000,'*k'); 
title('Wavelength (L)') 
ylabel('[km]'); xlabel('Mode');
xticklabels({'1',' ','2',' ','3',' ','4',' ','5'})
legend('Oladeji','Maarten','Jeffrey','Analytical')
print('L.png','-r300','-dpng')

figure
plot(n,CO(n),'b'); hold on
plot(n,CM(n),'r'); 
plot(n,CJ(n),'g'); 
plot(n,C(n),'*k'); 
title('Phase-speed (C)') 
ylabel('[m/s]'); xlabel('Mode');
xticklabels({'1',' ','2',' ','3',' ','4',' ','5'})
legend('Oladeji','Maarten','Jeffrey','Analytical')
print('C.png','-r300','-dpng')

figure
plot(n,CgO(n),'b'); hold on
plot(n,CgM(n),'r'); 
plot(n,CgJ(n),'g'); 
plot(n,Cg(n),'*k'); 
title('Group-speed (C_g)') 
ylabel('[m/s]'); xlabel('Mode'); 
xticklabels({'1',' ','2',' ','3',' ','4',' ','5'})
legend('Oladeji','Maarten','Jeffrey','Analytical')
print('Cg.png','-r300','-dpng')

fprintf('\nSolutions:') 
fprintf('\nMode 1 wavelength (L)')
fprintf('\nAnalytical = %4.2f[km]',L(1)/1000)
fprintf('\nOladeji = %4.2f[km]',LO(1)/1000)
fprintf('\nMaarten = %4.2f[km]',LM(1)/1000)
fprintf('\nJeffrey = %4.2f[km]\n',LJ(1)/1000)
fprintf('\nMode 1 group-speed (Cg)')
%fprintf('\n(Analytical = %4.4f',CgO(1))
fprintf('\nOladeji = %4.4f[m/s]',CgO(1))
fprintf('\nMaarten = %4.4f[m/s]',CgM(1))
fprintf('\nJeffrey = %4.4f[m/s]\n',CgJ(1))

fprintf('\nMode 1 phase-speed (C)')
fprintf('\nAnalytical = %4.4f[m/s]',C(1))
fprintf('\nOladeji = %4.4f[m/s]',CO(1))
fprintf('\nMaarten = %4.4f[m/s]',CM(1))
fprintf('\nJeffrey = %4.4f[m/s]\n',CJ(1))

% Move all figures to /data
system(['mv *.png ' figpath]); 
