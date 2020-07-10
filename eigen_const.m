clc; clear; close all; warning off 
%EIGEN_CONST solves the EVP for constant stratification.  
% EIGEN_CONST is an interactive script that solved and 
% plots the solution to the Sturm-Liouville 
%
% Created: Miguel Solano, May 14, 2020.

%% Defaults
set(groot,'defaultLineLineWidth',2);

% Paths 
addpath /data/msolano/Matlab/ocean_physics
addpath /data/msolano/matfiles
addpath /data/msolano/toolbox/InternalModes
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
N0 = N2Func(0);   % Constant stratification  
n = 2*64;     % Number of modes and vertical grid points
zf = linspace(zIn(1),zIn(2),n); % Output grid 
dz = diff(zf); 
zc = zf(1:end-1) + dz/2;
nz = numel(zc); 
H = nansum(dz);

% Compute modes using Jeffrey's toolkit
imf = InternalModes(rhoFunc,zIn,zf,lat); %,...
imc = InternalModes(rhoFunc,zIn,zc,lat); %,...

f =  imf.f0;  %  Coriolis frequency 

% Plot density 
figure
plot(rhoFunc(zIn),zIn); 
xlabel('\rho'); ylabel('Depth[m]')
text(1035,-1000,['N_0 = ' num2str(N0)])
print('const_rho.png','-r300','-dpng') 

%% Compute Eigenfunctions
% Analytical solutions 
Weig = zeros(nz+1,nz+1); 
Ueig = zeros(nz,nz); 
for i = 1:nz  % loop over nmodes=nz;
    Weig(:,i) = sin(zf*pi*i/H); 
end
Ueig = compute_ueig(Weig,dz');


%% Maarten (sturm_liouville_hyd_normalize.m) ; 
[~,~,~,WeigM,UeigM] = sturm_liouville_hyd_normalize(omega,N0*ones(size(zf)),dz(1),f);
WeigM = WeigM/(max(max(WeigM)));

%% Oladeji
N2 = N2Func(zc); 
rho = rhoFunc(zc); 
[~,~,~,WeigO,UeigO] = compute_eigen(rho',zf',f,omega); 

%% Jeffrey (InternalModes.m)
%imc.normalization = 'uMax'; 
%[UeigJ,~,~,~] = imc.ModesAtFrequency(omega); 

imf.normalization = 'wMax'; 
[~,WeigJ,~,~] = imf.ModesAtFrequency(omega); 
UeigJ = compute_ueig(WeigJ,dz'); 
umax = max(UeigJ(2,:)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Oladeji
figure 
subplot(121)
plot(UeigO(:,1:5),zc); hold on
plot(Ueig(:,1:5),zc,'+k','MarkerSize',2);
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]'); xlim([-umax umax])

subplot(122)
plot(WeigO(:,1:5),zf); hold on
plot(Weig(:,1:5),zf,'+k','MarkerSize',2);
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
print('const_eigen_ola.png','-r300','-dpng') 


% Maarten 
figure 
subplot(121)
plot(UeigM(:,2:6),fliplr(zc)); hold on
plot(Ueig(:,1:5),zc,'+k','MarkerSize',2);
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]'); xlim([-umax umax]) 
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
   
subplot(122)
plot(-WeigM(:,2:6),fliplr(zc)); hold on
plot(Weig(:,1:5),zf,'+k','MarkerSize',2);
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
print('const_eigen_maar.png','-r300','-dpng') 


% Jeffrey Early  
figure 
subplot(121)
plot(UeigJ(:,1:5),zc); hold on
plot(Ueig(:,1:5),zc,'+k','MarkerSize',2);
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]'); xlim([-umax umax])
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
   
subplot(122)
plot(WeigJ(:,1:5),zf); hold on
plot(Weig(:,1:5),zf,'+k','MarkerSize',2);
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
%legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
print('const_eigen_jeff.png','-r300','-dpng') 


% Move all figures to /data
system(['mv *.png ' figpath]); 
