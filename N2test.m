clc; clear; close all; 
% Test specifying N2 for Early's InternalModes 
% 
% Here I'm checking that the stratification frequency (N2) given 
% by InternalModes is the same before and after masking negative 
% N2 values to 1e-10

%% Paths
addpath /data/msolano/Matlab     % All functions!
addpath /data/msolano/forOladeji
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
dirout = '/data/msolano/figures/modes_hycom/';

% Load profile (North Atlantic).
load(['profile_loc2.mat'])

% Dimensions and constants
[nz,nt] = size(profile.uiso); 

% Constants
g = 9.81;                        % gravity 
omega = 12.1408331767/24/3600;   % M2 frequency 
lat = profile.latitude;          % Latitude 
%rho0 = 1025;                     % constant sea-water density 

% Time Vector (September 1-15, 2016)
t = datenum(2016,9,1):datenum(0,0,0,1,0,0):datenum(2016,9,15); 
nt = numel(t); 

% Vertical grid (surface to bottom) 
zf = profile.zf_mean;  % cell faces 
zc = profile.zc_mean;  % cell centers
dz = diff(zf);         % layer thickness 
H = nansum(dz);        % total depth 
N = numel(dz);         % number of layers 


%% Density 
% Note: the density profile at this location is only stable at the 
% surface, i.e., diff(rho_mean)<0 between the first and second grid points. 
rho_mean = profile.rho_mean+1000;  % Mean sea-water density (sigma2 to rho) 
rhof_mean = interp1(zc,rho_mean,zf,'linear','extrap'); % rho at faces

% Initilized Early's InternalModes
im1 = InternalModes(rhof_mean,zf,zf,lat); %,'method','spectral');
im2 = InternalModes(rhof_mean,zf,zf,lat,'method','spectral');
im3 = InternalModes(rhof_mean,zf,zf,lat,'method','finiteDifference');

% Eigenfunctions
[~,Weig1,~,~] = im1.ModesAtFrequency(omega); 
Ueig1 = compute_ueig(Weig1,dz); 

[~,Weig2,~,~] = im1.ModesAtFrequency(omega); 
Ueig2 = compute_ueig(Weig2,dz); 

[~,Weig3,~,~] = im1.ModesAtFrequency(omega); 
Ueig3 = compute_ueig(Weig3,dz); 

rho0 = im1.rho0; 
f = im1.f0;

% N2 from Oladeji
S = compute_eigen(rho_mean,rho0,zf,f,omega)

% Mask N2 when negative to very small value (1e-10) 
%N2 = im.N2; 
%N2(N2<0)=1e-10;

[N2,rho_meanmod] = N2mask(rho_mean,rho0,im3.N2,diff(zc)); 
rhof_meanmod = interp1(zc,rho_meanmod,zf,'linear','extrap');
im4 = InternalModes(rhof_meanmod,zf,zf,lat,'method','finiteDifference');


%% Plots
% Stratification 
figure
subplot(121)
plot(rho_mean,zc,'k'); hold on 
plot(im1.rho,im1.z,'b');
plot(im2.rho,im2.z,'r');
plot(im3.rho,im3.z,'g');
plot(im4.rho,im4.z,'m');
title('Density'); xlabel('\rho [kg/m^3]'); ylabel('Depth [m]')
legend('HYCOM','WKB','Spectral','FD','FD (masked)','Location','SouthWest')
ylim([-300 0]);

subplot(122)
semilogx(S.N2,zf,'k'); hold on
semilogx(im1.N2,im1.z,'b');
semilogx(im2.N2,im2.z,'r');
semilogx(im3.N2,im3.z,'g');
semilogx(im4.N2,im4.z,'m');
title('Stratification'); xlabel('N^2'); ylabel('Depth [m]')
ylim([-300 0]);
legend('Reference','WKB','Spectral','FD','FD (masked)','Location','SouthWest')
print('stratification.png','-r300','-dpng')

% Eigenfunctions


system(['rm ' dirout '*.png']);
system(['mv *.png ' dirout]);

%% EoF
