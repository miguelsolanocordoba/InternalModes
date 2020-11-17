clc; clear; close all; 
%EIGEN_EXP solved the EVP for exponential stratification.  
% 
% Created: Miguel Solano, May 14, 2020.


%% Initialize: paths, constants and variables
set(groot,'defaultLineLineWidth',2);

% Paths 
addpath /data/msolano/matfiles
addpath /data/msolano/toolbox/chebfun
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/

figpath = '/data/msolano/figures/eigen_analytical/'; % figure bin

% Constants
g = 9.81;                       % gravity 
rho0 = 1025;                    % constant sea-water density  
omega = 12.1408331767/24/3600;  % M2 frequency

%% Constant stratification 
lat = 45;  % latitude 

% Initialize profiles using pre-defined function
% NOTE: rho and N2 are FUNCTION HANDLES (rho@(z))
[rhoFunc,N2Func,zIn] = ...
     InternalModes.StratificationProfileWithName('exponential');


% Variables
N0 = N2Func(0);   % Constant stratification  
n = 2*64;     % Number of modes and vertical grid points
z = linspace(zIn(1),zIn(2),n); % Output grid 
H = zIn(1); 

% Compute modes using Jeffrey's toolkit
im = InternalModes(rhoFunc,zIn,z,lat); %,...
       %   'nModes',n,'method','finiteDifference'); 

f =  im.f0;  %  Coriolis frequency 

% Plot density and stratification 
figure; 
plot(rhoFunc(z),z); 
xlabel('\rho'); ylabel('Depth[m]')
print('exp_rho.png','-r300','-dpng') 

figure; 
plot(N2Func(z),z); 
xlabel('N^2'); ylabel('Depth[m]')
print('exp_N2.png','-r300','-dpng') 


%% Jeffrey (InternalModes.m)
im.normalization = 'uMax'; 
[UeigJ,~,~,~] = im.ModesAtFrequency(omega); 

im.normalization = 'wMax'; 
[~,WeigJ,~,~] = im.ModesAtFrequency(omega); 

% Plot eigenvalues
figure 
subplot(121)
plot(UeigJ(:,1:5),z)
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]')

subplot(122)
plot(WeigJ(:,1:5),z)
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
print('exp_eigen_jeff.png','-r300','-dpng') 


%% Maarten (sturm_liouville_hyd_normalize.m) 
dz = z(2) -  z(1); 
Nb = flipud(sqrt(N2Func(z))');
[~,~,~,Weig,Ueig] = sturm_liouville_hyd_normalize(omega,Nb,dz,f);

% Plot eigenvalues
zf = z(1:end-1) + diff(z)/2;
N =  numel(zf); 
AA = repmat(-sum(Ueig.^2.*dz,1)./H,[N 1]).^(1/2); 
normU = max(max(abs(Ueig))); 
normW = max(max(abs(Weig))); 

figure 
subplot(121)
plot(-Ueig(:,2:6)./normU,flipud(zf'))
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]')

subplot(122)
plot(-Weig(:,2:6)./normW,flipud(zf'))
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
print('exp_eigen_maar.png','-r300','-dpng') 


%% Plot comparison: Jeffrey vs Maarten
figure 
subplot(121)
plot(UeigJ(:,1:5),z); hold on
plot(Ueig(:,2:6)./normU,flipud(zf'),'+k','MarkerSize',2)
title('Horizontal Eigenvalues')
xlabel('U_{eig}'); ylabel('Depth [m]')

subplot(122)
plot(-WeigJ(:,1:5),z); hold on
plot(Weig(:,2:6)./normW,flipud(zf'),'+k','MarkerSize',2)
title('Vertical Eigenvalues')
xlabel('W_{eig}'); ylabel('Depth [m]')
legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
print('exp_eigen_comp.png','-r300','-dpng')

% Move all figures to /data
system(['mv *.png ' figpath]);
