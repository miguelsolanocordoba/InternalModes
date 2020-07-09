clc; clear; close all; 
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
%zc = zf(1:end-1)/2 + zf(2:end)/2;
nz = numel(zc); 
H = zIn(1); 

% Compute modes using Jeffrey's toolkit
imf = InternalModes(rhoFunc,zIn,zf,lat,'normalization','wMax'); %,...
imc = InternalModes(rhoFunc,zIn,zc,lat,'normalization','wMax'); %,...
       %   'nModes',n,'method','finiteDifference'); 

f =  imf.f0;  %  Coriolis frequency 

% Plot density 
plotf = input('Plot density? Y/N [N]','s'); 
if plotf=='Y'
  figure
  plot(rhoFunc(zIn),zIn); 
  xlabel('\rho'); ylabel('Depth[m]')
  text(1035,-1000,['N_0 = ' num2str(N0)])
  print('const_rho.png','-r300','-dpng') 
end

%% Maarten (sturm_liouville_hyd_normalize.m) 
%dz = z(2) -  z(1); 
%Nb = N0*ones(size(z));
%[~,~,~,Weig,Ueig] = sturm_liouville_hyd_normalize(omega,Nb,dz,f);

N2 = N2Func(zc); 
rho = rhoFunc(zc); 
[~,~,~,Weig,Ueig] = compute_eigen(rho',zf',f,omega); 
%[Weig,Ueig] = comp_eigen_ola(dz,N2,f,omega); 

% Plot eigenvalues
zc2 = zf(1:end-1) + diff(zf)/2;
%N =  numel(zf); 
%AA = repmat(-sum(Ueig.^2.*dz,1)./H,[N 1]).^(1/2); 
%normU = max(max(abs(Ueig))); 
%normW = max(max(abs(Weig))); 

plotf = input('Plot Olas modes? Y/N [N]','s'); 
if plotf=='Y'
   figure 
   subplot(121)
   plot(Ueig(:,1:5),zc)
   %plot(Ueig(:,2:6)./normU,flipud(zf'))
   title('Horizontal Eigenvalues')
   xlabel('U_{eig}'); ylabel('Depth [m]')
   
   subplot(122)
   plot(Weig(:,1:5),zf)
   title('Vertical Eigenvalues')
   xlabel('W_{eig}'); ylabel('Depth [m]')
   legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
   print('const_eigen_ola.png','-r300','-dpng') 
end

%% Jeffrey (InternalModes.m)
%im.normalization = 'uMax'; 
%[UeigJ,~,~,~] = im.ModesAtFrequency(omega); 

im.normalization = 'wMax'; 
[~,WeigJ,~,~] = imf.ModesAtFrequency(omega); 
UeigJ = compute_ueig(Weig,dz'); 


% Plot eigenvalues
plotf = input('Plot Jeffs modes? Y/N [N]','s'); 
if plotf=='Y'
   figure 
   subplot(121)
   plot(UeigJ(:,1:5),zc)
   title('Horizontal Eigenvalues')
   xlabel('U_{eig}'); ylabel('Depth [m]')
   legend('n=1','n=2','n=3','n=4','n=5','Location','NorthWest')
   
   subplot(122)
   plot(WeigJ(:,1:5),zf)
   title('Vertical Eigenvalues')
   xlabel('W_{eig}'); ylabel('Depth [m]')
   legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
   print('const_eigen_jeff.png','-r300','-dpng') 
end


%% Plot comparison: Jeffrey vs Maarten
plotf = input('Plot compare modes? Y/N [N]','s'); 
if plotf=='Y'
   figure 
   subplot(121)
   plot(UeigJ(:,1:5),zc); hold on
   %plot(Ueig(:,2:6)./normU,flipud(zf'),'+k','MarkerSize',2)
   plot(Ueig(:,1:5),zc,'+k','MarkerSize',2)
   title('Horizontal Eigenvalues')
   xlabel('U_{eig}'); ylabel('Depth [m]')
   
   subplot(122)
   plot(WeigJ(:,1:5),zf); hold on
   %plot(Weig(:,2:6)./normW,flipud(zf'),'+k','MarkerSize',2)
   plot(Weig(:,1:5),zf,'+k','MarkerSize',2)
   title('Vertical Eigenvalues')
   xlabel('W_{eig}'); ylabel('Depth [m]')
   legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
   print('const_eigen_comp.png','-r300','-dpng')
end

% Move all figures to /data
system(['mv *.png ' figpath]);
