clc; clear; close all; 
%EIGEN_CONTS solves the EVP for constant stratification.  
% 
% Created: Miguel Solano, May 14, 2020.


%% Initialize: paths, constants and variables
set(groot,'defaultLineLineWidth',2);

% Paths 
addpath /data/msolano/matfiles
addpath /data/msolano/toolbox/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/

figpath = '/data/msolano/figures/modes_analytical/'; % figure bin

% Constants
g = 9.81;                       % gravity 
rho0 = 1025;                    % constant sea-water density  
omega = 12.1408331767/24/3600;  % M2 frequency

%% Constant stratification 
lat = 45;  % latitude 

% Initialize profiles using pre-defined function
% NOTE: rho and N2 are FUNCTION HANDLES (rho@(z))
[rhoFunc,N2Func,zIn] = ...
     InternalModes.StratificationProfileWithName('constant');

% Variables
N0 = N2Func(0);   % Constant stratification  
n = 2*64;     % Number of modes and vertical grid points
z = linspace(zIn(1),zIn(2),n); % Output grid 
H = zIn(1); 

% Compute modes using Jeffrey's toolkit
im = InternalModes(rhoFunc,zIn,z,lat); %,...
       %   'nModes',n,'method','finiteDifference'); 

f =  im.f0;  %  Coriolis frequency 

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

N2 = N2Func(z); 
dz = diff(z)'; 
[Weig,Ueig] = comp_eigen_ola(dz,N2,f,omega); 

% Plot eigenvalues
zc = z(1:end-1) + diff(z)/2;
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
   plot(Weig(:,1:5),z)
   title('Vertical Eigenvalues')
   xlabel('W_{eig}'); ylabel('Depth [m]')
   legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
   print('const_eigen_ola.png','-r300','-dpng') 
end

%% Jeffrey (InternalModes.m)
%im.normalization = 'uMax'; 
%[UeigJ,~,~,~] = im.ModesAtFrequency(omega); 

im.normalization = 'wMax'; 
[~,WeigJ,~,~] = im.ModesAtFrequency(omega); 

H = nansum(dz); N = numel(dz); 
%dW2 = WeigJ(2:end,:) - WeigJ(1:end-1,:); 
dW2 = diff(WeigJ); 
dzu = repmat(dz,[1,size(dz,1)+1]); 
Ueig1 = dW2./dzu; 
AA = repmat(sum(Ueig1.^2.*dzu,1)./H,[N 1]).^(1/2); 
AA(AA==0) = Inf; 
UeigJ = Ueig1./AA; 
UeigJ(:,UeigJ(N,:)<0) = -UeigJ(:,UeigJ(N,:)<0);  

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
   plot(WeigJ(:,1:5),z)
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
   plot(WeigJ(:,1:5),z); hold on
   %plot(Weig(:,2:6)./normW,flipud(zf'),'+k','MarkerSize',2)
   plot(Weig(:,1:5),z,'+k','MarkerSize',2)
   title('Vertical Eigenvalues')
   xlabel('W_{eig}'); ylabel('Depth [m]')
   legend('n=1','n=2','n=3','n=4','n=5','Location','NorthEast')
   print('const_eigen_comp.png','-r300','-dpng')
end

% Move all figures to /data
system(['mv *.png ' figpath]);
