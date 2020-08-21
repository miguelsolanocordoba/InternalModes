clc; clear; close all; 
%%COMPUTE_EIGTILE solves the EVP for a HYCOM tile
% COMPUTE_EIGTILE computes the eigenfunctions and eigenvalues in an 
% area NX x NY, from a HYCOM tile. The eigenfunctions/values are then
% saved in .mat format for later use. Currently 5 methods are used: 
% 
% Case 1: WKB (Early) 
% Case 2: FD2 (Early) 
% Case 3: FD2 (Oladeji U) 
% Case 4: FD2 (Oladeji W) 
% Case 5: Unifrm (dz=25m)
% 
% Created July 13, 2020 by M. Solano

%% Defaults
set(0,'defaultAxesFontSize',6)
warning('off'); 

addpath /data/msolano/Matlab/ocean_physics
addpath /data/msolano/toolbox/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/
figpath = '/data/msolano/figures/test';

%% Read HYCOM variables
hycom = read_hycom(); % !!!edit read_hycom.m to change tile!!! 
[nx,ny,nz,nt] = size(hycom.rho); % tile size 
nx = 100; % over-ride x-size 
ny = 100; % over-ride y-size
dzi = 25; % layer thickness for uniform grid [m]

% Constants
g = 9.806;                      % gravity 
omega = 12.1408331767/24/3600;  % M2 frequency
nfiles = nx*ny; 
nmodes = 5; 

% Variables to plot
Ueig1 = zeros(nx,ny,41,nmodes); 
Ueig2 = zeros(nx,ny,41,nmodes); 
Ueig3 = zeros(nx,ny,41,nmodes); 
Ueig4 = zeros(nx,ny,41,nmodes); 
Ueig5 = zeros(nx,ny,200,nmodes); 

k1 = zeros(nx,ny,nmodes); 
k2 = zeros(nx,ny,nmodes); 
k3 = zeros(nx,ny,nmodes); 
k4 = zeros(nx,ny,nmodes); 
k5 = zeros(nx,ny,nmodes); 


fprintf('\nComputing eigenvalues...\n')
count = 0; 
counti = 0;
for ii = 1:nx
    counti = counti + 1; 
    countj = 0; 
    for j = 1:ny
	countj = countj + 1; 
	count = count + 1; 
	fprintf('%d/%d\n',count,nfiles)
	   
        f = coriolis(hycom.lat(ii,j)); % Coriolis frequency
    
        % Compute grid points from layer thickness
        [~,~,nz,nt] = size(hycom.rho); 
        zc1 = zeros(nz,nt); zf1 = zeros(nz+1,nt); nzm = zeros(nt,1); 
        for l = 1:nt
            zf1(2:end,l) = squeeze(cumsum(hycom.dz(ii,j,:,l),3)); 
            zc1(:,l) = 0.5.*(zf1(1:nz,l) + zf1(2:nz+1,l));
	    thk0 = diff(zf1(:,l))>0.01; % number of layers <1cm
	    nzm(l) = sum(thk0);
        end
        nz = min(nzm);  % number of layers (nz) cut to minimum # of non-zero layers
        zf = -flipud(zf1(1:nz+1,:)); % flip and negative down (convention)
        zc = -flipud(zc1(1:nz,:)); 
        
        % Compute mean layer	
        zf_mean = mean(zf,2); zc_mean = mean(zc,2); 
        dz = diff(zf_mean); 
        H = nansum(dz); % total depth  
        N = numel(dz);  % # of layers 
        
        % Compute mean density 
        rho = flipud(squeeze(hycom.rho(ii,j,1:nz,:))); 
        rho_mean = flipud(squeeze(mean(hycom.rho(ii,j,1:nz,:),4)));
	rhof_mean = zeros(nz+1,1); 
	for l = 1:nt
	    rhof = interp1(zc(:,l),rho(:,l),zf(:,l),'linear','extrap');
	end
	rhof_mean = mean(rhof,2);  % mean density at faces 

%*****  Compute eigenfunctions 
%% Case 1: WKB (Early) 
        im = InternalModes(rhof_mean,zf_mean,zf_mean,hycom.lat(ii,j));
        [~,Weig1,~,k] = im.ModesAtFrequency(omega);
        Ueig1(ii,j,1:nz,:) = compute_ueig(Weig1(:,1:nmodes),dz);
	k1 = k(:,1:nmodes); clear k

%% Case 2: FD (Early) 
        imFD = InternalModes(rhof_mean,zf_mean,zf_mean,hycom.lat(ii,j),...
	                     'method','finiteDifference',...
			     'orderOfAccuracy',2);
	[~,Weig2,~,k] = imFD.ModesAtFrequency(omega); 
        Ueig2(ii,j,1:nz,:) = compute_ueig(Weig2(:,1:nmodes),dz); 
	k2 = k(:,1:nmodes); clear k
	
%% Case 3: FD (Oladeji W)
        S3 = compute_eigen(rho_mean,zf_mean,f,omega);
	Ueig3(ii,j,1:nz,:) = S3.Ueig(:,1:nmodes); 
	k3(ii,j,:) = S3.k(1:nmodes);
        	
%% Case 4: FD (Oladeji U)
        S4 = compute_eigenU(rho_mean,zf_mean,f,omega);
	Ueig4(ii,j,1:nz,:) = S4.Ueig(:,1:nmodes); 
	k4(ii,j,:) = S4.k(1:nmodes);

%% Case 5: Uniform (Maarten)
	nzM = round(H/dzi); % number of layers
        zfM = linspace(-H,0,nzM+1)';
        dz1 = diff(zfM);
        zcM = zfM(1:end-1) + dz1/2;
        N2M = interp1(zf_mean,im.N2,zfM,'linear','extrap');
	N2M(N2M<0) = 1e-10; % Mask unstable (negative) stratification
        [~,~,L5,~,Ueig6t] = sturm_liouville_hyd_normalize(omega,sqrt(N2M),dz1(1),f);
        Ueig6t(:,2:2:end) = -Ueig6t(:,2:2:end);

        Ueig5(ii,j,1:nzM,:) = Ueig6t(:,2:nmodes+1);
	k5(ii,j,:) = 2*pi./(L5(2:nmodes+1));
    end
end

lon = hycom.lon(1:nx,1:ny); 
lat = hycom.lat(1:nx,1:ny); 
depth = hycom.h(1:nx,1:ny); 

%% Save output 
save('eigentile.mat','lon','lat','depth','Ueig1','Ueig2','Ueig3','Ueig4',...
                     'Ueig5','k1','k2','k3','k4','k5')
system(['mv eigentile.mat ' figpath '/']);

