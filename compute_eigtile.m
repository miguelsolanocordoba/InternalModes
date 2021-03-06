clc; clear; close all; 
%%COMPUTE_EIGTILE solves the EVP for a HYCOM tile
% COMPUTE_EIGTILE computes the eigenfunctions and eigenvalues in an 
% area NX x NY, from a HYCOM tile. The eigenfunctions/values are then
% saved in .mat format for later use. Currently 5 methods are used: 
% 
% Case 1: Spectral (Early) 
% Case 2: FD2 (Early) 
% Case 3: FD2 (Oladeji U) 
% Case 4: FD2 (Oladeji W) 
% Case 5: Unifrm (dz=25m)
% 
% Created July 13, 2020 by M. Solano

%% Defaults
set(0,'defaultAxesFontSize',6)
warning('off'); 

addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/
figpath = '/data/msolano/matfiles';

%% *** INPUT *** %% 
model = 'GLBc0.04';  % Global HYCOM
exptn = 190; % Experiment number
xtile = 15;  % blki 
ytile = 25;  % blkj
dzi = 25;    % layer thickness for uniform grid [m] ***
nmodes = 10; % number of modes saved 

% Convert to character strings
xtilestr = num2str(xtile); 
ytilestr = num2str(ytile); 
exptnstr = num2str(exptn);  

% Read hycom variables
hycom = read_hycom(model,exptn,xtile,ytile); % !!!edit read_hycom.m to change tile!!! 
[nx,ny,nz,nt] = size(hycom.rho); % tile size 
nfiles = nx*ny; 

% Constants
g = 9.806;                      % gravity 
omega = 12.1408331767/24/3600;  % M2 frequency

% Pre-allocation
%Ueig1 = zeros(nx,ny,41,nmodes); 
%Ueig2 = zeros(nx,ny,41,nmodes); 
%Ueig3 = zeros(nx,ny,41,nmodes); 
%Ueig4 = zeros(nx,ny,41,nmodes); 
Ueig = zeros(nx,ny,265,nmodes); 
Weig = zeros(nx,ny,265,nmodes); 

% To calculate the number of layers for the uniform case:
% nzM*dzi > max(depth)
% nzM > max(depth)/dzi ~ 265 for amazon (dzi=25m; max(h)=6400m)

%k1 = zeros(nx,ny,nmodes); 
%k2 = zeros(nx,ny,nmodes); 
%k3 = zeros(nx,ny,nmodes); 
%k4 = zeros(nx,ny,nmodes); 
k = zeros(nx,ny,nmodes); 


fprintf('\nComputing eigenvalues...\n')
count = 0; 
counti = 0;
for ii = 1:nx
    counti = counti + 1; 
    countj = 0; 
    for j = 1:ny
        % Counters
	countj = countj + 1; 
	count = count + 1; 
	fprintf('%d/%d\n',count,nfiles)

	% Skip if on land
	if hycom.h(ii,j)>1e4 
	   continue  % ignore masked points
        elseif hycom.h(ii,j)<dzi
	   continue  % ignore very shallow zones (<dzi)
	end

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
        zf_mean = mean(zf,2); 
	zc_mean = mean(zc,2); 
        dz = diff(zf_mean); 
        H = nansum(dz); % total depth  
        N = numel(dz);  % # of layers 
        
        % Interpolate density to time-mean layer and compute mean
        rho = flipud(squeeze(hycom.rho(ii,j,1:nz,:))); 
	rhoint = zeros(nz,nt);
	for i = 1:nt
	    rhoint(:,i) = interp1(zc(:,i),rho(:,i),zc_mean,'linear','extrap'); 
        end
        %rho_mean = sort(mean(rhoint,2),'descend'); % mean sea-water density
        rho_meant = mean(rhoint,2)+1000; % mean sea-water density
	rhof_meant = interp1(zc_mean,rho_meant,zf_mean,'linear','extrap');

	% Mask rho
        im = InternalModes(rhof_meant,zf_mean,zf_mean,hycom.lat(ii,j),...
	'method','finiteDifference');
	rho0 = im.rho0;
	[N2,rho_mean] = N2mask(rho_meant,im.rho0,im.N2,diff(zc_mean));
%	rhof_mean = interp1(zc_mean,rho_mean,zf_mean,'linear','extrap');

	% Compute perturbation pressure
%        pert = compute_pertpress(rhoint,zc,zf); % perturbation pressure

%*****  Compute eigenfunctions 
%% Case 1: Spectral (Early) 
% Input/output at faces for Weig. Ueig computed at centers. 
%        im = InternalModes(rhof_mean,zf_mean,zf_mean,hycom.lat(ii,j),...
%	                   'method','spectral');
%	rho0 = im.rho0; 
%        [~,Weig1,~,k] = im.ModesAtFrequency(omega);
%        Ueig1(ii,j,1:nz,:) = compute_ueig(Weig1(:,1:nmodes),dz);
%	k1(ii,j,:) = k(:,1:nmodes); clear k
%
%%% Case 2: FD (Early) 
%% Input/output at faces for Weig. Ueig computed at centers. 
%        imFD = InternalModes(rhof_mean,zf_mean,zf_mean,hycom.lat(ii,j),...
%	                     'method','finiteDifference',...
%			     'orderOfAccuracy',2);
%	[~,Weig2,~,k] = imFD.ModesAtFrequency(omega); 
%        Ueig2(ii,j,1:nz,:) = compute_ueig(Weig2(:,1:nmodes),dz); 
%	k2(ii,j,:) = k(:,1:nmodes); clear k
%        %%% Case 3: FD (Oladeji W)
%% See compute_eigen for details
%        S3 = compute_eigen(rho_mean,rho0,zf_mean,f,omega);
%	Ueig3(ii,j,1:nz,:) = S3.Ueig(:,1:nmodes); 
%	k3(ii,j,:) = S3.k(1:nmodes);
%        	
%%% Case 4: FD (Oladeji U)
%% See compute_eigenU for details
        S4 = compute_eigenU(rho_mean,rho0,zf_mean,f,omega);
%	Ueig4(ii,j,1:nz,:) = S4.Ueig(:,1:nmodes); 
%	k4(ii,j,:) = S4.k(1:nmodes);

%% Case 5: Uniform (Maarten)
% Note: Uniform case uses stratification from Oladeji's function (already masked). 
	nzM = round(H/dzi); % number of layers
	nzM = max([nzM 20]); % Require minimum number points (10)
        zfM = linspace(-H,0,nzM+1)';
        dz1 = diff(zfM);
        zcM = zfM(1:end-1) + dz1/2;
        N2M = interp1(zc_mean,S4.N2,zfM,'linear','extrap');
        [~,~,L5,Weig6t,Ueig6t] = sturm_liouville_hyd_normalize(omega,sqrt(N2M),dz1(1),f);
        Ueig6t(:,2:2:end) = -Ueig6t(:,2:2:end);

        Weig(ii,j,1:nzM,:) = Weig6t(:,2:nmodes+1);
        Ueig(ii,j,1:nzM,:) = Ueig6t(:,2:nmodes+1);
	k(ii,j,:) = 2*pi./(L5(1:nmodes));

    end
end

lon = hycom.lon(1:nx,1:ny); 
lat = hycom.lat(1:nx,1:ny); 
depth = hycom.h(1:nx,1:ny); 

%% Save output 
save([figpath '/' model '_' exptnstr '_' xtilestr '_' ytilestr '_eig.mat'],'Weig','Ueig','k')

%save([figpath '/eigentile_AMZN1942.mat'],'lon','lat','depth','Ueig1','Ueig2',...
%     'Ueig3','Ueig4','Ueig5','k1','k2','k3','k4','k5')

