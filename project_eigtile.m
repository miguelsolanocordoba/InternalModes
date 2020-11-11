clc; clear; close all; 
%%PROJECT_EIGTILE projects eigenfunctions onto velocity/pressure
% Used to project the saved eigenfunctions from compute_eigtile 
% into velocity/pressure. Makes plots. 
%
% Currently 5 methods are used: 
% 
% Case 1: WKB (Early) 
% Case 2: FD2 (Early) 
% Case 3: FD2 (Oladeji U) 
% Case 4: FD2 (Oladeji W) 
% Case 5: Unifrm (dz=25m)
% 
% Created August 21, 2020 by M. Solano

%% Defaults
set(0,'defaultAxesFontSize',6)
warning('off'); 

addpath /data/msolano/Matlab/ocean_physics
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/
datapath = '/data/msolano/matfiles';

%% *** INPUT *** %% 
model = 'GLBc0.04';  % Global HYCOM
exptn = 221; % Experiment number
xtile = 15;  % blki
ytile = 25;  % blkj 
dzi = 25;    % layer thickness for uniform grid [m] ***
nmodes = 10; % number of modes saved 

% Convert to character strings
xtilestr = num2str(xtile);
ytilestr = num2str(ytile);
exptnstr = num2str(exptn);

% Read hycom variables
hycom = read_hycom(exptn,xtile,ytile); % !!!edit read_hycom.m to change tile!!! 
[nx,ny,nz,nt] = size(hycom.rho); % tile size 
nfiles = nx*ny;

% Constants
g = 9.806;                      % gravity 
omega = 12.1408331767/24/3600;  % M2 frequency

% Pre-allocation
r2u = zeros(nx,ny,nmodes); r2v = zeros(nx,ny,nmodes); r2p = zeros(nx,ny,nmodes); 
%r2u1 = zeros(nx,ny,nmodes); r2v1 = zeros(nx,ny,nmodes); r2p1 = zeros(nx,ny,nmodes); 
%r2u2 = zeros(nx,ny,nmodes); r2v2 = zeros(nx,ny,nmodes); r2p2 = zeros(nx,ny,nmodes); 
%r2u3 = zeros(nx,ny,nmodes); r2v3 = zeros(nx,ny,nmodes); r2p3 = zeros(nx,ny,nmodes); 
%r2u4 = zeros(nx,ny,nmodes); r2v4 = zeros(nx,ny,nmodes); r2p4 = zeros(nx,ny,nmodes); 
%r2u5 = zeros(nx,ny,nmodes); r2v5 = zeros(nx,ny,nmodes); r2p5 = zeros(nx,ny,nmodes); 

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
	   
        if hycom.h(ii,j)>1e4
           continue  % ignore masked points
        elseif hycom.h(ii,j)<25
           continue  % ignore very shallow zones
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

        % flip variables upside-down 
        uiso = flipud(squeeze(hycom.uiso(ii,j,1:nz,:))); 
        viso = flipud(squeeze(hycom.viso(ii,j,1:nz,:))); 
        rho = flipud(squeeze(hycom.rho(ii,j,1:nz,:))); 
        rho_mean = flipud(squeeze(mean(hycom.rho(ii,j,1:nz,:),4)));
	
	rhof_mean = zeros(nz+1,1); 
	for l = 1:nt
	    rhof = interp1(zc(:,l),rho(:,l),zf(:,l),'linear','extrap');
	end
	rhof_mean = mean(rhof,2);  % mean density at faces 
    
        % filter u/v and interpolate to constant zc_mean 
        [uiso,viso] = filter_uviso(uiso,viso); % filtering (9-15 hours) 
	for i = 1:nt
	    uiso(:,i) = interp1(zc(:,i),uiso(:,i),zc_mean,'linear','extrap'); 
	    viso(:,i) = interp1(zc(:,i),viso(:,i),zc_mean,'linear','extrap'); 
	end

	% compute perturbation pressure (interpolated) and filter
        pert = compute_pertpress(rho,zc,zf); % perturbation pressure
        pertp = filter_pertp(pert);                       % filter pressure
    
%****** Project velocities
%%%      Case 1: Spectral (Early) 
%        [r2u1(ii,j,:),~] = compute_modfit(uiso,dz,squeeze(Ueig1(ii,j,1:nz,:)),nmodes);
%        [r2v1(ii,j,:),~] = compute_modfit(viso,dz,squeeze(Ueig1(ii,j,1:nz,:)),nmodes);
%        [r2p1(ii,j,:),~] = compute_modfit(pertp,dz,squeeze(Ueig1(ii,j,1:nz,:)),nmodes);
%
%%%      Case 2: FD (Early) 
%        [r2u2(ii,j,:),~] = compute_modfit(uiso,dz,squeeze(Ueig2(ii,j,1:nz,:)),nmodes);
%        [r2v2(ii,j,:),~] = compute_modfit(viso,dz,squeeze(Ueig2(ii,j,1:nz,:)),nmodes);
%        [r2p2(ii,j,:),~] = compute_modfit(pertp,dz,squeeze(Ueig2(ii,j,1:nz,:)),nmodes);
%
%%%      Case 3: FD (Early) 
%        [r2u3(ii,j,:),~] = compute_modfit(uiso,dz,squeeze(Ueig3(ii,j,1:nz,:)),nmodes);
%        [r2v3(ii,j,:),~] = compute_modfit(viso,dz,squeeze(Ueig3(ii,j,1:nz,:)),nmodes);
%        [r2p3(ii,j,:),~] = compute_modfit(pertp,dz,squeeze(Ueig3(ii,j,1:nz,:)),nmodes);
%
%%%      Case 4: FD (Early) 
%        [r2u4(ii,j,:),~] = compute_modfit(uiso,dz,squeeze(Ueig4(ii,j,1:nz,:)),nmodes);
%        [r2v4(ii,j,:),~] = compute_modfit(viso,dz,squeeze(Ueig4(ii,j,1:nz,:)),nmodes);
%        [r2p4(ii,j,:),~] = compute_modfit(pertp,dz,squeeze(Ueig4(ii,j,1:nz,:)),nmodes);

%%      Case 5: FD (Early) 
        % Interpolate to uniform grid
        nzM = round(H/dzi);
	nzM = max([nzM 10]); % Require minimum number points (10)
	zfM = linspace(-H,0,nzM+1)';
        dz1 = diff(zfM);
        zcM = zfM(1:end-1) + dz1/2;

        ufiltM = zeros(nzM,nt); vfiltM = zeros(nzM,nt); pfiltM = zeros(nzM,nt);
        for i = 1:nt
            ufiltM(:,i) = interp1(zc_mean,uiso(:,i),zcM,'linear','extrap');
            vfiltM(:,i) = interp1(zc_mean,viso(:,i),zcM,'linear','extrap');
            pfiltM(:,i) = interp1(zc_mean,pertp(:,i),zcM,'linear','extrap');
        end
        [r2u5(ii,j,:),~] = compute_modfit(ufiltM,dz1,squeeze(Ueig5(ii,j,1:nzM,:)),nmodes);
        [r2v5(ii,j,:),~] = compute_modfit(vfiltM,dz1,squeeze(Ueig5(ii,j,1:nzM,:)),nmodes);
        [r2p5(ii,j,:),~] = compute_modfit(pfiltM,dz1,squeeze(Ueig5(ii,j,1:nzM,:)),nmodes);

    end
end


%% Save output 
save([datapath '/' model '_' xtilestr '_' ytilestr '_r2.mat'],'r2u','r2v','r2p')
%save([datapath '/eigentile_AMZN1942_r2.mat'],'r2u1','r2u2','r2u3','r2u4','r2u5',...
%         'r2v1','r2v2','r2v3','r2v4','r2v5','r2p1','r2p2','r2p3','r2p4','r2p5')

%%%%%%%%%%%%%%%%%%%%%%%% EoF 
