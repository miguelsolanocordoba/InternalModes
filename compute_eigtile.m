clc; clear; close all; 
%%COMPUTE_EIGTILE solves the EVP for a HYCOM tile
% 
% 
% Created July 13, 2020

%% Defaults
warning('off'); 

addpath /data/msolano/Matlab/ocean_physics
addpath /data/msolano/toolbox/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
figpath = '/data/msolano/figures/test';

%% Read HYCOM variables
hycom = read_hycom(); 
[nx,ny,nz,nt] = size(hycom.rho); 

% Constants
g = 9.806;                      % gravity 
omega = 12.1408331767/24/3600;  % M2 frequency
%nfiles = nx*ny; 
nfiles = 20*20;
nmodes = 5; 

% Variables to plot
CO = zeros(nx,ny,nmodes); 
CgO = zeros(nx,ny,nmodes); 
LO = zeros(nx,ny,nmodes); 
CE = zeros(nx,ny,nmodes); 
CgE = zeros(nx,ny,nmodes); 
LE = zeros(nx,ny,nmodes); 

fprintf('\nComputing eigenvalues...\n')
count = 0; 
for i = 1:10
    for j = 1:10
	count = count + 1; 
	fprintf('%d/%d\n',count,nfiles)
	   
        f = coriolis(hycom.lat(i,j)); % Coriolis frequency
    
        %% Pre-processing
        [nx,ny,nz,nt] = size(hycom.rho); 

        % compute grid 
        zc1 = zeros(nz,nt); zf1 = zeros(nz+1,nt); 
        for l = 1:nt
            zf1(2:end,l) = squeeze(cumsum(hycom.dz(i,j,:,l),3)); 
            zc1(:,l) = 0.5.*(zf1(1:nz,l) + zf1(2:nz+1,l));
        end
        thk0=diff(zf1(:,1))>0.1; % mask layers < 10cm 
        nz = sum(thk0);          % number of layers > 0.1m
        
        zf = -flipud(zf1(1:nz+1,:)); clear zf1;
        zc = -flipud(zc1(1:nz,:)); clear zc1
    
        zf_mean = mean(zf,2); zc_mean = mean(zc,2); 
        dz = diff(zf_mean); 
        H = nansum(dz); 
        N = numel(dz); 
        
        % flip variables upside-down 
        uiso = flipud(squeeze(hycom.uiso(i,j,1:nz,:))); 
        viso = flipud(squeeze(hycom.viso(i,j,1:nz,:))); 
        rho = flipud(squeeze(hycom.rho(i,j,1:nz,:))); 
        rho_mean = flipud(squeeze(mean(hycom.rho(i,j,1:nz,:),4)));
	
	rhof_mean = zeros(nz+1,1); 
	for l = 1:nt
	    rhof = interp1(zc(:,l),rho(:,l),zf(:,l),'linear','extrap');
	end
	rhof_mean = mean(rhof,2);  % mean density at faces 
    
        % compute perturbation pressure and apply band-pass filter 
        [uiso,viso] = filter_uviso(uiso,viso); % filtering (9-15 hours) 
        pert = comp_pertpress(rho,dz);         % perturbation pressure
        pertp = filter_pertp(pert);            % filter pressure
    
        %% Compute eigenfunctions 
        % Oladeji
        [C,Cg,L,~,~] = compute_eigen(rho_mean,zf_mean,f,omega);
        CO(i,j,:) = C(1:5); 
        CgO(i,j,:) = Cg(1:5); 
	LO(i,j,:) = L(1:5); 
        	
        % Early 
        im = InternalModes(rhof_mean,zf_mean,zc_mean,hycom.lat(i,j)); %,...
	im.normalization = 'wMax';
	im.method = 'finiteDifference'; 
%             'method','finiteDifference','orderOfAccuracy','4');
        [~,~,h,K] = im.ModesAtFrequency(omega); 
        CgE(i,j,:) = sqrt(g*h(1:5)); 
        LE(i,j,:) = 2*pi./K(1:5);

    end
end


%% Plots
lon = hycom.lon(1:10,1:10); 
lat = hycom.lat(1:10,1:10); 

figure
subplot(321)
pcolor(lon,lat,squeeze(LO(1:10,1:10,1))./1000); shading interp
title('Mode 1: wave-length (Oladeji)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(322)
pcolor(lon,lat,squeeze(LE(1:10,1:10,1))./1000); shading interp
title('Mode 1: wave-length (Early)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(323)
pcolor(lon,lat,squeeze(LO(1:10,1:10,1))./1000); shading interp
title('Mode 2: wave-length (Oladeji)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(324)
pcolor(lon,lat,squeeze(LE(1:10,1:10,1))./1000); shading interp
title('Mode 2: wave-length (Early)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(325)
pcolor(lon,lat,squeeze(LO(1:10,1:10,1))./1000); shading interp
title('Mode 3: wave-length (Oladeji)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(326)
pcolor(lon,lat,squeeze(LE(1:10,1:10,1))./1000); shading interp
title('Mode 3: wave-length (Early)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

print('L.png','-r300','-dpng')

figure
subplot(321)
pcolor(lon,lat,CgO(1:10,1:10,1)); shading interp
title('Mode 1: group velocity (Oladeji)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(322)
pcolor(lon,lat,CgE(1:10,1:10,1)); shading interp
title('Mode 1: group velocity (Early)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(323)
pcolor(lon,lat,CgO(1:10,1:10,2)); shading interp
title('Mode 2: group velocity (Oladeji)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(324)
pcolor(lon,lat,CgE(1:10,1:10,2)); shading interp
title('Mode 2: group velocity (Early)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(325)
pcolor(lon,lat,CgO(1:10,1:10,3)); shading interp
title('Mode 3: group velocity (Oladeji)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

subplot(326)
pcolor(lon,lat,CgE(1:10,1:10,3)); shading interp
title('Mode 3: group velocity (Early)') 
xlabel('longitude'); ylabel('latitude') 
cb = colorbar; %caxis([0 1000]) 

print('Cg.png','-r300','-dpng')

system(['mv *.png ' figpath '/'])
