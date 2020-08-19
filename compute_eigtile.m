clc; clear; close all; 
%%COMPUTE_EIGTILE solves the EVP for a HYCOM tile
% 
% 
% Created July 13, 2020

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
hycom = read_hycom(); 
[nx,ny,nz,nt] = size(hycom.rho); 
nx = 80; % over-ride x-size 
ny = 80; % over-ride y-size

% Constants
g = 9.806;                      % gravity 
omega = 12.1408331767/24/3600;  % M2 frequency
nfiles = nx*ny/4; 
nmodes = 5; 

% Variables to plot
r2M = zeros(nx/2,ny/2,nmodes); 
r2W = zeros(nx/2,ny/2,nmodes); 
r2J = zeros(nx/2,ny/2,nmodes); 
r2O = zeros(nx/2,ny/2,nmodes); 
%CJ = zeros(nx,ny,nmodes); 
%CgJ = zeros(nx,ny,nmodes); 
%LJ = zeros(nx,ny,nmodes); 
%CO = zeros(nx,ny,nmodes); 
%CgO = zeros(nx,ny,nmodes); 
%LO = zeros(nx,ny,nmodes); 

fprintf('\nComputing eigenvalues...\n')
count = 0; 
counti = 0;
for ii = 1:2:nx
    counti = counti + 1; 
    countj = 0; 
    for j = 1:2:ny
	countj = countj + 1; 
	count = count + 1; 
	fprintf('%d/%d\n',count,nfiles)
	   
        f = coriolis(hycom.lat(ii,j)); % Coriolis frequency
    
        %% Pre-processing
        [~,~,nz,nt] = size(hycom.rho); 

        % compute grid 
        zc1 = zeros(nz,nt); zf1 = zeros(nz+1,nt); 
        for l = 1:nt
            zf1(2:end,l) = squeeze(cumsum(hycom.dz(ii,j,:,l),3)); 
            zc1(:,l) = 0.5.*(zf1(1:nz,l) + zf1(2:nz+1,l));
        end
        thk0=diff(zf1(:,1))>0.1; % mask layers < 10cm 
        nz = sum(thk0)-1;          % number of layers > 0.1m
        
        zf = -flipud(zf1(1:nz+1,:)); clear zf1;
        zc = -flipud(zc1(1:nz,:)); clear zc1
    
        zf_mean = mean(zf,2); zc_mean = mean(zc,2); 
        dz = diff(zf_mean); 
        H = nansum(dz); 
        N = numel(dz); 
        
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
    
        % compute perturbation pressure and apply band-pass filter 
        [uiso,viso] = filter_uviso(uiso,viso); % filtering (9-15 hours) 
%        pert = comp_pertpress(rho,dz);         % perturbation pressure
%        pertp = filter_pertp(pert);            % filter pressure
    
%%      Compute eigenfunctions 

%%      Case 1: FD Oladeji
        S = compute_eigen(rho_mean,zf_mean,f,omega);
        [r2O(counti,countj,:),~] = compute_modfit(uiso,dz,S.Ueig,5); 	
        	
%%      Case 2: FD (Early) 
        imFD = InternalModes(rhof_mean,zf_mean,zf_mean,hycom.lat(ii,j),...
	                     'method','finiteDifference',...
			     'orderOfAccuracy',2); 
%	imFD.method = 'finiteDifference'; 
        [~,WeigFD,~,~] = imFD.ModesAtFrequency(omega); 
	UeigFD = compute_ueig(WeigFD,dz); 
	  
        [r2J(counti,countj,:),~] = compute_modfit(uiso,dz,UeigFD,nmodes); 	
	 
%%      Case 3: WKB (Early)
        im = InternalModes(rhof_mean,zf_mean,zf_mean,hycom.lat(ii,j)); %,...
        [~,Weig,h,k] = im.ModesAtFrequency(omega);
        Ueig = compute_ueig(Weig,dz);

        [r2W(counti,countj,:),~] = compute_modfit(uiso,dz,Ueig,nmodes);

%%      Case 4: Uniform (Maarten)
        dzi = 25; % uniform res = 25 meters
	nzM = round(H/dzi); % number of layers
        zfM = linspace(-H,0,nzM+1)';
        dz1 = diff(zfM);
        zcM = zfM(1:end-1) + dz1/2;
        N2M = interp1(zf_mean,im.N2,zfM,'linear','extrap');
        [C6,Cg6,L6,Weig6,Ueig6t] = sturm_liouville_hyd_normalize(omega,sqrt(N2M),dz1(1),f);
        Ueig6t(:,2:2:end) = -Ueig6t(:,2:2:end);
        UeigM = Ueig6t(:,2:end);
        
        % Interpolate to uniform grid
        ufiltM = zeros(nzM,nt);
        for i = 1:nt
            ufiltM(:,i) = interp1(zc_mean,uiso(:,i),zcM,'linear','extrap');
        end

	[r2M(counti,countj,:),~] = compute_modfit(ufiltM,dz1,UeigM,nmodes);

    end
end

lon = hycom.lon(1:2:nx,1:2:ny); 
lat = hycom.lat(1:2:nx,1:2:ny); 
depth = hycom.h(1:2:nx,1:2:ny); 

%% Save output 
save('eigtile_whole2.mat','lon','lat','depth','r2W','r2J','r2O','r2M')
system(['mv eigtile_whole2.mat ' figpath '/']);

% Compute R^2 averages
r2J_mean = squeeze(real(mean(r2J,[1 2])));
r2O_mean = squeeze(real(mean(r2O,[1 2])));
r2M_mean = squeeze(real(mean(r2M,[1 2])));
r2W_mean = squeeze(real(mean(r2W,[1 2])));


%% Plots
figure
bar([r2W_mean r2J_mean r2O_mean r2M_mean])
title('Regression Coefficient (R^2)'); ylim([0 1])
set(gca,'FontSize',10); xticklabels({'1','1-2','1-3','1-4','1-5'});
legend('WKB','FD (Early)','FD (Oladeji)','Uniform','Location','NorthWest'); %,...
print('-dpng','-r300','r2_mean.png')

figure
subplot(221)
pcolor(lon,lat,squeeze(r2O(:,:,5))); shading interp
title('FD Oladeji')    
xlabel('longitude'); ylabel('latitude')
cb = colorbar; caxis([0.6 1]) 
subplot(222)
pcolor(lon,lat,real(squeeze(r2J(:,:,5)))); shading interp
title('FD Early')    
xlabel('longitude'); ylabel('latitude')
cb = colorbar; caxis([0.6 1]); 
subplot(223)
pcolor(lon,lat,real(squeeze(r2W(:,:,5)))); shading interp
title('WKB Early')
xlabel('longitude'); ylabel('latitude')
cb = colorbar; caxis([0.6 1]);
subplot(224)
pcolor(lon,lat,squeeze(r2M(:,:,5))); shading interp
title('Uniform')
xlabel('longitude'); ylabel('latitude')
cb = colorbar; caxis([0.6 1]); 

print('R2_maps.png','-r300','-dpng')

%figure
%subplot(131)
%pcolor(lon,lat,squeeze(LO(1:nx,1:ny,1))./1000); shading interp
%title('Mode 1 wavelength (L)') 
%xlabel('longitude'); ylabel('latitude') 
%cb = colorbar; %caxis([0 1000]) 
%subplot(132)
%pcolor(lon,lat,squeeze(LO(1:nx,1:ny,2))./1000); shading interp
%title('Mode 3 wavelength (L)') 
%xlabel('longitude'); ylabel('latitude') 
%cb = colorbar; %caxis([0 1000]) 
%subplot(133)
%pcolor(lon,lat,squeeze(LO(1:nx,1:ny,3))./1000); shading interp
%title('Mode 3 wavelength (L)') 
%xlabel('longitude'); ylabel('latitude') 
%cb = colorbar; %caxis([0 1000]) 
%subplot(154)
%pcolor(lon,lat,squeeze(LO(1:nx,1:ny,4))./1000); shading interp
%title('Mode 4 wavelength (L)') 
%xlabel('longitude'); ylabel('latitude') 
%cb = colorbar; %caxis([0 1000]) 
%subplot(155)
%pcolor(lon,lat,squeeze(LO(1:nx,1:ny,5))./1000); shading interp
%title('Mode 5 wavelength (L)') 
%xlabel('longitude'); ylabel('latitude') 
%cb = colorbar; %caxis([0 1000]) 
%
%print('L.png','-r300','-dpng')


system(['mv *.png ' figpath '/'])
