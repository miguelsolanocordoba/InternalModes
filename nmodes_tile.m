clear; clc; close all;
%%NMODES_TILE computes number of resolvable modes
% 
% 
% Created July 13, 2020

%% Paths
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/
datapath = '/data/msolano/matfiles';
figpath = '/data/msolano/figures/modes_hycom/tile2D';

% Defaults
set(groot,'defaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',6)

%% *** INPUT *** %% 
model = 'GLBc0.04';  % Global HYCOM
exptn = 190;         % Experiment number
xtile = 18;          % blki 
ytile = 40;          % blkj
nmodes = 10;       % number of modes

% Convert to character strings
xtilestr = num2str(xtile);
ytilestr = num2str(ytile);
exptnstr = num2str(exptn);

% File name (used for figures and matfiles)
fname = [model '_' exptnstr '_' xtilestr '_' ytilestr]
load([datapath '/' fname '_eig.mat'])

% Read HYCOM tile 
hycom = read_hycom(model,exptn,xtile,ytile); 
lon = hycom.lon; 
lat = hycom.lat; 
depth = hycom.h; 

% Find one point
lont = 314; 
latt = 4; 
[~,indy] = min(abs(lont-lon(1,:))); 
[~,indx] = min(abs(latt-lat(:,1)));

return

%% Number of resolvable modes
dzi = 25; 
[nx,ny,nz,nmode] = size(Weig);
nt=337;
num_modes = zeros(nx,ny); 
nzp = zeros(nx,ny); 

%for i = 1:nx
%   for j = 1:ny
for i = indx
   for j = indy

      if hycom.h(i,j) > 1e4
          continue
      end

      % Grid
      nz = 41;
      zc1 = zeros(nz,nt); zf1 = zeros(nz+1,nt); nzm = zeros(nt,1);
      for l = 1:nt
          zf1(2:end,l) = squeeze(cumsum(hycom.dz(i,j,:,l),3));
          zc1(:,l) = 0.5.*(zf1(1:nz,l) + zf1(2:nz+1,l));
          thk0 = diff(zf1(:,l))>0.01; % number of layers <1cm
          nzm(l) = sum(thk0);
      end 
      nz = min(nzm);  % number of layers (nz) cut to minimum # of non-zero layers
      nzp(i,j) = nz; 
      zf = -flipud(zf1(1:nz+1,:)); % flip and negative down (convention)
      zc = -flipud(zc1(1:nz,:));

      % Compute mean layer
      zf_mean = mean(zf,2);
      zc_mean = mean(zc,2);
      dz = diff(zf_mean);
      H = nansum(dz); % total depth
      N = numel(dz);  % # of layers

      % Uniform grid
      nzM = round(H/dzi); % number of layers
      nzM = max([nzM 20]);
      zfM = linspace(-H,0,nzM+1)';
      zcM = zfM(1:end-1) + dzi/2;
      %nz = min([numel(zfM) 213]);
      nz = min([numel(zfM) 265]);

%% Compute number of resolvable modes (num_modes & num_modesU)
      % Weig
      Weig1 = zeros(nz,nmodes); 
      Weig1 = squeeze(Weig(i,j,1:nz,1:nmodes)); 
      Ueig1 = zeros(nz,nmodes); 
      Ueig1 = squeeze(Ueig(i,j,1:nz,1:nmodes)); 
      nres = 1; % Start checking from mode 2
      resflag = true;
      nzH = numel(zf_mean);
      A = zeros(nz-1,1);
      while resflag
         nres = nres + 1;
         A = Weig1(1:end-1,nres+1).*Weig1(2:end,nres+1);
         ind = find(A<0);
         wzeros = 0.5*(zfM(ind)+zfM(ind+1));
         %wzeros = 0.5*(zcM(ind)+zcM(ind+1));
	 wzeros = [-H;wzeros;0]; % add end-points (surface/bottom)
       
         for ii = 1:numel(wzeros)-1
            if sum(zf_mean>wzeros(ii) & zf_mean<wzeros(ii+1))==0;
                 resflag = false;
		 num_modes(i,j) = nres;
            end
         end


	 if nres==9
            num_modes(i,j) = nmodes; 
	    break
	 end
       
      end %while

      a = 2.35;
      figure
      plot(Ueig1(:,nres),zfM); hold on
      plot(Ueig1(:,nres+1),zfM); hold on
      plot(zeros(numel(wzeros),1),wzeros,'*b')
      plot(linspace(-a,a,nzH),ones(1,nzH).*zf_mean(1),'k','LineWidth',1)
      plot(linspace(-a,a,nzH),ones(1,nzH).*zf_mean(2),'k','LineWidth',1)
      plot(linspace(-a,a,nzH),ones(1,nzH).*zf_mean(3),'k','LineWidth',1)
      plot(linspace(-a,a,nzH),ones(1,nzH).*zf_mean(4),'k','LineWidth',1)
      plot(linspace(-a,a,nzH),ones(1,nzH).*zf_mean(5),'k','LineWidth',1)
      plot(zeros(nzH,1),zf_mean,'+k')
      plot(zeros(nzH,1),zf_mean,'k','LineWidth',0.5)
      title(['Resolvable modes= ' num2str(num_modes(i,j))])
      xlabel('U_{eig}'); ylabel('Depth [m]'); ylim([-H 0])
      pbaspect([1 2 1]);
      legend(['U_{eig} Mode ' num2str(nres)],['U_{eig} Mode ' num2str(nres+1)],...
       'inflection points','HYCOM layers','Location','EastOutside')
      print('-dpng','-r300','Ueig.png')

      %% Ueig
%      Ueig1 = zeros(nz,nmodes);
%      Ueig1 = squeeze(Ueig(i,j,1:nz,1:nmodes));
%      nres = 1; % Start checking from mode 2
%      resflag = true;
%      nzH = numel(zf_mean);
%      A = zeros(nz-1,1);
%      while resflag
%         nres = nres + 1;
%         A = Ueig1(1:end-1,nres+1).*Ueig1(2:end,nres+1);
%         ind = find(A<0);
%         %wzeros = 0.5*(zfM(ind)+zfM(ind+1));
%         wzeros = 0.5*(zcM(ind)+zcM(ind+1));
%         wzeros = [-H;wzeros;0]; % add end-points (surface/bottom)
%         %wzeros = wzeros(1:end-1) - dzi;
%         wzeros = wzeros(1:end-1);
%
%         for ii = 1:numel(wzeros)-1
%            if sum(zf_mean>wzeros(ii) & zf_mean<wzeros(ii+1))==0;
%                 resflag = false;
%                 num_modesU(i,j) = nres;
%            end
%         end
%
%         if nres==9
%            num_modesU(i,j) = nmodes;
%            break
%         end
%
%      end %while

%       figure
%      plot(Weig1(:,1),zfM); hold on
%      plot(Weig1(:,2),zfM); hold on
%      plot(Ueig1(:,nres+1),zfM); hold on
%      plot(Ueig1(:,9),zfM); hold on
%      plot(zeros(nzH,1),zf_mean,'*k')
%      plot(zeros(numel(wzeros),1),wzeros,'*b')
%      plot(zeros(nzH,1),zf_mean,'k','LineWidth',2)
%      title(['num modes= ' num2str(num_modesU(i,j))])
%      xlabel('U_{eig}'); ylabel('Depth [m]'); ylim([-H 0])
%      pbaspect([1 2 1]);
%      legend('Mode 1','Mode 2','Mode 3','Grid','zeros','Location','SouthWest')
%      print('-dpng','-r300','num_modesU.png')


   end % ny
end  % nx

%save([datapath '/num_modes_' fname '.mat'],'num_modes','nzp')

return
%% Figure

contvec = [-3000 -2000 -1000 -500];

% Bathymetry
figure; colormap('jet')
pcolor(lon,lat,-hycom.h); shading flat; hold on
colorbar; caxis([-5000 0])
[C,h]=contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
clabel(C,h);
plot(lon(78,110),lat(78,110),'*k')
title('Bathymetry')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300',['bathy_' fname '.png'])

% Number of layers
figure; colormap(jet(15))
pcolor(lon,lat,nzp); shading flat; hold on
cb = colorbar;
set(cb,'xtick',[26.5:3:40.5],'xticklabel',{'26','29','31','34','37','40'})
caxis([26 41]);
[C,h]=contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
clabel(C,h);
title('Layers >1mm')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300',['layer_' fname '.png'])

%figure; colormap(jet(12))
%pcolor(lon,lat,nzp); shading flat; hold on
%colorbar; caxis([27 38]);
%contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
%title('Layers >1mm')
%xlabel('Longitude'); ylabel('Latitude')

% Number of modes resolved
figure; colormap(jet(nmodes+1));
pcolor(lon,lat,num_modes); shading flat; hold on
cb = colorbar;
set(cb,'xtick',[0.5:1:nmodes+0.5],'xticklabel',{'0','1','2','3','4','5','6','7','8','9','10'})
caxis([0 nmodes+1])
[C,h]=contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
clabel(C,h);
plot(lon(78,110),lat(78,110),'*k','MarkerSize',10)
plot(lon(78,110),lat(78,110),'ok','MarkerSize',10)
title('Number of resolvable modes')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300',['num_modes_' fname '.png'])

%figure; colormap(jet(5));
%pcolor(lon,lat,num_modes); shading flat; hold on
%colorbar; caxis([1 5])
%contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
%title('Number of resolvable modes')
%xlabel('Longitude'); ylabel('Latitude')

return
% R2 maps
figure; colormap('jet')
%pcolor(lon,lat,r2_uniform); shading flat; hold on
%pcolor(lon,lat,squeeze(sqrt(r2u5(:,:,5).^2 + r2v5(:,:,5).^2))); shading flat; hold on
pcolor(lon,lat,squeeze(r2u5(:,:,5))); shading flat; hold on
colorbar; caxis([0.75 1])
contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
title('u regression coefficient (R^2)')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300','R2u_amzn1942.png')

figure; colormap('jet')
%pcolor(lon,lat,r2_uniform); shading flat; hold on
%pcolor(lon,lat,squeeze(sqrt(r2u5(:,:,5).^2 + r2v5(:,:,5).^2))); shading flat; hold on
pcolor(lon,lat,squeeze(r2v5(:,:,5))); shading flat; hold on
colorbar; caxis([0.75 1])
contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
title('v regression coefficient (R^2)')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300','R2v_amzn1942.png')

figure; colormap('jet')
%pcolor(lon,lat,r2_uniform); shading flat; hold on
pcolor(lon,lat,squeeze(sqrt(r2u5(:,:,5).^2 + r2v5(:,:,5).^2))); shading flat; hold on
colorbar; caxis([0.75 1])
contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
title('u/v regression coefficient (R^2)')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300','R2uv_amzn1942.png')

figure; colormap('jet')
pcolor(lon,lat,squeeze(r2p5(:,:,5))); shading flat; hold on
colorbar; caxis([0.8 1])
contour(lon,lat,-hycom.h,[contvec],'k','LineWidth',1.5);
title('p regression coefficient (R^2)')
xlabel('Longitude'); ylabel('Latitude')
print('-dpng','-r300','R2p_amzn1942.png')

%% Move all figures to dirout
system(['rm ' figpath '/*.png']);
system(['mv *.png ' figpath]);

%% EoF 
