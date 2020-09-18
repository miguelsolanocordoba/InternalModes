clc; clear; close all; 
%%UVFIT_COMP compararison to fit modal velocities based on the EVP.
%  
% Created: Miguel Solano, May 31, 2020.


%% Defaults
warning('off');
set(groot,'defaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',6)


%% Paths
addpath /data/msolano/forOladeji % Matfiles (data)
addpath /data/msolano/Matlab    v% All functions!
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/

dirout = '/data/msolano/figures/modes_hycom/';

% plotting options (1=yes, 0=no)
plotini = 0; % Stratification and filtering 
ploteig = 1; % Eigenvalues (Ueig,Weig) 
plotfit = 1; % Velocity and fit (pcolor and time series)
plotsta = 1; % Statistics (R2 and S2)
fntsz = 6;   % legend font size 


%% Load variables
load profile_loc3.mat     

uiso = profile.uiso; 
viso = profile.viso; 
ufilt = profile.ufilt; 
vfilt = profile.vfilt; 
ufiltint = profile.ufiltint; 
vfiltint = profile.vfiltint; 

% Dimensions and constants
[nz,nt] = size(uiso); 

% Constants
g = 9.806;                        % gravity 
omega = 12.1408331767/24/3600;   % M2 frequency 

%% Cases
%cname = {'Spectral','FD (Early)','FD (Oladeji)','Uniform'};
cname = {'FD (Interp)','Uniform (Interp)','FD (Eta)','Uniform (Eta)'};


%% Initialize
% Time Vector (September 1-15, 2016)
t = datenum(2016,9,1):datenum(0,0,0,1,0,0):datenum(2016,9,15); 
nt = numel(t); 

% Vertical grid (surface to bottom) 
dz = diff(profile.zf_mean);         % layer thickness 
H = nansum(dz);        % total depth 
N = numel(dz);         % number of layers 
lat = profile.latitude;              % Latitude 

% Interpolate density to mean layer
rhoint = zeros(nz,nt);
rho = zeros(nz,nt);
for i = 1:nt
    rho(:,i) = profile.rho(:,i) + 1000; 
    rhoint(:,i) = interp1(profile.zc(:,i),profile.rho(:,i)+1000,...
                          profile.zc_mean,'linear','extrap');
end
rho_mean = mean(rhoint,2);
rhof_mean = interp1(profile.zc_mean,rho_mean,profile.zf_mean,'linear','extrap');

% Compute perturbation pressure
pert = compute_pertpress(rhoint,profile.zc,profile.zf); 
pertp = filter_pertp(pert); 

pert2 = compute_pertpress2(rho,profile.zc,profile.zf); 
pertp2 = filter_pertp(pert2); 

pert(:,1)-pert2(:,1)

%% PLOT: perturbation pressure
[x,y] = meshgrid(t,profile.zc_mean); 
fntsz = 8; 
ymin = -H; 

figure
plot(t,pert(end,:),'k'); hold on
plot(t,pert2(end,:),'b');
text(t(10),100,'SURFACE')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('P" [N/m^2]')
legend('Interp','Eta','Location','SouthEast')
print('p_surface.png','-r300','-dpng')

figure
plot(t,pert(1,:),'k'); hold on
plot(t,pert2(1,:),'b');
text(t(10),-300,'BOTTOM')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('P" [N/m^2]')
legend('Interp','Eta','Location','NorthEast')
print('p_bottom.png','-r300','-dpng')

figure; colormap('jet') 
subplot(211)
pcolor(x,y,pert); shading flat; colorbar; 
title('Perturbation Pressure (Interp)') 
datetick('x','dd'); xlim([t(1) t(end)]);
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(212)
pcolor(x,y,pert2); shading flat; colorbar; 
title('Perturbation Pressure (Eta)') 
datetick('x','dd'); xlim([t(1) t(end)]);
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);
print('p_contour.png','-r300','-dpng')

%system(['rm ' dirout '*.png']);
system(['mv *.png ' dirout]);

%% Compute eigenfunctions

% Initilized Early's InternalModes
zf = profile.zf_mean; 
zc = profile.zc_mean; 

imspectral = InternalModes(rhof_mean,zf,zf,lat,'method','spectral');
imf = InternalModes(rhof_mean,zf,zf,lat,'method','finiteDifference');

f = imf.f0; 
rho0 = imf.rho0; 

%% Case 1: SPECTRAL
% Eigenfunctions
[~,Weig1,h1,k1] = imspectral.ModesAtFrequency(omega); % omega-const EVP 
Ueig1 = compute_ueig(Weig1,dz);                  % Normalized Ueig

% Eigenvalues
L1 = 2*pi./k1;                                   % wave-length 
C1 = omega./k1;                                  % phase-speed
Cg1 = (omega^2-f^2)./(omega*k1);                 % group-speed 

%% Case 2: FD (Early) 
% Eigenfunctions
[~,Weig2,h2,k2] = imf.ModesAtFrequency(omega);   % omega-const EVP
Ueig2 = compute_ueig(Weig2,dz);                  % Normalized Ueig

% Eigenvalues
L2 = 2*pi./k2;                                   % wave-length 
C2 = omega./k2;                                  % phase-speed
Cg2 = (omega^2-f^2)./(omega*k2);                 % group-speed 
%Ce2 = sqrt(omega^2-f^2)./(k2.^2);                % eigen-speed
Ce2 = sqrt(g*h2);


%% Case 3: FD Weig (Oladeji) 
% Eigenfunctions
S = compute_eigen(rho_mean,rho0,zf,f,omega);
Ueig3 = S.Ueig;

% Eigenvalues
L3 = S.L;
C3 = S.C;
Cg3 = S.Cg;
k3 = S.k;
Ce3 = sqrt(omega^2-f^2)./(k3.^2);                % eigen-speed

%% Case 4: Same as (5) but using ufiltint
% Uniform grid
dzi = 25;                       % 25-m spacing (dz)
nzM = round(H/dzi);
zfM = linspace(-H,0,nzM+1)';
dz1 = diff(zfM);
zcM = zfM(1:end-1) + dz1/2;

% Interpolate stratification to uniform grid
N2M = interp1(zf,S.N2,zfM,'linear','extrap');

% Compute eigenfunctions
[C4,Cg4,L4,Weig4,Ueig4t] = sturm_liouville_hyd_normalize(omega,sqrt(N2M),dz1(1),f);
Ueig4t(:,2:2:end) = -Ueig4t(:,2:2:end);
Ueig4 = Ueig4t(:,2:end);

% Interpolate to uniform grid
pfiltM = zeros(nzM,nt);
pfiltM2 = zeros(nzM,nt);
for i = 1:nt
    pfiltM(:,i) = interp1(profile.zc(:,i),pertp(:,i),zcM,'linear','extrap');
    pfiltM2(:,i) = interp1(profile.zc(:,i),pertp2(:,i),zcM,'linear','extrap');
end

%% Mode fitting and statistics
% Compute the fit
nmodes = 5; 

% CASE 1
[r2u1,p1] = compute_modfit(pertp,dz,Ueig3,nmodes); 

% CASE 2
[r2u2,p2] = compute_modfit(pfiltM,dz1,Ueig4,nmodes); 

% CASE 3
[r2u3,p3] = compute_modfit(pertp2,dz,Ueig3,nmodes); 

% CASE 4
[r2u4,p4] = compute_modfit(pfiltM2,dz1,Ueig4,nmodes); 


% Fit of all modes (1 to nmodes)
pfit1 = sum(p1,3); 
pfit2 = sum(p2,3); 
pfit3 = sum(p3,3); 
pfit4 = sum(p4,3); 



%% PLOTS
ymin = -H; 

%% plotini: stratification, filtered velocities
if plotini

% Stratification
figure
subplot(121)
plot(im.rho,im.z,'r'); 
title('Density'); xlabel('\sigma_2 [kg/m^3'); ylabel('Depth [m]')
ylim([-H 0]);
subplot(122)
semilogx(sqrt(im.N2),im.z,'r'); 
title('Stratification'); xlabel('N [1/sec]'); ylabel('Depth [m]')
ylim([-H 0]);
print('stratification.png','-r300','-dpng')

% Pressure (time-series)
figure
subplot(211)
plot(t,pert(end,:),'k'); hold on
plot(t,pertp(end,:),'b'); 
text(t(280),150,'SURFACE')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('P" [N/m^2]') 
subplot(212)
plot(t,pert(1,:),'k'); hold on
plot(t,pertp(1,:),'b'); 
text(t(10),5,'BOTTOM')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('P" [N/m^2]') 
print('pfiltering.png','-r300','-dpng') 

end % if plotini


%% ploteig: eigenfunctions (Ueig, Weig)
if ploteig
figure;
subplot(151)
plot(Ueig1(:,1),zc,'k'); hold on
plot(Ueig2(:,1),zc,'b');
plot(Ueig3(:,1),zc,'--r');
plot(Ueig4(:,1),zcM,'--g');
title('Mode 1'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
legend(cname{1},cname{2},cname{3},cname{4},'Box','on',...
       'Position',[0.04 0.67 0.06 0.1],'FontSize',5)

subplot(152)
plot(Ueig1(:,2),zc,'k'); hold on
plot(Ueig2(:,2),zc,'b');
plot(Ueig3(:,2),zc,'--r');
plot(Ueig4(:,2),zcM,'--g');
title('Mode 2'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(153)
plot(Ueig1(:,3),zc,'k'); hold on
plot(Ueig2(:,3),zc,'b');
plot(Ueig3(:,3),zc,'--r');
plot(Ueig4(:,3),zcM,'--g');
title('Mode 3'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(154)
plot(Ueig1(:,4),zc,'k'); hold on
plot(Ueig2(:,4),zc,'b');
plot(Ueig3(:,4),zc,'--r');
plot(Ueig4(:,4),zcM,'--g');
title('Mode 4'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(155)
plot(Ueig1(:,5),zc,'k'); hold on
plot(Ueig2(:,5),zc,'b');
plot(Ueig3(:,5),zc,'--r');
plot(Ueig4(:,5),zcM,'--g');
title('Mode 5'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

print('Ueig.png','-r300','-dpng')

end % if ploteig

%% plotfit: fitted velocities vs original velocities
if plotfit 

[x,y] = meshgrid(t,profile.zc_mean); 
[x1,y1] = meshgrid(t,zcM); 
cmin = -10; cmax = 10; 
ymin = -H; 

% U (color contours)
figure
subplot(511)
pcolor(x,y,pertp); shading interp; colorbar; caxis([cmin cmax])
title('Band-passed P_{pert}');
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(512)
pcolor(x,y,pfit1); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{1}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(513)
pcolor(x1,y1,pfit2); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{2}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,pfit3); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{3}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x1,y1,pfit4); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{4}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

print('pfit_vs_pert.png','-r300','-dpng')

% V (color contours)

% Time-series
% surface
figure 
plot(t,pertp(1,:),'k','LineWidth',2); hold on
plot(t,pfit1(1,:),'b'); 
plot(t,pfit2(1,:),'r') 
plot(t,pfit4(1,:),'g')
datetick('x','dd'); ylabel('p [N/m^2]'); pbaspect([5 1 1])
legend('p_{pert}',cname{1},cname{2},cname{4},'Location','EastOutside')
print('p_bot.png','-r300','-dpng') 

% middepth
figure
plot(t,pertp(8,:),'k','LineWidth',2); hold on
plot(t,pfit1(8,:),'b'); 
plot(t,pfit2(8,:),'r') 
plot(t,pfit4(8,:),'g')
datetick('x','dd'); ylabel('p [N/m^2]'); pbaspect([5 1 1])
legend('p_{pert}',cname{1},cname{2},cname{4},'Location','EastOutside')
print('p_mid.png','-r300','-dpng') 

% bottom
figure
plot(t,pertp(end,:),'k','LineWidth',2); hold on
plot(t,pfit1(end,:),'b'); 
plot(t,pfit2(end,:),'r') 
plot(t,pfit4(end,:),'g')
datetick('x','dd'); ylabel('p [N/m^2]');  pbaspect([5 1 1])
legend('p_{pert}',cname{1},cname{2},cname{4},'Location','EastOutside')
print('p_sur.png','-r300','-dpng') 

end % if plotfit


%% plotsta: coefficient of determination R^2 and S^2
if plotsta

figure
bar([r2u1 r2u2 r2u3 r2u4])
title('R^2 for P_{fit}'); xlabel('Modes'); ylabel('R^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',8); xticklabels({'1','1-2','1-3','1-4','1-5'})
legend(cname{1},cname{2},cname{3},cname{4},'Location','NorthWest')

print('r2p.png','-r300','-dpng') 

end  % if plotsta


%% Move all figures to dirout
%system(['rm ' dirout '*.png']);
system(['mv *.png ' dirout]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
