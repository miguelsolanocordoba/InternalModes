clc; clear; close all; 
%%UVFIT_COMP computes harmonic fits (u,v,p) from HYCOM profiles.
% UVFIT_COMP is used to compare the eigenfunctions computed using 
% Jeffrey Early's InternalModes solver with other finite difference
% (FD) schemes. 
%
% INPUT OPTIONS: (=1 plot, =0 no plot) 
% plotini => plots density, stratification freq. and filtering.
% ploteig => plots horizontal eigenvalues (Ueig)
% plotfit => plots harmonic fit (pcolor and time-series)
% plotsta => plots statistics
% 
% Created: Miguel Solano, May 31, 2020.

%% Defaults
warning('off');
set(groot,'defaultLineLineWidth',1.5);
set(0,'defaultAxesFontSize',6)

%% Paths
addpath /data/msolano/forOladeji % Matfiles (data)
addpath /data/msolano/Matlab    v% All functions!
addpath /data/msolano/toolbox/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions

%% Input Options
% Location (=1 North Atlantic, =2 South Pacific)
loc = 2;  locstr = num2str(loc); 

% plotting options (1=yes, 0=no)
plotini = 1; % Stratification and filtering 
ploteig = 1; % Eigenvalues (Ueig,Weig) 
plotfit = 1; % Velocity and fit (pcolor and time series)
plotsta = 1; % Statistics (R2 and S2)
fntsz = 6;   % legend font size 

dirout = '/data/msolano/figures/modes_hycom/';

%% Load variables
load(['profile_loc' locstr '.mat'])
%load profile_loc1.mat     

uiso = profile.uiso; 
viso = profile.viso; 
ufilt = profile.ufilt; 
vfilt = profile.vfilt; 
ufiltint = profile.ufiltint; 
vfiltint = profile.vfiltint; 

% Dimensions and constants
[nz,nt] = size(uiso); 

% Constants
g = 9.81;                        % gravity 
omega = 12.1408331767/24/3600;   % M2 frequency 

%% Cases
cname = {'Early W','Early U','Early Int','FD2 W','FD2 U','FD2 Int'};


%% Initialize
% Time Vector (September 1-15, 2016)
t = datenum(2016,9,1):datenum(0,0,0,1,0,0):datenum(2016,9,15); 
nt = numel(t); 

% Vertical grid (surface to bottom) 
zf = profile.zf_mean;  % cell faces 
zc = profile.zc_mean;  % cell centers
dz = diff(zf);         % layer thickness 
H = nansum(dz);        % total depth 
N = numel(dz);         % number of layers 

% Density and latitude
rho_mean = profile.rho_mean;  % Mean density 
lat = profile.latitude;  % Latitude 
rhof_mean = interp1(zc,rho_mean,zf,'linear','extrap');

% Input/Output grids
zIn = zc; 
zOut = zf; 

% Initilized Early's InternalModes
imWKB = InternalModes(rhof_mean,zOut,zOut,lat)
im = InternalModes(rhof_mean,zOut,zOut,lat,'method','finiteDifference',...
                   'orderOfAccuracy',4)
imU = InternalModes(rho_mean,zIn,zIn,lat,'method','finiteDifference',...)
                   'orderOfAccuracy',4)

% Compute perturbation pressure
%pertp = comp_pertpress(profile.rho,dz); 

%% Compute eigenmodes
% Case 1: compute Ueig from Weig (Early)
[~,Weig1,~,~] = im.ModesAtFrequency(omega); % w-const EVP 
Ueig1 = compute_ueig(Weig1,dz); % Normalized Ueig

% Case 2: compute Ueig directly (Early) 
[UeigJ,~,~,~] = imU.ModesAtFrequency(omega); % w-const EVP 
UeigJ(isnan(UeigJ))=0; AA=repmat(sum(UeigJ.^2.*dz,1)./H,[N 1]).^(1/2);
AA(AA==0)=Inf; Ueig2=UeigJ./AA;
Ueig2(:,Ueig2(N,:)<0)=-Ueig2(:,Ueig2(N,:)<0);

% Case 3: Same as (1) but using ufiltint 
%[~,Weig3,~,~] = imWKB.ModesAtFrequency(omega); % w-const EVP 
%Ueig3 = compute_ueig(Weig1,dz); % Normalized Ueig
Ueig3 = Ueig1;

% Case 4: compute Ueig from Weig (Oladeji) 
[~,~,~,~,~,Ueig4] = comp_eigen_ola_old(rho_mean,zf,im.f0,omega); 

% Case 5: compute Ueig directly (Oladeji) 
[~,~,~,~,Ueig5] = comp_eigen_ola_new(rho_mean,zf,im.f0,omega); 

% Case 6: Same as (5) but using ufiltint
Ueig6 = Ueig5;

%% Mode fitting and statistics
% Compute the fit
nmodes = 5; 

% CASE 1
[SS_tot,SS_res,SS_fit,r2u1,s2u1,u1] = compute_modfit2(ufilt,dz,Ueig1,nmodes); 
[SS_tot,SS_res,SS_fit,r2v1,s2v1,v1] = compute_modfit2(vfilt,dz,Ueig1,nmodes); 

% CASE 2
[SS_tot,SS_res,SS_fit,r2u2,s2u2,u2] = compute_modfit2(ufilt,dz,Ueig2,nmodes); 
[SS_tot,SS_res,SS_fit,r2v2,s2v2,v2] = compute_modfit2(vfilt,dz,Ueig2,nmodes); 

% CASE 3
[SS_tot,SS_res,SS_fit,r2u3,s2u3,u3] = compute_modfit2(ufiltint,dz,Ueig3,nmodes); 
[SS_tot,SS_res,SS_fit,r2v3,s2v3,v3] = compute_modfit2(vfiltint,dz,Ueig3,nmodes); 

% CASE 4
[SS_tot,SS_res,SS_fit,r2u4,s2u4,u4] = compute_modfit2(ufilt,dz,Ueig4,nmodes); 
[SS_tot,SS_res,SS_fit,r2v4,s2v4,v4] = compute_modfit2(vfilt,dz,Ueig4,nmodes); 

% CASE 5
[SS_tot,SS_res,SS_fit,r2u5,s2u5,u5] = compute_modfit2(ufilt,dz,Ueig5,nmodes); 
[SS_tot,SS_res,SS_fit,r2v5,s2v5,v5] = compute_modfit2(vfilt,dz,Ueig5,nmodes); 

% CASE 6
[SS_tot,SS_res,SS_fit,r2u6,s2u6,u6] = compute_modfit2(ufiltint,dz,Ueig6,nmodes); 
[SS_tot,SS_res,SS_fit,r2v6,s2v6,v6] = compute_modfit2(vfiltint,dz,Ueig6,nmodes); 

% Fit of all modes (1 to nmodes)
ufit1 = sum(u1,3); vfit1 = sum(v1,3); 
ufit2 = sum(u2,3); vfit2 = sum(v2,3); 
ufit3 = sum(u3,3); vfit3 = sum(v3,3); 
ufit4 = sum(u4,3); vfit4 = sum(v4,3); 
ufit5 = sum(u5,3); vfit5 = sum(v5,3); 
ufit6 = sum(u6,3); vfit6 = sum(v6,3); 

% Coefficient of determination
%[r2u1,s2u1] = compute_r2m(ufilt,u1,dz,H); 
%[r2v1,s2v1] = compute_r2m(vfilt,v1,dz,H); 
%
%[r2u2,s2u2] = compute_r2m(ufilt,u2,dz,H); 
%[r2v2,s2v2] = compute_r2m(vfilt,v2,dz,H); 
%
%[r2u3,s2u3] = compute_r2m(ufilt,u3,dz,H); 
%[r2v3,s2v3] = compute_r2m(vfilt,v3,dz,H); 
%
%[r2u4,s2u4] = compute_r2m(ufilt,u4,dz,H); 
%[r2v4,s2v4] = compute_r2m(vfilt,v4,dz,H); 
%
%[r2u5,s2u5] = compute_r2m(ufilt,u5,dz,H); 
%[r2v5,s2v5] = compute_r2m(vfilt,v5,dz,H); 
%
%[r2u6,s2u6] = compute_r2m(ufilt,u6,dz,H); 
%[r2v6,s2v6] = compute_r2m(vfilt,v6,dz,H); 


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

% Filtered velocities
figure
subplot(411)
plot(t,uiso(end,:),'k'); hold on
plot(t,ufilt(end,:),'b'); 
text(t(10),-0.4,'SURFACE')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('u [m/s]') 
%legend('original','filtered','Location','NorthEast')
subplot(412)
plot(t,viso(end,:),'k'); hold on
plot(t,vfilt(end,:),'b'); 
text(t(10),0.18,'SURFACE')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('v [m/s]') 
subplot(413)
plot(t,uiso(1,:),'k'); hold on
plot(t,ufilt(1,:),'b'); 
text(t(10),0.05,'BOTTOM')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('u [m/s]') 
subplot(414)
plot(t,viso(1,:),'k'); hold on
plot(t,vfilt(1,:),'b'); 
text(t(300),0.03,'BOTTOM')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('v [m/s]') 
print('filtering.png','-r300','-dpng') 

end % if plotini


%% ploteig: eigenfunctions (Ueig, Weig)
if ploteig
figure;
subplot(151)
plot(Ueig3(:,1),zc,'k','LineWidth',2.5); hold on
plot(Ueig1(:,1),zc,'y'); hold on
plot(Ueig2(:,1),zc,'b');
plot(Ueig4(:,1),zc,'--r');
plot(Ueig5(:,1),zc,'--g');
title('Mode 1'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
legend(cname{3},cname{1},cname{2},cname{4},cname{5},'Box','on',...
       'Position',[0.04 0.67 0.06 0.1],'FontSize',5)

subplot(152)
plot(Ueig3(:,2),zc,'k','LineWidth',2.5); hold on
plot(Ueig1(:,2),zc,'y'); 
plot(Ueig2(:,2),zc,'b');
plot(Ueig4(:,2),zc,'--r');
plot(Ueig5(:,2),zc,'--g');
title('Mode 2'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(153)
plot(Ueig3(:,3),zc,'k','LineWidth',2.5); hold on
plot(Ueig1(:,3),zc,'y');
plot(Ueig2(:,3),zc,'b');
plot(Ueig4(:,3),zc,'--r');
plot(Ueig5(:,3),zc,'--g');
title('Mode 3'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(154)
plot(Ueig3(:,4),zc,'k','LineWidth',2.5); hold on
plot(Ueig1(:,4),zc,'y');
plot(Ueig2(:,4),zc,'b');
plot(Ueig4(:,4),zc,'--r');
plot(Ueig5(:,4),zc,'--g');
title('Mode 4'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(155)
plot(Ueig3(:,5),zc,'k','LineWidth',2.5); hold on
plot(Ueig1(:,5),zc,'y'); 
plot(Ueig2(:,5),zc,'b');
plot(Ueig4(:,5),zc,'--r');
plot(Ueig5(:,5),zc,'--g');
title('Mode 5'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

print('Ueig.png','-r300','-dpng')

end % if ploteig

%% plotfit: fitted velocities vs original velocities
if plotfit 

[x,y] = meshgrid(t,zIn); 
cmin = -0.02; cmax = 0.02; 
ymin = -H; 

% Residuals: U (color contours)
figure
subplot(511)
pcolor(x,y,ufilt); shading interp; colorbar; caxis([-0.05 0.05])
title('Band-passed u_{iso}');
datetick('x','dd'); xlim([t(1) t(end)]);
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(512)
pcolor(x,y,ufit1-ufilt); shading interp; colorbar; caxis([cmin cmax])
title(['Residual (u_{fit}-u_{iso}) for ' cname{1}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(513)
pcolor(x,y,ufit2-ufilt); shading interp; colorbar; caxis([cmin cmax])
title(['Residual (u_{fit}-u_{iso}) for ' cname{2}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,ufit4-ufilt); shading interp; colorbar; caxis([cmin cmax])
title(['Residual (u_{fit}-u_{iso}) for ' cname{4}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x,y,ufit5-ufilt); shading interp; colorbar; caxis([cmin cmax])
title(['Residual (u_{fit}-u_{iso}) for ' cname{5}])
datetick('x','dd'); xlim([t(1) t(end)]); xlabel('Days since Sept. 1, 2016')
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

print('ufit_vs_uiso.png','-r300','-dpng')

% V (color contours)
figure
subplot(511)
pcolor(x,y,vfilt); shading interp; colorbar; caxis([cmin cmax])
title('Band-passed V_{iso}'); 
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(512)
pcolor(x,y,vfit1-vfilt); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for ' cname{1}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(513)
pcolor(x,y,vfit2-vfilt); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for ' cname{2}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,vfit4-vfilt); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for ' cname{4}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x,y,vfit5-vfilt); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for ' cname{5}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

print('vfit_vs_viso.png','-r300','-dpng')

% V (color contours)
figure
subplot(511)
pcolor(x,y,vfilt); shading interp; colorbar; caxis([cmin cmax])
title('Band-passed V_{iso}');
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(512)
pcolor(x,y,squeeze(v1(:,:,1))); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for mode 1'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(513)
pcolor(x,y,sum(v1(:,:,1:2),3)); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for modes 1-2'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,sum(v1(:,:,1:3),3)); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} after modes 1-3'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x,y,vfit1); shading interp; colorbar; caxis([cmin cmax])
title(['V_{fit} for modes 1-5'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

print('vfit_check.png','-r300','-dpng')

% Time-series
% surface
figure
subplot(211)
plot(t,ufilt(30,:),'k','LineWidth',2.5); hold on
plot(t,squeeze(u1(30,:,1)),'b');
plot(t,sum(u1(30,:,1:2),3),'r')
plot(t,sum(u1(30,:,1:3),3),'g')
datetick('x','dd'); ylabel('u [m/s]'); pbaspect([5 1 1])
legend('u_{iso}','Mode 1','Modes 1-2','Modes 1-3','Location','EastOutside')
subplot(212)
plot(t,vfilt(30,:),'k','LineWidth',2.5); hold on
plot(t,squeeze(v1(30,:,1)),'b');
plot(t,sum(v1(30,:,1:2),3),'r')
plot(t,sum(v1(30,:,1:3),3),'g')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}','Mode 1','Modes 1-2','Modes 1-3','Location','EastOutside')
print('uv_sur_check.png','-r300','-dpng')

figure 
subplot(211)
plot(t,ufilt(1,:),'k','LineWidth',2.5); hold on
plot(t,ufit1(1,:),'b'); 
plot(t,ufit2(1,:),'r') 
plot(t,ufit4(1,:),'g')
datetick('x','dd'); ylabel('u [m/s]'); pbaspect([5 1 1])
legend('u_{iso}',cname{1},cname{2},cname{4},'Location','EastOutside')
subplot(212)
plot(t,vfilt(1,:),'k','LineWidth',2.5); hold on
plot(t,vfit1(1,:),'b'); 
plot(t,vfit2(1,:),'r') 
plot(t,vfit4(1,:),'g')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}',cname{1},cname{2},cname{4},'Location','EastOutside')
print('uv_bot.png','-r300','-dpng') 

% middepth
figure
subplot(211)
plot(t,ufilt(8,:),'k','LineWidth',2.5); hold on
plot(t,ufit1(8,:),'b'); 
plot(t,ufit2(8,:),'r') 
plot(t,ufit4(8,:),'g')
datetick('x','dd'); ylabel('u [m/s]'); pbaspect([5 1 1])
legend('u_{iso}',cname{1},cname{2},cname{4},'Location','EastOutside')

subplot(212)
plot(t,vfilt(8,:),'k','LineWidth',2.5); hold on
plot(t,vfit1(8,:),'b'); 
plot(t,vfit2(8,:),'r') 
plot(t,vfit4(8,:),'g')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}',cname{1},cname{2},cname{4},'Location','EastOutside')
print('uv_mid.png','-r300','-dpng') 

% bottom
figure
subplot(211)
plot(t,ufilt(end,:),'k','LineWidth',2.5); hold on
plot(t,ufit1(end,:),'b'); 
plot(t,ufit2(end,:),'r') 
plot(t,ufit4(end,:),'g')
datetick('x','dd'); ylabel('u [m/s]');  pbaspect([5 1 1])
legend('u_{iso}',cname{1},cname{2},cname{4},'Location','EastOutside')
subplot(212)
plot(t,vfilt(end,:),'k','LineWidth',2.5); hold on
plot(t,vfit1(end,:),'b'); 
plot(t,vfit2(end,:),'r') 
plot(t,vfit4(end,:),'g')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}',cname{1},cname{2},cname{4},'Location','EastOutside')
print('uv_sur.png','-r300','-dpng') 

end % if plotfit


%% plotsta: coefficient of determination R^2 and S^2
if plotsta

figure
subplot(121)
bar([r2u1 r2u2 r2u3 r2u4 r2u5 r2u6])
title(['R^2 (u) at location ' locstr]); 
xlabel('Modes'); ylabel('R^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',fntsz); xticklabels({'1','1-2','1-3','1-4','1-5'});

subplot(122)
bar([r2v1 r2v2 r2v3 r2v4 r2v5 r2v6])
title(['R^2 (v) at location ' locstr]); 
xlabel('Modes'); ylabel('R^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',fntsz); xticklabels({'1','1-2','1-3','1-4','1-5'});
legend(cname{1},cname{2},cname{3},cname{4},cname{5},cname{6},...
       'Position',[0.16 0.57 0.06 0.05],'FontSize',5,'Box','on')
print('r2.png','-r300','-dpng') 


figure
subplot(121)
bar([s2u1 s2u2 s2u3 s2u4 s2u5 s2u6])
title('SS_{fit}/SS_{tot} for u_{fit}'); xlabel('Modes'); ylabel('S^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',fntsz); xticklabels({'1','1-2','1-3','1-4','1-5'});
%legend(cname{1},cname{2},cname{3},cname{4},cname{5},cname{6},...
%       'Position',[0.16 0.57 0.06 0.05],'FontSize',5,'Box','on')

subplot(122)
bar([s2v1 s2v2 s2v3 s2v4 s2v5 s2v6])
title('SS_{fit}/SS_{tot} for v_{fit}'); xlabel('Modes'); ylabel('S^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',fntsz); xticklabels({'1','1-2','1-3','1-4','1-5'});
legend(cname{1},cname{2},cname{3},cname{4},cname{5},cname{6},...
       'Position',[0.36 0.57 0.06 0.05],'FontSize',5,'Box','on')

print('r2check.png','-r300','-dpng') 


figure
subplot(211)
plot(t,SS_tot(:,2),'k'); hold on 
plot(t,SS_res(:,2),'b'); 
plot(t,SS_fit(:,2),'r'); 
title('SS (u) before fitting mode 3')
legend('SS_{tot}','SS_{res}','SS_{fit}','Location','North')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('SS') 

subplot(212)
plot(t,SS_tot(:,3),'k'); hold on 
plot(t,SS_res(:,3),'b'); 
plot(t,SS_fit(:,3),'r'); 
title('SS (u) after fitting mode 3')
legend('SS_{tot}','SS_{res}','SS_{fit}','Location','North')
pbaspect([5 1 1]); datetick('x','dd'); ylabel('SS') 

print('stat_check.png','-r300','-dpng')

end  % if plotsta


%% Move all figures to dirout
system(['rm ' dirout '*.png']);
system(['mv *.png ' dirout]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
