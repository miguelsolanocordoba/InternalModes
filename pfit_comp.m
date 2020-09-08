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

dirout = '/data/msolano/figures/modes_hycom/';

% plotting options (1=yes, 0=no)
plotini = 1; % Stratification and filtering 
ploteig = 1; % Eigenvalues (Ueig,Weig) 
plotfit = 1; % Velocity and fit (pcolor and time series)
plotsta = 1; % Statistics (R2 and S2)
fntsz = 6;   % legend font size 


%% Load variables
load profile_loc2.mat     

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
cname = {'Early W','Early U','Early Int','Ola W','Ola U','Ola Int'};


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
im = InternalModes(rhof_mean,zOut,zOut,lat,'method','finiteDifference',...
                   'orderOfAccuracy',4)
imU = InternalModes(rho_mean,zIn,zIn,lat,'method','finiteDifference',...)
                   'orderOfAccuracy',2)

% Compute perturbation pressure
pert = comp_pertpress(profile.rho,dz); 
%pert = compute_pertpress(profile.rho,profile.zc,profile.zc_mean,profile.zf_mean); 
pertp = filter_pertp(pert); 
%pertpint = comp_pertpressint(profile.rho,diff(profile.zc)); 

%clear profile

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
Ueig3 = Ueig1;

% Case 4: compute Ueig from Weig (Oladeji) 
[~,Ueig4] = comp_eigen_ola(dz,im.N2,im.f0,omega); 
%[~,~,~,~,Ueig4] = comp_eigen_ola_old(rho_mean,zf,im.f0,omega); 

% Case 5: compute Ueig directly (Oladeji) 
[~,~,~,~,Ueig5] = comp_eigen_ola2(dz,im.N2,im.f0,omega,0); 
%[~,~,~,~,Ueig5] = comp_eigen_ola_new(rho_mean,zf,im.f0,omega); 

% Case 6: Same as (5) but using ufiltint
Ueig6 = Ueig5;


%% Mode fitting and statistics
% Compute the fit
nmodes = 5; 

% CASE 1
[~,~,~,r2u1,s2u1,p1] = compute_modfit2(pertp,dz,Ueig1,nmodes); 

% CASE 2
[~,~,~,r2u2,s2u2,p2] = compute_modfit2(pertp,dz,Ueig2,nmodes); 

% CASE 3
[~,~,~,r2u3,s2u3,p3] = compute_modfit2(pertp,dz,Ueig3,nmodes); 

% CASE 4
[~,~,~,r2u4,s2u4,p4] = compute_modfit2(pertp,dz,Ueig4,nmodes); 

% CASE 5
[~,~,~,r2u5,s2u5,p5] = compute_modfit2(pertp,dz,Ueig5,nmodes); 

% CASE 6
[~,~,~,r2u6,s2u6,p6] = compute_modfit2(pertp,dz,Ueig6,nmodes); 

% Fit of all modes (1 to nmodes)
pfit1 = sum(p1,3); 
pfit2 = sum(p2,3); 
%pfit3 = sum(p3,3); 
pfit4 = sum(p4,3); 
pfit5 = sum(p5,3); 
%pfit6 = sum(p6,3); 

% Coefficient of determination
%[r2u1,s2u1] = compute_r2m(pertp,p1,dz,H); 

%[r2u2,s2u2] = compute_r2m(pertp,p2,dz,H); 

%[r2u3,s2u3] = compute_r2m(pertp,p3,dz,H); 

%[r2u4,s2u4] = compute_r2m(pertp,p4,dz,H); 

%[r2u5,s2u5] = compute_r2m(pertp,p5,dz,H); 

%[r2u6,s2u6] = compute_r2m(pertp,p6,dz,H); 



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
plot(Ueig4(:,1),zc,'--r');
plot(Ueig5(:,1),zc,'--g');
title('Mode 1'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
legend(cname{1},cname{2},cname{4},cname{5},'Box','on',...
       'Position',[0.04 0.67 0.06 0.1],'FontSize',5)

subplot(152)
plot(Ueig1(:,2),zc,'k'); hold on
plot(Ueig2(:,2),zc,'b');
plot(Ueig4(:,2),zc,'--r');
plot(Ueig5(:,2),zc,'--g');
title('Mode 2'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(153)
plot(Ueig1(:,3),zc,'k'); hold on
plot(Ueig2(:,3),zc,'b');
plot(Ueig4(:,3),zc,'--r');
plot(Ueig5(:,3),zc,'--g');
title('Mode 3'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(154)
plot(Ueig1(:,4),zc,'k'); hold on
plot(Ueig2(:,4),zc,'b');
plot(Ueig4(:,4),zc,'--r');
plot(Ueig5(:,4),zc,'--g');
title('Mode 4'); xlabel('U_{eig}'); ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);

subplot(155)
plot(Ueig1(:,5),zc,'k'); hold on
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
pcolor(x,y,pfit2); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{2}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,pfit4); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{4}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x,y,pfit5); shading interp; colorbar; caxis([cmin cmax])
title(['P_{fit} for ' cname{5}])
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
bar([r2u1 r2u2 r2u4 r2u5])
title('R^2 for P_{fit}'); xlabel('Modes'); ylabel('R^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',8); xticklabels({'1','1-2','1-3','1-4','1-5'})
legend(cname{1},cname{2},cname{4},cname{5},'Location','NorthWest')

print('r2p.png','-r300','-dpng') 

end  % if plotsta


%% Move all figures to dirout
system(['rm ' dirout '*.png']);
system(['mv *.png ' dirout]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
