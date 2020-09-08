clc; clear; close all; 
%%UVFIT_COMP computes harmonic fits (u,v) from HYCOM profiles.
% UVFIT_COMP is used to compare the eigenfunctions computed using 
% Jeffrey Early's InternalModes solver with other finite difference
% (FD) schemes. 
%
% INPUT OPTIONS: (=1 plot, =0 no plot) 
% plotini => plots density, stratification freq. and filtering.
% ploteig => plots horizontal eigenvalues ()
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
addpath /data/msolano/Matlab     % All functions!
addpath /data/msolano/toolbox/GLOceanKit/Matlab/InternalModes
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/BSpline
addpath /data/msolano/toolbox/GLNumericalModelingKit/Matlab/Distributions
addpath /home/mbui/Temp/forMiguel/funcs/

%% Input Options
% Location (=1 North Atlantic, =2 South Pacific)
loc = 3;  locstr = num2str(loc); 

% plotting options (1=yes, 0=no)
plotini = 1; % Stratification and filtering 
ploteig = 1; % Eigenvalues (,Weig) 
plotfit = 1; % Velocity and fit (pcolor and time series)
plotsta = 1; % Statistics (R2 and S2)
fntsz = 6;   % legend font size 

dirout = '/data/msolano/figures/modes_hycom/';

%% Load variables
load(['profile_loc' locstr '.mat'])

uiso = profile.uiso; 
viso = profile.viso; 
ufilt = profile.ufilt; 
vfilt = profile.vfilt; 
ufiltint = profile.ufiltint; 
vfiltint = profile.vfiltint; 

%% Cases
cname = {'WKB','FD (Early)','FD (Early U)','FD (Oladeji)','FD (Oladeji U)','Uniform'};

%% Initialize
% Dimensions and constants
[nz,nt] = size(uiso); 

% Constants
g = 9.81;                        % gravity 
omega = 12.1408331767/24/3600;   % M2 frequency 

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
%rho_mean = sort(profile.rho_mean,'descend')+1000;  % Mean density 
rho_mean = profile.rho_mean+1000;  % Mean sea-water density 
lat = profile.latitude;  % Latitude 
rhof_mean = interp1(zc,rho_mean,zf,'linear','extrap');

% Initilized Early's InternalModes
imWKB = InternalModes(rhof_mean,zf,zf,lat);
imf = InternalModes(rhof_mean,zf,zf,lat,'method','finiteDifference');
imU = InternalModes(rho_mean,zc,zc,lat,'method','finiteDifference');

f = imU.f0; 

N2 = imU.N2; 
N2(end)=N2(end-1);
imU2 = InternalModes(N2,zc,zc,lat,'method','finiteDifference',...
                     'N2',1,'rho0',imf.rho0);

imU.N2-imU2.N2


%% Mask N2 and re-compute rho_mean
% Method 1 (N2 at faces, rho at centers)
dzc = diff(zc); 
[N2mod,rho_meanmod] = N2mask(rho_mean,imf.N2,dzc);
rhof_meanmod = interp1(zc,rho_meanmod,zf,'linear','extrap');
im = InternalModes(rhof_meanmod,zf,zf,lat);
%
%% Method 2 (N2 at centers, rho at faces) 
%[N2mod2,rhof_meanmod2] = N2maskf(rhof_mean,imU.N2,dz);
%rho_meanmod2 = interp1(zf,rhof_meanmod2,zc,'linear','extrap');
%im2 = InternalModes(rhof_meanmod2,zf,zf,lat);
%
%  
figure
subplot(121)
plot(rho_mean,zc,'k'); hold on 
plot(rho_meanmod,zc,'b');
title('Density'); xlabel('\rho [kg/m^3]'); ylabel('Depth [m]')
legend('Original','Modified')
ylim([-H 0]);
subplot(122)
semilogx(imf.N2,zf,'k'); hold on
semilogx(N2mod,zf,'b');
semilogx(im.N2,zf,'r');
title('Stratification'); xlabel('N^2'); ylabel('Depth [m]')
ylim([-H 0]);
legend('Original','Masked','Early','Location','SouthWest')
print('stratification1.png','-r300','-dpng')
%
%figure
%subplot(121)
%plot(rho_mean,zc,'k'); hold on
%plot(rhof_meanmod2,zf,'b');
%title('Density'); xlabel('\rho [kg/m^3]'); ylabel('Depth [m]')
%legend('Original','Modified')
%ylim([-H 0]);
%subplot(122)
%semilogx(imf.N2,zf,'k'); hold on
%semilogx(N2mod2,zc,'b');
%semilogx(im2.N2,zf,'r');
%title('Stratification'); xlabel('N^2'); ylabel('Depth [m]')
%ylim([-H 0]);
%legend('Original','Masked','Early','Location','SouthWest')
%print('stratification2.png','-r300','-dpng')
%
%system(['rm ' dirout '*.png']);
%system(['mv *.png ' dirout]);


%% Compute eigenmodes
% *All cases are interpolated to time-mean layer
% *All cases (except 6) are on native grid (HYCOM)

%% Case 1: WKB solution (Early)
% Eigenfunctions
[~,Weig1,h1,k1] = imWKB.ModesAtFrequency(omega); % omega-const EVP 
Ueig1 = compute_ueig(Weig1,dz);                  % Normalized Ueig

% Eigenvalues
L1 = 2*pi./k1;                                   % wave-length 
C1 = omega./k1;                                  % phase-speed
Cg1 = (omega^2-f^2)./(omega*k1);                 % group-speed 

%% Case 2: FD Weig (Early) 
% Eigenfunctions
[~,Weig2,h2,k2] = imf.ModesAtFrequency(omega);   % omega-const EVP
Ueig2 = compute_ueig(Weig2,dz);                  % Normalized Ueig

% Eigenvalues
L2 = 2*pi./k2;                                   % wave-length 
C2 = omega./k2;                                  % phase-speed
Cg2 = (omega^2-f^2)./(omega*k2);                 % group-speed 
%Ce2 = sqrt(omega^2-f^2)./(k2.^2);                % eigen-speed
Ce2 = sqrt(g*h2);

%% Case 3: FD  (Early) 
% Eigenfunctions
imU.normalization = 'wMax';
[t,~,h3,k3] = imU.ModesAtFrequency(omega);   % omega-const EVP

dzu = repmat(dz,[1,size(dz,1)]);
AA = repmat(sum(t.^2.*dzu,1)./H,[N 1]).^(1/2);
AA(AA==0) = Inf;
 = Ueigt./AA;
(:,Ueig(N,:)<0) = -Ueig(:,Ueig(N,:)<0);
Ueig3 = Ueig;

% Eigenvalues
L3 = 2*pi./k3;                                   % wave-length 
C3 = omega./k3;                                  % phase-speed
Cg3 = (omega^2-f^2)./(omega*k3);                 % group-speed 

%% Case 4: FD Weig (Oladeji) 
% Eigenfunctions
S = compute_eigen(rho_mean,zf,f,omega); 
Ueig4 = S.Ueig;

% Eigenvalues
L4 = S.L; 
C4 = S.C;
Cg4 = S.Cg;
k4 = S.k;
Ce4 = sqrt(omega^2-f^2)./(k4.^2);                % eigen-speed

%% Case 5: FD  (Oladeji) 
% Eigenfunctions
S = compute_eigenU(rho_mean,zf,f,omega); 
Ueig5 = S.Ueig;

% Eigenvalues
L5 = S.L; 
C5 = S.C;
Cg5 = S.Cg;

%% Case 6: Same as (5) but using ufiltint
%nzM = 100; % number of layers (uniform)
dzi = 25; 
nzM = round(H/dzi); 
zfM = linspace(-H,0,nzM+1)'; 
dz1 = diff(zfM); 
zcM = zfM(1:end-1) + dz1/2;
%N2M = interp1(zf,imf.N2,zfM,'linear','extrap'); 
N2M = interp1(zc,S.N2,zfM,'linear','extrap'); 
[C6,Cg6,L6,Weig6,Ueig6t] = sturm_liouville_hyd_normalize(omega,sqrt(N2M),dz1(1),f);
%[C6,Cg6,L6,Weig6,Ueig6t] = sturm_liouville_hyd_normalize(omega,N2M,dz1(1),f);
Ueig6t(:,2:2:end) = -Ueig6t(:,2:2:end);
Ueig6 = Ueig6t(:,2:end); 

% Interpolate to uniform grid
ufiltM = zeros(nzM,nt); 
vfiltM = zeros(nzM,nt); 
for i = 1:nt
    ufiltM(:,i) = interp1(profile.zc(:,i),profile.ufilt(:,i),zcM,'linear','extrap');  
    vfiltM(:,i) = interp1(profile.zc(:,i),profile.vfilt(:,i),zcM,'linear','extrap');  
end

%% Print some info 
fprintf('\nSolutions:')

fprintf('\nMode 1 wavelength (L)')
fprintf('\nOladeji = %4.2f[km]',L4(1)/1000)
fprintf('\nJeffrey = %4.2f[km]',L3(1)/1000)
fprintf('\nMaarten = %4.2f[km]\n',L6(1)/1000)

fprintf('\nMode 1 group-speed (Cg)')
fprintf('\nOladeji = %4.4f[m/s]',Cg4(1))
fprintf('\nJeffrey = %4.4f[m/s]',Cg3(1))
fprintf('\nMaarten = %4.4f[m/s]\n',Cg6(1))

fprintf('\nMode 1 phase-speed (C)')
fprintf('\nOladeji = %4.4f[m/s]',C4(1))
fprintf('\nJeffrey = %4.4f[m/s]',C3(1))
fprintf('\nMaarten = %4.4f[m/s]\n',C6(1))

fprintf('\nMode 1 eigen-speed (C)')
fprintf('\nOladeji = %4.4f[m/s]',Ce4(1))
fprintf('\nJeffrey = %4.4f[m/s]\n',Ce2(1))


%% Mode fitting and statistics
% Compute the fit
nmodes = 5; 

% CASE 1
[r2u1,u1] = compute_modfit(ufiltint,dz,Ueig1,nmodes); 
[r2v1,v1] = compute_modfit(vfiltint,dz,Ueig1,nmodes); 

% CASE 2
[r2u2,u2] = compute_modfit(ufiltint,dz,Ueig2,nmodes); 
[r2v2,v2] = compute_modfit(vfiltint,dz,Ueig2,nmodes); 

% CASE 3
[r2u3,u3] = compute_modfit(ufiltint,dz,Ueig3,nmodes); 
[r2v3,v3] = compute_modfit(vfiltint,dz,Ueig3,nmodes); 

% CASE 4
[r2u4,u4] = compute_modfit(ufiltint,dz,Ueig4,nmodes); 
[r2v4,v4] = compute_modfit(vfiltint,dz,Ueig4,nmodes); 

% CASE 5
[r2u5,u5] = compute_modfit(ufiltint,dz,Ueig5,nmodes); 
[r2v5,v5] = compute_modfit(vfiltint,dz,Ueig5,nmodes); 

% CASE 6
[r2u6,u6] = compute_modfit(ufiltM,dz1,Ueig6,nmodes); 
[r2v6,v6] = compute_modfit(vfiltM,dz1,Ueig6,nmodes); 

% Fit of all modes (1 to nmodes)
ufit1 = sum(u1,3); vfit1 = sum(v1,3); 
ufit2 = sum(u2,3); vfit2 = sum(v2,3); 
ufit3 = sum(u3,3); vfit3 = sum(v3,3); 
ufit4 = sum(u4,3); vfit4 = sum(v4,3); 
ufit5 = sum(u5,3); vfit5 = sum(v5,3); 
ufit6 = sum(u6,3); vfit6 = sum(v6,3); 


%% PLOTS
ymin = -H; 

%% plotini: stratification, filtered velocities
if plotini

% Stratification
figure
subplot(121)
plot(imf.rho+1000,imf.z,'r'); 
title('Density'); xlabel('\rho [kg/m^3]'); ylabel('Depth [m]')
ylim([-H 0]);
subplot(122)
semilogx(imf.N2,imf.z,'k'); hold on
semilogx(S.N2,zc,'b'); 
title('Stratification'); xlabel('N [1/sec]'); ylabel('Depth [m]')
ylim([-500 0]);
legend('Early (im.N2)','Oladeji (S.N2)','Location','SouthEast')
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


%% ploteig: eigenfunctions (, Weig)
if ploteig
%ymin = -100;
%figure;
%subplot(161)
%semilogx(sqrt(imf.N2),imf.z,'k'); hold on
%semilogx(sqrt(imU.N2),imU.z,'b'); 
%%plot(imf.N2,imf.z,'k'); hold on
%%semilogx(sqrt(N2M),zfM,'b'); hold on
%title('Stratification'); xlabel('N^2 [1/s^2]'); ylabel('Depth [m]')
%axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]); 
%legend('no mask','with mask','Box','on',...
%       'Position',[0.045 0.27 0.06 0.1],'FontSize',5)
%
%subplot(162)
%%plot(1(:,1),zc,'k','LineWidth',2.5); hold on
%plot(2(:,1)-Ueig3(:,1),zc,'b'); hold on
%plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
%title('Mode 1'); xlabel('U_{eig}');% ylabel('Depth [m]');
%yticklabels('off')
%axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
%%legend(cname{1},cname{2},cname{4},cname{6},'Box','on',...
%
%subplot(163)
%%plot(1(:,2),zc,'k','LineWidth',2.5); hold on
%plot(2(:,2)-Ueig3(:,2),zc,'b'); hold on
%plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
%title('Mode 2'); xlabel('U_{eig}'); %ylabel('Depth [m]');
%axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
%yticklabels([])
%
%subplot(164)
%%plot(1(:,3),zc,'k','LineWidth',2.5); hold on
%plot(2(:,3)-Ueig3(:,3),zc,'b'); hold on
%plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
%title('Mode 3'); xlabel('U_{eig}'); %ylabel('Depth [m]');
%axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
%yticklabels([])
%
%subplot(165)
%%plot(1(:,4),zc,'k','LineWidth',2.5); hold on
%plot(2(:,4)-Ueig3(:,4),zc,'b'); hold on
%plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
%title('Mode 4'); xlabel('U_{eig}');% ylabel('Depth [m]');
%axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
%yticklabels([])
%
%subplot(166)
%%plot(1(:,5),zc,'k','LineWidth',2.5); hold on
%plot(2(:,5)-Ueig3(:,5),zc,'b'); hold on
%plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
%title('Mode 5'); xlabel('U_{eig}');% ylabel('Depth [m]');
%axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
%yticklabels([])

figure;
subplot(161)
semilogx(imf.N2,imf.z,'k'); hold on
semilogx(N2M,zfM,'b'); 
%plot(imf.N2,imf.z,'k'); hold on
%semilogx(sqrt(N2M),zfM,'b'); hold on
title('Stratification'); xlabel('N^2 [1/s^2]'); ylabel('Depth [m]')
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
ylim([ymin 0]); yticklabels([])

subplot(162)
plot(Ueig1(:,1),zc,'k','LineWidth',2.5); hold on
plot(Ueig2(:,1),zc,'b'); hold on
plot(Ueig3(:,1),zc,'r');
plot(Ueig4(:,1),zc,'--r');
plot(Ueig6(:,1),zcM,'--g');
plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
title('Mode 1'); xlabel('U_{eig}');% ylabel('Depth [m]');
yticklabels('off')
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
legend(cname{1},cname{2},cname{4},cname{6},'Box','on',...
       'Position',[0.035 0.67 0.06 0.1],'FontSize',5)

subplot(163)
plot(Ueig1(:,2),zc,'k','LineWidth',2.5); hold on
plot(Ueig2(:,2),zc,'b'); 
plot(Ueig3(:,2),zc,'r');
plot(Ueig4(:,2),zc,'--r');
plot(Ueig6(:,2),zcM,'--g');
plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
title('Mode 2'); xlabel('U_{eig}'); %ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
yticklabels([])

subplot(164)
plot(Ueig1(:,3),zc,'k','LineWidth',2.5); hold on
plot(Ueig2(:,3),zc,'b');
plot(Ueig3(:,3),zc,'r');
plot(Ueig4(:,3),zc,'--r');
plot(Ueig6(:,3),zcM,'--g');
plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
title('Mode 3'); xlabel('U_{eig}'); %ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
yticklabels([])

subplot(165)
plot(Ueig1(:,4),zc,'k','LineWidth',2.5); hold on
plot(Ueig2(:,4),zc,'b');
plot(Ueig3(:,4),zc,'r');
plot(Ueig4(:,4),zc,'--r');
plot(Ueig6(:,4),zcM,'--g');
plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
title('Mode 4'); xlabel('U_{eig}');% ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
yticklabels([])

subplot(166)
plot(Ueig1(:,5),zc,'k','LineWidth',2.5); hold on
plot(Ueig2(:,5),zc,'b'); 
plot(Ueig3(:,5),zc,'r');
plot(Ueig4(:,5),zc,'--r');
plot(Ueig6(:,5),zcM,'--g');
plot(zeros(2,1),[-H 0],'--k','LineWidth',1.0);
title('Mode 5'); xlabel('U_{eig}');% ylabel('Depth [m]');
axis tight; ylim([ymin 0]); set(gca,'FontSize',6); pbaspect([1 4 1]);
yticklabels([])

print('.png','-r300','-dpng')

end % if ploteig

%% plotfit: fitted velocities vs original velocities
if plotfit 

[x,y] = meshgrid(t,zc); 
[x1,y1] = meshgrid(t,zcM); 
cmin = -0.04; cmax = 0.04; 
ymin = -H; 

% Residuals: U (color contours)
figure; colormap('jet')
subplot(511)
pcolor(x,y,ufilt); shading flat; colorbar; caxis([cmin cmax])
title('Band-passed eastward velocity (u_{iso})');
datetick('x','dd'); xlim([t(1) t(end)]);
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);
set(gca,'xticklabel',[],'TickLength',[0.01 0.01]);

subplot(512)
pcolor(x,y,ufit1); shading flat; colorbar; caxis([cmin cmax])
title([cname{1}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);
set(gca,'xticklabel',[],'TickLength',[0.01 0.01]);

subplot(513)
pcolor(x,y,ufit2); shading flat; colorbar; caxis([cmin cmax])
title([cname{2}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);
set(gca,'xticklabel',[],'TickLength',[0.01 0.01]);

subplot(514)
pcolor(x,y,ufit4); shading flat; colorbar; caxis([cmin cmax])
title([cname{4}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);
set(gca,'xticklabel',[],'TickLength',[0.01 0.01]);

subplot(515)
pcolor(x1,y1,ufit6); shading flat; colorbar; caxis([cmin cmax])
title([cname{6}])
datetick('x','dd'); xlim([t(1) t(end)]); 
xlabel('Days since Sept. 1, 2016','FontSize',8)
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);
set(gca,'TickLength',[0.01 0.01],'TickDir','both');

print('ufit_vs_uiso.png','-r300','-dpng')

% V (color contours)
figure; colormap('jet')
subplot(511)
pcolor(x,y,vfilt); shading flat; colorbar; caxis([cmin cmax])
title('Band-passed northward velocity (v_{iso})'); 
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(512)
pcolor(x,y,vfit1); shading flat; colorbar; caxis([cmin cmax])
title([cname{1}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(513)
pcolor(x,y,vfit2); shading flat; colorbar; caxis([cmin cmax])
title([cname{2}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,vfit4); shading flat; colorbar; caxis([cmin cmax])
title([cname{4}])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x1,y1,vfit6); shading flat; colorbar; caxis([cmin cmax])
title([cname{6}])
datetick('x','dd'); xlim([t(1) t(end)])
xlabel('Days since Sept. 1, 2016','FontSize',8)
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

print('vfit_vs_viso.png','-r300','-dpng')

% V (color contours)
figure; colormap('jet')
subplot(511)
pcolor(x,y,vfilt); shading flat; colorbar; caxis([cmin cmax])
title('Band-passed V_{iso}');
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(512)
pcolor(x,y,squeeze(v1(:,:,1))); shading flat; colorbar; caxis([cmin cmax])
title(['V_{fit} for mode 1'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(513)
pcolor(x,y,sum(v1(:,:,1:2),3)); shading flat; colorbar; caxis([cmin cmax])
title(['V_{fit} for modes 1-2'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(514)
pcolor(x,y,sum(v1(:,:,1:3),3)); shading flat; colorbar; caxis([cmin cmax])
title(['V_{fit} after modes 1-3'])
datetick('x','dd'); xlim([t(1) t(end)])
ylabel('Depth [m]'); ylim([ymin 0]); set(gca,'FontSize',fntsz);

subplot(515)
pcolor(x,y,vfit1); shading flat; colorbar; caxis([cmin cmax])
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
plot(t,ufit5(1,:),'m')
datetick('x','dd'); ylabel('u [m/s]'); pbaspect([5 1 1])
legend('u_{iso}',cname{1},cname{2},cname{4},cname{5},'Location','EastOutside')
subplot(212)
plot(t,vfilt(1,:),'k','LineWidth',2.5); hold on
plot(t,vfit1(1,:),'b'); 
plot(t,vfit2(1,:),'r') 
plot(t,vfit4(1,:),'g')
plot(t,vfit5(1,:),'m')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}',cname{1},cname{2},cname{4},cname{5},'Location','EastOutside')
print('uv_bot.png','-r300','-dpng') 

% middepth
figure
subplot(211)
plot(t,ufilt(8,:),'k','LineWidth',2.5); hold on
plot(t,ufit1(8,:),'b'); 
plot(t,ufit2(8,:),'r') 
plot(t,ufit4(8,:),'g')
plot(t,ufit5(8,:),'m')
datetick('x','dd'); ylabel('u [m/s]'); pbaspect([5 1 1])
legend('u_{iso}',cname{1},cname{2},cname{4},cname{5},'Location','EastOutside')

subplot(212)
plot(t,vfilt(8,:),'k','LineWidth',2.5); hold on
plot(t,vfit1(8,:),'b'); 
plot(t,vfit2(8,:),'r') 
plot(t,vfit4(8,:),'g')
plot(t,vfit5(8,:),'m')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}',cname{1},cname{2},cname{4},cname{5},'Location','EastOutside')
print('uv_mid.png','-r300','-dpng') 

% bottom
figure
subplot(211)
plot(t,ufilt(end,:),'k','LineWidth',2.5); hold on
plot(t,ufit1(end,:),'b'); 
plot(t,ufit2(end,:),'r') 
plot(t,ufit4(end,:),'g')
plot(t,ufit5(end,:),'m')
datetick('x','dd'); ylabel('u [m/s]');  pbaspect([5 1 1])
legend('u_{iso}',cname{1},cname{2},cname{4},cname{5},'Location','EastOutside')
subplot(212)
plot(t,vfilt(end,:),'k','LineWidth',2.5); hold on
plot(t,vfit1(end,:),'b'); 
plot(t,vfit2(end,:),'r') 
plot(t,vfit4(end,:),'g')
plot(t,vfit5(end,:),'m')
datetick('x','dd'); ylabel('v [m/s]'); pbaspect([5 1 1])
legend('v_{iso}',cname{1},cname{2},cname{4},cname{5},'Location','EastOutside')
print('uv_sur.png','-r300','-dpng') 

end % if plotfit


%% plotsta: coefficient of determination R^2 and S^2
if plotsta

%figure
%plot(t,SS_tot(:,5),'k','LineWidth',2); hold on 
%plot(t,SS_res1(:,5),'b')
%plot(t,SS_res2(:,5),'r')
%plot(t,SS_res4(:,5),'g')
%plot(t,SS_res6(:,5),'m')
%datetick('x','dd'); ylabel('SS'); pbaspect([5 1 1])
%legend('SS_{total}',cname{1},cname{2},cname{4},cname{6})
%print('SS.png','-r300','-dpng') 


figure
subplot(121)
bar([r2u1 r2u2 r2u4 r2u6])
title(['R^2 (u) at location ' locstr]); 
xlabel('Modes'); ylabel('R^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',fntsz); xticklabels({'1','1-2','1-3','1-4','1-5'});
legend(cname{1},cname{2},cname{4},cname{5},cname{6},...
       'Position',[0.17 0.59 0.06 0.05],'FontSize',5,'Box','off')

subplot(122)
bar([r2v1 r2v2 r2v4 r2v6])
title(['R^2 (v) at location ' locstr]); 
xlabel('Modes'); ylabel('R^2'); pbaspect([1.5 1 1]) 
set(gca,'FontSize',fntsz); xticklabels({'1','1-2','1-3','1-4','1-5'});
legend(cname{1},cname{2},cname{4},cname{6},...
       'Position',[0.6 0.59 0.06 0.05],'FontSize',5,'Box','off')
print('r2.png','-r300','-dpng') 


end  % if plotsta


%% Move all figures to dirout
system(['rm ' dirout '*.png']);
system(['mv *.png ' dirout]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
