clear; clc; close all
%%NISKINE_SAVEPROF saves HYCOM profiles for use with InternalModes
% 
% Created: June 22, 2020 by M. Solano 

% Format 
IEEE = 'ieee-be';
addpath /data/msolano/Matlab

%% Experiment and tile number 
% loc1 > runnum=221;  blki=27; blkj=45; a=100; b=75;
% loc2 > runnum=190;  blki=15; blkj=25; a=40;  b=40;
% loc2 > runnum=190;  blki=15; blkj=25; a=40;  b=40;

runnum = 190;
runnumstr = num2str(runnum);
blki=15;
blkj=25;

% Output point 
a = 40; b = 40; 

% Directories
%dirin = '/data2/mbui/for_keshav/tiles/'; % loc1
dirin = '/data2/msolano/forEmanuel/hycom/GLBc0.04/expt_19.0/'; 
dirout = '/data/msolano/forOladeji/';
saveflag = 0;

fprintf('*** Running niskine_saveprof ***\n')
fprintf('Input directory: %s\n',dirin)
fprintf('Output directory: %s\n',dirout)
fprintf('iTile = %d\n',blki)
fprintf('jTile = %d\n',blkj)

% Grid file data
depfile = [dirin 'griddata/depth_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
lonfile = [dirin 'griddata/plon_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
latfile = [dirin 'griddata/plat_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];

% Variables 
fname1 = [dirin 'u_iso/u_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
fname2 = [dirin 'v_iso/v_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
fname3 = [dirin 'thknss/thknss_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
fname4 = [dirin 'sig/sig_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];

% Dimensions
nbf = 3;  % halo/padding 

nx=150; ny=200; nz=41; nt=624;
nxb=nx+nbf*2;
nyb=ny+nbf*2;

lenrec2 = nxb*nyb+2;
t = datenum(2016,9,1):datenum(0,0,0,1,0,0):datenum(2016,9,15);
nt = numel(t); 

% Load grid
fiddep = fopen(depfile,'r',IEEE);
fidlon = fopen(lonfile,'r',IEEE);
fidlat = fopen(latfile,'r',IEEE);

depdata = fread(fiddep,lenrec2,'single');
londata = fread(fidlon,lenrec2,'single');
latdata = fread(fidlat,lenrec2,'single');

lon1 = []; lat1 = []; depth1 = [];

depth1(:,:) = permute(reshape(depdata(2:end-1),[nxb nyb]),[2 1]);
lon1(:,:) = permute(reshape(londata(2:end-1),[nxb nyb]),[2 1]);
lat1(:,:) = permute(reshape(latdata(2:end-1),[nxb nyb]),[2 1]);

fclose(fiddep);
fclose(fidlon);
fclose(fidlat);

%% Read variables
% Open files
fid1 = fopen(fname1,'r',IEEE);
fid2 = fopen(fname2,'r',IEEE);
fid3 = fopen(fname3,'r',IEEE);
fid4 = fopen(fname4,'r',IEEE);

%% Read variables: u_iso, v_iso, sig, thknss
uiso1 = [];
viso1 = [];
thknss1 = [];
sig1 = [];

% extract layer thickness, u,v in space and time
fprintf('\nExtracing variables...\n') 
for i=1:nt
    
    fprintf('%d/%d\n',i,nt)	
    
    for k=1:nz
        alldata1 = fread(fid1,lenrec2,'single');
        alldata2 = fread(fid2,lenrec2,'single');
        alldata3 = fread(fid3,lenrec2,'single');
        alldata4 = fread(fid4,lenrec2,'single');
        
        uiso1(:,:,k,i)   = permute(reshape(alldata1(2:end-1),[nxb nyb]),[2 1]);
        viso1(:,:,k,i)   = permute(reshape(alldata2(2:end-1),[nxb nyb]),[2 1]);
        thknss1(:,:,k,i) = permute(reshape(alldata3(2:end-1),[nxb nyb]),[2 1]);
        sig1(:,:,k,i)    = permute(reshape(alldata4(2:end-1),[nxb nyb]),[2 1]);
        
    end
end
    
fprintf('\nClosing files...\n')
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);


%% Save output for use in InternalModes
fprintf('\nPreperaring output...\n') 
% Dimensions
[nx,ny,nz,nt] = size(uiso1); 


% Grid data 
depth = depth1(a,b); clear depth1
lon = lon1(a,b);     clear lon1
lat = lat1(a,b);     clear lat1

% Compute cell centers(zc) and faces(zf)
zf1 = zeros(nz+1,nt); 
zc1 = zeros(nz,nt); 

for i = 1:nt
    zf1(2:end,i) = squeeze(cumsum(thknss1(a,b,:,i),3)); 
    zc1(:,i) = 0.5.*(zf1(1:nz,i) + zf1(2:nz+1,i));
end

% if layer thickness < 0.1m, do not save! 
thk0=diff(zf1(:,1))>0.1;
nz = sum(thk0);  % number of layers > 0.1m

zf = -flipud(zf1(1:nz+1,:)); clear zf1;  
zc = -flipud(zc1(1:nz,:)); clear zc1

zf_mean = mean(zf,2); % time-mean cell face
zc_mean = mean(zc,2); % time-mean cell center

% Save main variables (u_iso,v_iso,zf,zc)
uiso = flipud(squeeze(uiso1(a,b,1:nz,:)));            % Velocity (u) 
viso = flipud(squeeze(viso1(a,b,1:nz,:)));            % Velocity (v) 
rho = flipud(squeeze(sig1(a,b,1:nz,:)));              % Density 
rho_mean = flipud(squeeze(mean(sig1(a,b,1:nz,:),4))); % Time-mean density 

% Filter and interpolate baroclinic velocities
[ufilt,vfilt] = filter_uviso(uiso,viso); % Band-pass velocities

ufiltint = zeros(size(ufilt)); 
vfiltint = zeros(size(vfilt)); 
for i=1:nt 
    ufiltint(:,i) = interp1(zc(:,i),ufilt(:,i),zc_mean,'linear','extrap'); 
    vfiltint(:,i) = interp1(zc(:,i),vfilt(:,i),zc_mean,'linear','extrap'); 
end


% Save into profile structure 
fprintf('\nSaving profile...\n') 
profile = struct('latitude',lat,'longitude',lon,'depth',depth,...
                 'uiso',uiso,'viso',viso,'ufilt',ufilt,'vfilt',vfilt,...
		 'ufiltint',ufiltint,'vfiltint',vfiltint,'rho',rho,...
		 'zc',zc,'zf_mean',zf_mean,...
                 'zc_mean',zc_mean,'rho_mean',rho_mean)

save([dirout 'profile_loc2.mat'],'profile','-v7.3');

fprintf('SUCCESS!!!\nProfile saved to: %s\n',dirout)
