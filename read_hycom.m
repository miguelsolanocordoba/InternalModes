function hycom = read_hycom()
%%READ_HYCOM reads in HYCOM's output (.BinF)
% HYCOM = READ_HYCOM reads HYCOM binaries (.BinF) and saves output 
% into a Matlab structure:
%
% hycom.time  % time (in datenum format)
% hycom.lon   % longitude 
% hycom.lat   % latitude 
% hycom.h     % depth 
% hycom.dz    % layer thickness
% hycom.uiso  % baroclinic velocity (u) 
% hycom.viso  % baroclinic velocity (v) 
% hycom.rho   % density  
% 
% Created: July 10, 2020 by M. Solano 

% Format 
IEEE = 'ieee-be';
addpath /data/msolano/Matlab

%% Experiment and tile number 
% loc1 > runnum=221;  blki=27; blkj=45;
% loc2 > runnum=190;  blki=15; blkj=25; 

runnum = 221;
runnumstr = num2str(runnum);
blki=27;
blkj=45;

% Directories
dirin = '/data2/mbui/for_keshav/tiles/'; % loc1
%dirin = '/data2/msolano/forEmanuel/hycom/GLBc0.04/expt_19.0/'; 

fprintf('\nReading HYCOM files (read_hycom)\n')
fprintf('Input directory: %s\n',dirin)
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
fname5 = [dirin 'temp/T_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
fname6 = [dirin 'sal/S_' num2str(runnum) '_blk_' ...
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

lon = []; lat = []; depth = [];

depth(:,:) = permute(reshape(depdata(2:end-1),[nxb nyb]),[2 1]);
lon(:,:) = permute(reshape(londata(2:end-1),[nxb nyb]),[2 1]);
lat(:,:) = permute(reshape(latdata(2:end-1),[nxb nyb]),[2 1]);

fclose(fiddep);
fclose(fidlon);
fclose(fidlat);

%% Read variables
% Open files
fid1 = fopen(fname1,'r',IEEE);
fid2 = fopen(fname2,'r',IEEE);
fid3 = fopen(fname3,'r',IEEE);
fid4 = fopen(fname4,'r',IEEE);
fid5 = fopen(fname5,'r',IEEE);
fid6 = fopen(fname6,'r',IEEE);

fname5

%% Read variables: u_iso, v_iso, sig, thknss
uiso = [];
viso = [];
thknss = [];
sig = [];
sal = [];
temp = [];

% extract layer thickness, u,v in space and time
fprintf('\nReading HYCOM output: \n') 
for i=1:nt
    
    fprintf('%d/%d\n',i,nt)	
    
    for k=1:nz
        alldata1 = fread(fid1,lenrec2,'single');
        alldata2 = fread(fid2,lenrec2,'single');
        alldata3 = fread(fid3,lenrec2,'single');
        alldata4 = fread(fid4,lenrec2,'single');
        alldata5 = fread(fid5,lenrec2,'single');
        alldata6 = fread(fid6,lenrec2,'single');
        
        uiso(:,:,k,i)   = permute(reshape(alldata1(2:end-1),[nxb nyb]),[2 1]);
        viso(:,:,k,i)   = permute(reshape(alldata2(2:end-1),[nxb nyb]),[2 1]);
        thknss(:,:,k,i) = permute(reshape(alldata3(2:end-1),[nxb nyb]),[2 1]);
        sig(:,:,k,i)    = permute(reshape(alldata4(2:end-1),[nxb nyb]),[2 1]);
        temp(:,:,k,i)    = permute(reshape(alldata5(2:end-1),[nxb nyb]),[2 1]);
        sal(:,:,k,i)    = permute(reshape(alldata6(2:end-1),[nxb nyb]),[2 1]);
        
    end
end
    
fprintf('\nDone reading variables!\n')
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf]; 
a = [nbf+1:ny+nbf]; 

%% Save output to hycom (structure) 
hycom.time = t;               % time (in datenum format)
hycom.lon = lon(a,b);         % longitude 
hycom.lat = lat(a,b);         % latitude 
hycom.h   = depth(a,b);       % depth 
hycom.dz  = thknss(a,b,:,:);  % layer thickness
hycom.uiso = uiso(a,b,:,:);   % baroclinic velocity (u) 
hycom.viso = viso(a,b,:,:);   % baroclinic velocity (v) 
hycom.rho = sig(a,b,:,:);     % density  
hycom.salt = sal(a,b,:,:);    % salinity 
hycom.temp = temp(a,b,:,:);   % temperature 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
