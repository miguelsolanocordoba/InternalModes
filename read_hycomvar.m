function hycom = read_hycom(runum,blki,blkj,var)
%%READ_HYCOMVAR reads in a variable from HYCOM's output (.BinF)
% HYCOM = READ_HYCOMVAR reads HYCOM binaries (.BinF) and saves output 
% into a Matlab structure:
%
% hycom.time  % time (in datenum format)
% hycom.lon   % longitude 
% hycom.lat   % latitude 
% hycom.h     % depth 
% hycom.dz    % layer thickness
%
% VARIABLES
% hycom.uiso  % baroclinic velocity (u) 
% hycom.viso  % baroclinic velocity (v) 
% hycom.rho   % density  
% hycom.temp  % temperature 
% hycom.salt  % salinity
% hycom.ubar  % eastward barotropic velocity
% hycom.vbar  % northward barotropic velocity
% hycom.ssh   % sea surface height 
% 
% Created: October 3, 2020 by M. Solano 

% Format 
IEEE = 'ieee-be';
addpath /data/msolano/Matlab

%% Experiment and tile number 
% North Atlantic > runnum=221;  blki=27; blkj=45;
% South Pacific > runnum=190;  blki=15; blkj=25;
% Amazon (1) > runnum=190;  blki=19; blkj=40;
% Amazon (2) > runnum=190;  blki=19; blkj=41;
% Amazon (3) > runnum=190;  blki=18; blkj=40;
% Amazon (4) > runnum=190;  blki=18; blkj=41;

%runnum = 190; 
dirin = '/data2/msolano/hycom/GLBc0.04/expt_19.0/'; % EXPT_19.0
%runnum =221; dirin = '/data2/mbui/for_keshav/tiles/';             % EXPT_22.1
runnumstr = num2str(runnum);
%blki=19;
%blkj=40;

% Print some info
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
fnamedz = [dirin 'thknss/thknss_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];

% Variables 
switch var
   case 'uiso'
fname = [dirin 'u_iso/u_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 3; 
   case 'viso'
fname = [dirin 'v_iso/v_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 3; 
   case 'rho'
fname = [dirin 'sig/sig_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 3; 
   case 'temp'
fname = [dirin 'temp/T_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 3; 
   case 'salt'
fname = [dirin 'sal/S_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 3; 
   case 'ssh'
fname = [dirin 'srfhgt/srfhgt_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 2; 
   case 'steric'
fname = [dirin 'steric/steric_' num2str(runnum) '_blk_' ...
           num2str(blki) '_' num2str(blkj) '.BinF'];
      vardim = 2; 
   otherwise 
      disp('VAR must be: uiso, viso, rho, temp, salt, ssh, steric')
      return
end


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
fiddz = fopen(fnamedz,'r',IEEE);
fid = fopen(fname,'r',IEEE);

%% Read selected variable
vars = [];
thknss = [];

% extract layer thickness, u,v in space and time
fprintf('\nReading HYCOM output: \n') 
for i=1:nt
    fprintf('%d/%d\n',i,nt)	
    for k=1:nz
        alldatadz = fread(fiddz,lenrec2,'single');
        thknss(:,:,k,i) = permute(reshape(alldatadz(2:end-1),[nxb nyb]),[2 1]);
    end
end
   
if vardim==3
   for i=1:nt
       for k=1:nz
           alldata = fread(fid,lenrec2,'single');
           vars(:,:,k,i)   = permute(reshape(alldata(2:end-1),[nxb nyb]),[2 1]);
       end
   end
elseif vardim==2
   for i=1:nt
       alldata = fread(fid,lenrec2,'single');
       vars(:,:,i)   = permute(reshape(alldata(2:end-1),[nxb nyb]),[2 1]);
   end
end


fprintf('\nDone reading variables!\n')
fclose(fid);
fclose(fiddz);

% Don't save halos (nbf) 
b = [nbf+1:nx+nbf]; 
a = [nbf+1:ny+nbf]; 

%% Save output to hycom (structure) 
hycom.time = t;               % time (in datenum format)
hycom.lon = lon(a,b);         % longitude 
hycom.lat = lat(a,b);         % latitude 
hycom.h   = depth(a,b);       % depth 
hycom.dz  = thknss(a,b,:,:);  % layer thickness
switch vardim
   case 2  % 2D variable
      hycom.vars = vars(a,b,:);    
   case 3  % 3D variable
      hycom.vars = vars(a,b,:,:);  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EoF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
