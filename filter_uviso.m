function [uiso,viso] = filter_uviso(u_iso,v_iso)
%%FILTER_UVISO filters velocities
%
%
% Created: June 1, 2020 by M. Solano

% Parameters
fc1 = 1/(15/24);      % lower cutoff frequency (cpd)
fc2 = 1/(9/24);     % upper cutoff frequency (cpd)
dt = 1/24;    % time interval (days)
N = 3;     % n-nd order of butterworth filter (1-9)
frc = 0.4;    % fraction of data to be added to the ends of original time series (padding)

% uiso
var2 = u_iso;

[ll,mm] = size(var2);
uiso = zeros(ll,mm);

for i = 1:ll
    var3 = squeeze(var2(i,:));
    var4 = bandpass_mpad(var3,fc1,fc2,dt,N,1,frc);
    uiso(i,:) = var4;        
end

% viso
var2 = v_iso;
[ll,mm] = size(var2);
viso = zeros(ll,mm);

for i = 1:ll
    var3 = squeeze(var2(i,:));
    var4 = bandpass_mpad(var3,fc1,fc2,dt,N,1,frc);
    viso(i,:) = var4;        
end
