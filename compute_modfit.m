function [r2o,fit] = compute_modfit(var,dz,Ueig,nmodes) 
%%COMPUTE_MODFIT Computes fits for velocities and pressure perturbations to each vertical mode
% dz represents the layer thickness vector
% var represents either velocity (baroclinic) or pressure perturbation
% Ueig represents hor. eig. function array
% Created: May 31, 2020 by M. Solano

% Dimensions
[nz,nt] = size(var);

% Depth 
H = nansum(dz); 

% Pre-allocate
var_fit  = zeros(nz,nt,nmodes);

% Main loop
rem = var; 
r2 = zeros(nt,nmodes); 
r2o = zeros(nmodes,1); 
SS_tot = zeros(nt,nmodes); 
SS_res = zeros(nt,nmodes); 

for n = 1:nmodes
    amp = sum(repmat(Ueig(:,n),[1,nt]).*rem.*repmat(dz,[1,nt]),1)./H; 
    var2 = repmat(amp,[nz,1]).*repmat(Ueig(:,n),[1,nt]); 
    var_fit(:,:,n) = var2; 
    rem = rem-var2;
    
%    for i = 1:nt
%	ymean = sum(var(:,i).*dz)./H;
%	SS_tot(i,n) = sum((var(:,i)).^2.*dz);
%	SS_res(i,n) = sum(((rem(:,i))).^2.*dz); 
%	SS_fit(i,n) = sum((var2(:,i)).^2.*dz); 
%
%	r2(i,n) = SS_res(i,n)./SS_tot(i,n); 
%    end

    SS_tot = sum(var.^2.*dz,1); 
    SS_res = sum(rem.^2.*dz,1); 
    r2o(n) = 1 - mean(SS_res./SS_tot); 
%    r2o(n) = 1-mean(r2(:,n)); 
end
fit = var_fit; 
