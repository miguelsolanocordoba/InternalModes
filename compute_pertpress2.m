function p = compute_pertpress2(rho,zc,zf) 
%%COMPUTE_PERTPRESS2 computes perturbation pressure from eta-method
%  
%
% Created: September 15, 2020 by M. Solano


g = 9.806; % gravity  

% Dimensions
[nz,nt] = size(rho); 

% Compute mean layer
zcmean = mean(zc,2); 
zfmean = mean(zf,2); 

% Compute mean density in time-mean layer
dz = diff(zfmean);   % time-mean layer thickness
rhoint = zeros(size(rho)); 
rhof = zeros(nz+1,nt); 
for i=1:nt
    rhoint(:,i) = interp1(zc(:,i),rho(:,i),zcmean);
%    rhof(:,i) = interp1(zc(:,i),rho(:,i),zf(:,i),'linear','extrap'); 
    rhof(:,i) = interp1(zc(:,i),rho(:,i),zfmean,'linear','extrap'); 
end
rho_mean = mean(rho,2); 
%rho_mean = mean(rhoint,2); 
rhof_mean = mean(rhof,2); 

% Interpolate mean density to faces and compute N2
%rhof_mean = interp1(zcmean,rho_mean,zfmean,'nearest','extrap'); 
N2 = compute_bvf(rhof_mean,1035.4312,zfmean); 

% Pre-allocate 
rhodiff = zeros(nz,nt); 
prho = zeros(nz,nt); 

% Compute density anomaly 
for i=1:nt
    eta = zc(:,i) - zcmean;
    rhodiff(:,i) = (rho_mean/g).*N2.*eta + rho(:,i) - rho_mean; 
    prho(:,i) = g*rhodiff(:,i).*dz; 
end

prho1 = zeros(size(prho)); 

for i=1:nz
    for j=1:nt
	prho1(i,j) = nansum(prho(1:i,j)); 
    end
end

% Depth-averaged pressure
hc = nansum(dz); 
prhodz = squeeze(sum(prho1.*repmat(dz,[1,nt]),1))./squeeze((repmat(hc,[1,nt])));

% Perturbation pressure
for i = 1:nz
    p(i,:) = squeeze(prho1(i,:)) - prhodz; 
end
