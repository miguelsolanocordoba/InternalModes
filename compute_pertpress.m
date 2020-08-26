function p = compute_pertpress(rho1,zc,zcmean,zfmean) 
%%COMPUTE_PERTPRESS computes perturbation pressure from density
%  
%
% Created: June 23, 2020 by M. Solano


g = 9.806; % gravity  

% Dimensions
[nz,nt] = size(rho1); 

% Compute mean density in time-mean layer
%z = mean(zc,2); % time-mean layer
dz = diff(zfmean);   % time-mean layer thickness
rho = zeros(size(rho1)); 
for i=1:nt
    rho(:,i) = interp1(zc(:,i),rho1(:,i),zcmean);
end
rho_mean = mean(rho,2); 

% Pre-allocate 
rhodiff = zeros(nz,nt); 
prho = zeros(nz,nt); 

% Compute density anomaly 
for i=1:nt
    rhodiff(:,i) = rho(:,i) - rho_mean; 
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


