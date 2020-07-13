function p = comp_pertpress(rho,dz) 
%%COMP_PERTPRESS computes perturbation pressure from density
% 
%
% Created: June 23, 2020 by M. Solano

g = 9.8; % gravity  

% Dimensions
[nz,nt] = size(rho); 
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


