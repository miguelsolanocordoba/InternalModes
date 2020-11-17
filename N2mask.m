function [N2mod,rhomod] = N2mask(rho,rho0,N2,dz)
%%N2MASK masks negative N2 to 1-10
%
% N2 at faces, rho at centers, dz = diff(zc);
%
% Created: Sept. 2, 2020 by M. Solano

g = 9.806; % gravity

% Number of layers (rho must be at centers)
nz = numel(rho);
N2mod = N2; 
N2mod(N2<0)=1e-10; % mask negative N2 to 1e-10
%N2mod = abs(N2);

rhomod = zeros(size(rho));
rhomod(nz) = rho(nz); % surface sea-water density: rho0 = rho(end)
for i = nz-1:-1:1
    rhomod(i) = rhomod(i+1) + rho0/g*dz(i)*N2mod(i+1);
end
