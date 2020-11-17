function [N2mod,rhomod] = N2maskf(rho,rho0,N2,dz)
%%N2MASK masks negative N2 to 1-10
%
% N2 at centers, rho at faces, dz = diff(zf)
%
% Created: Sept. 2, 2020 by M. Solano

g = 9.806;
nz = numel(N2); % number of layers 
N2mod = N2; 
N2mod(N2<0)=1e-10;

rhomod = zeros(size(rho));
rhomod(nz+1) = rho(nz+1);
for i = nz:-1:1
    rhomod(i) = rhomod(i+1) + rho0/g*dz(i)*N2mod(i);
end

%rhodiff = rhomod(3) - rhomod(2); 
%drhodz = rhodiff/dz(2);
%rhomod(1) = rhomod(2) + drhodz*dz(1);
