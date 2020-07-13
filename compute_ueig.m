function Ueig = compute_ueig(Weig,dz)
%%COMPUTE_UEIG(Weig,dz) computes Ueig=d(Weig)/dz (cell centers)
% 
%
% Created: May 31, 2020

H = nansum(dz);  % Depth 
N = numel(dz);   % Layers (faces) 

% Compute the derivative
dW2= (Weig(2:end,:)- Weig(1:end-1,:));
dzu = repmat(dz,[1,size(dz,1)-1]);
Ueig1 = dW2./dz;

% Normalize (see Buijsman et al. 2020)
AA=repmat(sum(Ueig1.^2.*dz,1)./H,[N 1]).^(1/2);
AA(AA==0)=Inf;
Ueig2=Ueig1./AA;
Ueig2(:,Ueig2(N,:)<0)=-Ueig2(:,Ueig2(N,:)<0); %% bottom

% Normalized horizontal eigenvalue at cell center 
Ueig = Ueig2; % OUTPUT
