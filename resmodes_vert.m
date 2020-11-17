function num_modes = resmodes_vert(rho,dz) 
%% RESMODES_VERT computes number of resolvable modes
%
%
% Created: October 19, 2020 by M. Solano


% Compute number of resolvable modes (nres)
Weig1 = zeros(nz,10);
Weig1 = squeeze(Weig(i,j,1:nz,1:10));
nres = 1; % Start checking from mode 2
resflag = true;
nzH = numel(zf_mean);
A = zeros(nz-1,1);
while resflag
   nres = nres + 1;
   A = Weig1(1:end-1,nres+1).*Weig1(2:end,nres+1);
   ind = find(A<0);
   wzeros = 0.5*(zfM(ind)+zfM(ind+1));
   wzeros = [-H;wzeros;0]; % add end-points (surface/bottom)

   for ii = 1:numel(wzeros)-1
      if sum(zf_mean>wzeros(ii) & zf_mean<wzeros(ii+1))==0;
           resflag = false;
           num_modes(i,j) = nres;
      end
   end

   if nres==9
      num_modes(i,j) = 10;
      break
   end

end %while
