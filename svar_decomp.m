function vardec=fsvar_decomp(irf_m,Ka,h)
%% compute the variance decomposition with the median response
% irf_m: median of irfb (bootstrapped impluse responses of the aggregate variables and functional factors)
% Ka: number of aggregate varibles
% h: horizons
% variance decomposition for Ka aggregate variables
% at horizon h computed using the median reponse
% The rth row and kth column of V(:,:,h) is the contribution of kth shock to 
% the rth aggregate variable
% The rth row and kth column of V_total(:,h) is the total forecast variance of
% rth aggregate variable
% vardec: Ka*Ka*h
% The rth row and the kth column of vardec(:,:,h) is the percentage of contribution
% kth shock to the forcast variance of the rth aggregate variable to the kth shock.

V=zeros(Ka,Ka,h);% ith row and jth column stores the variance of ith varable due to jth shock
V_total=zeros(Ka,h);% total variance of aggregate variables

V(:,:,1)=irf_m(:,:,1).^2;
V_temp(:,:,1)=irf_m(:,:,1)*irf_m(:,:,1)';
V_total(:,1)=diag(V_temp(:,:,1));

for i=2:h
	V(:,:,i)=V(:,:,i-1)+irf_m(:,:,i).^2;
	V_temp(:,:,i)=V_temp(:,:,i-1)+irf_m(:,:,i)*irf_m(:,:,i)';
	V_total(:,i)=diag(V_temp(:,:,i));
end

vardec=zeros(Ka+Kf,Ka+1,h);
% calculate the percentage variance of each variable due to the kth shock
for i=1:h
for k=1:Ka
	vardec(:,k,i)=V(:,k,i)./V_total(:,i);
end
end


	



