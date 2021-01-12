% This code computes the historical decomposition of aggregate variables 
% and functional factors
% based on notes, A Primer On Vector Autoregressions by Ambrogio Cesa-Bianchi
% The results contains Ka+2 columns
% The first column gives the initial values values of the HD
% The rest of the clolumns give the HD due to the Ka Shocks.
% This version of HD has initial condtions, which is different from Kalian
% textbook one (w/ initial condtions)

%% Historical Decompostion
function [hisdec]= hist_decomp(B,e,A,Z,T)
%    B: on impact impulse response
%    e: structral shocks
%    A: reduced-form parameters
%    Z: orginal time series 
%    T = 168

global r p Ka Kf

X = [];
for i=1:p
   X = [X Z(p+1-i:end-i,:)];
end
% extend A from VAR(p) for the extended VAR(1) model
A_big = [A; [eye((p-1)*r), zeros((p-1)*r,r)]]; 

%% compute the IRF from period 0 to period T-p, which is needed to get HD
irf0_big = [B; zeros((p-1)*r,Ka+1)];    
irfhis_big = zeros(p*r,Ka+1,T-p);
for k=1:Ka+1
    for ii=1:T-p
        i=ii-1;
        irfhis_big(:,k,ii) = A_big^i*irf0_big(:,k);
    end
end
irfhis = irfhis_big(1:r,:,:);

%% HD from y_p+1
hisdec=zeros(Ka+Kf,Ka+1,T-p); % Ka+2:one more clolumn for initial condition

% e is T-p dimension, with first elemenet at time p+1
% HD ts starting from period p+1, which is T-p horizons
% First calculate the HD for the (p+1)th observation, which does not 
% innner product.
for j=1:r  %variables
    for k=1:Ka+1  %shocks
    hisdec(j,k,1)=squeeze(irfhis(j,k,1))*e(1,k);
    end
end

%%Next calculate the HD for the rest of the obsevations, which require
%%inner products.
for j=1:r %variables
    for k=1:Ka+1  %shocks
        for i=p:T-p %horizons
    hisdec(j,k,i)=dot(squeeze(irfhis(j,k,1:i)),e(i:-1:1,k));
        end
    end
end

hisdec_int=zeros(r,1,T-p);

%% initial conditions for all 
for i=p:T-1 % horizons
        Temp=X(1,:)*((A_big')^(i-(p-1)));
        hisdec_int(:,1,i-p+1)=Temp(1:Ka+Kf);
end

hisdec=cat(2,hisdec_int,hisdec);

