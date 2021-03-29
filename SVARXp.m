function [A,Ax,ep,e,irf] = SVARXp(z,ext,p,h)

% this code computes a varx model with 2 lags of external variables 
% Before you run this program, you need to first run a reduced-form var for 
% external varibales to get rid of the serial correlations in ext.
% z: endogenous variables
% ext: external (exogenous) variables
% p: order of svar
% h: number of horizons to compute IRFs
% A: coefficient corresponding to z
% Ax: coefficient corresponding to ext
% ep: residuals
% e: structral shocks
% irf: impulse responses with h horizons

r=size(z,2);

y = z(p+1:end,:); 

x = [];
for i=1:p
   x = [x z(p+1-i:end-i,:)];
end

% add the external shocks with 2 lags 
x=[x ext(3:end,:) ext(2:end-1,:)];

Ap = y'*x*inv(x'*x);
ep = y - x*Ap'; 

[T,r] = size(z);
Sm = (ep-mean(ep))'*(ep-mean(ep))/T;

Ax=Ap(:,r*p+1:end); % coeffcient corresponding to the external variable
A=Ap(:,1:r*p);      % coeffcient corresponding to endogenous variables

% Recursive shock identification 
irf0 = chol(Sm,'lower');    % cholesky decomposition, initial impact matrix
e = ep*inv(irf0');          % identified structural shocks

% Impulse response analysis

% extend A for the extended VAR(1) model
A_big = [A; [eye((p-1)*r), zeros((p-1)*r,r)]]; 

% extend irf0 for the extended VAR(1) model
irf0_big = [irf0; zeros((p-1)*r,r)];    

irf_big = zeros(p*r,r,h+1);
for k=1:r
    for ii=1:h+1
        i=ii-1;
        irf_big(:,k,ii) = A_big^i*irf0_big(:,k);
    end
end

irf = irf_big(1:r,:,:);
