function [y_minus_Cfit_L2_2 dy_minus_Cfit_L2_2 H_est]= y_minus_Cfit_L2_2...
    (p,t,y,model,rescale_parameters,x_fixed)
% [y_minus_Cfit_L2_2 dy_minus_Cfit_L2_2 H_est]= y_minus_Cfit_L2_2...
%   (p,t,y,model,rescale_parameters,x_fixed)
% computes sum of squared differences || y - Cfit(p,t) ||_2^2 and
% its derivatives and Hessian approximation (J'.*J)
% Inputs: p - model parameters, vector or mtx if more curves
%       t - time points
%       y - measured data, vector or mtx if more curves
%       model - structure containing model information, see C_fit_FT.m
%       rescale_parameters - constant used to rescale parameters p
%       x_fixed - logical vector of parameters not to be optimized
% Outputs:
%       y_minus_Cfit_L2_2 ...sum of squared differences: || y - Cfit(p,t) ||_2^2
%       dy_minus_Cfit_L2_2 ...gradient: -2*dCfit/dp.'*(y - Cfit(p,t))
%       H_est ... Hessian approximation: 2*dCfit/dp.'*dCfit/dp

% memory demanding version
[Cfit,dCfit]=C_fit_FT(p,t,model,rescale_parameters,x_fixed);

no_of_derivatives=size(dCfit,2);
no_of_curves=size(Cfit,2);
N=length(t);

% stochastic sampling
% sampling_pct=0.2;
% sampling=ceil(N*rand(no_of_curves,ceil(N*sampling_pct)));
sampling=true(no_of_curves,N);
sampling_pct=1;


y_minus_Cfit=y-Cfit;
dy_minus_Cfit_L2_2=zeros(no_of_curves,length(x_fixed),class(t));
for n=1:no_of_curves
    dy_minus_Cfit_L2_2(n,~x_fixed)=-2*sum(repmat(...
        y_minus_Cfit(sampling(n,:),n),1,no_of_derivatives).*...
        dCfit(sampling(n,:),:,n),1)./sampling_pct;
end
y_minus_Cfit_L2_2=y_minus_Cfit(:).'*y_minus_Cfit(:);

if nargout>2
H_est=cell(no_of_curves,1);
    for n=1:no_of_curves
        H_est{n}=2*dCfit(sampling(n,:),:,n).'*dCfit(sampling(n,:),:,n);
    end
end