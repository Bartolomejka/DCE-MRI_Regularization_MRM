function [x_aif x_irf]=split_x_single_channel(x,rescale_parameters,no_of_aif_par)
% [x_aif x_irf]=split_x(x,x_fixed,x_start,no_of_aif_par,no_of_curves)
% Split vector x into subvecors x_aif (aif-model parameters) and
% x_irf (irf-model parameters)
% x has dimensions rows: no_of_curves, columns: (no_of_aif_par+no_of_irf_par)

x=x./repmat(rescale_parameters,size(x,1),1);  % rescale parameters

% separete AIF and TRF
x_aif=x(:,1:no_of_aif_par);
x_irf=x(:,(no_of_aif_par+1):end);

% % [x_aif x_irf]=split_x(x,x_fixed,x_start,no_of_aif_par,no_of_curves)
% % Split vector x into subvecors x_aif (aif-model parameters) and
% % x_irf (irf-model parameters)
% % x has dimensions rows: (no_of_aif_par+no_of_irf_par), columns: no_of_curves
% 
% x=x./repmat(rescale_parameters,1,size(x,2));  % rescale parameters
% 
% % separete AIF and TRF
% x_aif=x(1:no_of_aif_par,:);
% x_irf=x((no_of_aif_par+1):end,:);
end