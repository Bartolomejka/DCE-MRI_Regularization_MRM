function [Cfit,dCfit,aif_and_h]=C_fit_FT(x,t,model,rescale_parameters,x_fixed)
% [Cfit,dCfit,aif_and_h]=C_fit_FT(x,t,model,rescale_parameters,x_fixed)
% renders perfusion curve/curves Cfit = F^-1{AIF(x_aif,w)*IRF(x_irf,w)}
% and its/theirs derivatives: dCfit/dx = [dCfit/dx_aif; dCfit/dx_irf]
% together with AIF and IFR in time domain (costs 2 more FFTs), if needed
% Inputs: x - model parameters, vector or mtx if more curves: [x_aif;x_irf]
%       t - time points
%       model - structure containing model information:
%           .aif: .model - function handle to AIF-model, i.g.: @aif_koh_norm_FT
%                 .no_of_par - number of parameters expected
%           .irf: .model - function handle to IRF-model, i.g.: @TH_trunc_FT
%                 .no_of_par - number of parameters expected
%       rescale_parameters - constant used to rescale parameters p
%       x_fixed - logical vector of parameters not to be optimized

% separately evaluate AIF and IRF

no_of_timesamples_or=length(t);

% no_of_curves=length(y)/no_of_timesamples;
no_of_curves=size(x,1);

no_of_aif_par=model.aif.no_of_par;
no_of_irf_par=model.irf.no_of_par;

[x_aif x_irf]=split_x_single_channel(x,rescale_parameters,no_of_aif_par);

aif_rescale_parameters=rescale_parameters(and([true(1,no_of_aif_par) false(1,no_of_irf_par)],~x_fixed));
irf_rescale_parameters=rescale_parameters(and([false(1,no_of_aif_par) true(1,no_of_irf_par)],~x_fixed));

% time-rescale
t_rescale=1;
Ts=t(2)-t(1);
t=(t(1):Ts/t_rescale:t(end))';
% t=t(1:end-~mod(length(t),2)); % why? see fourier2time
% Ns=length(t);

%compute anti-aliasing limits
if strcmp(func2str(model.irf.model),'TH_FT')
    E=1-exp(-x_irf(:,4)./x_irf(:,1));
    Kep=x_irf(:,1).*E(:)./x_irf(:,3);
    % Ns=ceil(max(length(t) + 6*1./Kep/Ts));
    %     Ns=ceil(max(length(t) + 6*1./Kep/Ts+abs(x_irf(:,end))/Ts));
    Ns=length(t) + 10*1./Kep/Ts+abs(x_irf(:,end))/Ts; % ja davam 10 kvuli derivacim
elseif strcmp(func2str(model.irf.model),'TH_Sourbron')
    %     Ns=ceil(max(length(t) + 10*x_irf(:,4)/Ts+abs(x_irf(:,end)+x_irf(:,3))/Ts));
    %     Ns=ceil(max(length(t) + 10*x_irf(:,2)/Ts+abs(x_irf(:,end))/Ts)); % experimentalni jen z celkoveho T
    TTc=(x_irf(:,2)-x_irf(:,3));
    Ns=length(t) + 10*TTc./(1-exp(-TTc./x_irf(:,4)))/Ts+abs(x_irf(:,end))/Ts; % presne % ja davam 10 kvuli derivacim (PS)
else
    Ns=length(t)*3*ones(no_of_curves,1); % 3 - because of derivatives
    %     Ns=ceil(max(length(t) + 50*x_irf(:,2)/Ts+abs(x_irf(:,end))/Ts));
end
Ns=max(ceil(Ns),length(t));
no_of_timesamples=length(t);
Ns=2.^nextpow2(Ns); % nove 2017
% Ns=repmat(max(Ns),no_of_curves,1); % if speed up not used
[Ns_unique, m, ind] = unique(Ns);
% Ns_unique=Ns;ind=1; % IUPESM2018 comparison
Cfit=zeros(no_of_timesamples,no_of_curves,class(t));
if nargout>1
    dCfit=zeros(no_of_timesamples_or,nnz(~x_fixed),no_of_curves,class(t));
end
pom_ind=1:no_of_curves;
for n_Ns=1:length(Ns_unique)
    selected_Ns=ind==n_Ns;
    Ns=Ns_unique(n_Ns);
    Ns_limit=32768; %2^17
    if Ns>(Ns_limit+1) 
        warning(['Number of samples in Fourier domain is high:' num2str(Ns) '. Possible aliasing.'])
        Ns=Ns_limit; %2^17
    end
    no_of_curves=nnz(selected_Ns);
    ind_vector=pom_ind(selected_Ns);
    % Ns=Ns+~mod(Ns,2);
    % w=(0:floor(Ns/2))*2*pi/Ts/(Ns-1); % stara
    w=(0:(Ns/2-mod(Ns,2)/2))*2*pi/Ts/Ns; % nova 2017
    w=w(:);

    
    t_cut=t(end)-x_aif(selected_Ns,end); % kvuli aif_nonblind urychleni - neprovereno !!! asi neni dobre pro zaporne posuny v IRF
%     t_cut=t(end)-x_aif(selected_Ns,end)-x_irf(selected_Ns,end);
    % separete AIF and TRF - universal version
    
    % AIF
    if nargout<2
        % without derivatives
        aif=zeros(length(w),no_of_curves);
        h=zeros(length(w),no_of_curves);
        aif_ind=~x_fixed(1:no_of_aif_par);
        for n=1:no_of_curves
            aif(:,n)=model.aif.model(x_aif(ind_vector(n),:),w,t_cut(n),aif_ind);
            h(:,n)=model.irf.model(x_irf(ind_vector(n),:),w);
            Cfit(:,ind_vector(n))=fourier2time(h(:,n).*aif(:,n),no_of_timesamples);
        end
    else
        % with derivatives
        aif=zeros(length(w),no_of_curves);
        h=zeros(length(w),no_of_curves);
        % IRF
        aif_ind=~x_fixed(1:no_of_aif_par);
        irf_ind=~x_fixed((no_of_aif_par+1):end);
        for n=1:no_of_curves
            [aif(:,n) daif]=model.aif.model(x_aif(ind_vector(n),:),w,t_cut(n),aif_ind);
            [h(:,n) dh]=model.irf.model(x_irf(ind_vector(n),:),w,irf_ind);
            Cfit(:,ind_vector(n))=fourier2time(h(:,n).*aif(:,n),no_of_timesamples);
            % derivatives
            %         if any(x_fixed)
            %             daif(:,x_fixed(1:no_of_aif_par))=[]; % erase fixed aif-parameters
            %             dh(:,x_fixed((no_of_aif_par+1):end))=[]; % erase fixed parameter derivatives; 6 is insted of number of parameters of the irf-model
            %         end
            
            for m=1:size(daif,2) %dCfit/dx_AIF
                dCfit_daif=fourier2time(h(:,n).*daif(:,m),no_of_timesamples);
                dCfit_daif=dCfit_daif./aif_rescale_parameters(m);
                dCfit(:,m,ind_vector(n))=dCfit_daif;
            end
            
            M=length(aif_rescale_parameters);
            for m=1:size(dh,2) %dCfit/dx_IRF
                dCfit_dirf=fourier2time(dh(:,m).*aif(:,n),no_of_timesamples);
                dCfit_dirf=dCfit_dirf./irf_rescale_parameters(m);
                dCfit(:,m+M,ind_vector(n))=dCfit_dirf;
            end
            
        end        
        dCfit=dCfit*t_rescale;
    end
    
    if nargout>2
        for n=1:no_of_curves
            aif_and_h.aif(:,ind_vector(n))=fourier2time(aif(:,n),no_of_timesamples);
            aif_and_h.irf(:,ind_vector(n))=fourier2time(h(:,n),no_of_timesamples);
        end
    end
    
end
