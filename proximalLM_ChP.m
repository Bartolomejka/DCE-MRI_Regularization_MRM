%% minimisation: proximal Newton/gradient method with ChP-denoising substep - corrected
TV_regularization=true;%false;%
N_proximal=50; % no of main proximal-gradient iterations
gamma_general=10^-(2/3); % SUBMITTED PAPER - set in outer script

lambda_base=1e-3; % base step-length
lambda_min=eps*1e6;
lambda_max=lambda_base*1e6;

% profile on
time_spent=0;
tic
clear dy_minus_Cfit_L2_2 G x_iter x_iter_before SSD lambda L1_func_iter d func_eval rho



if strcmp('TH_Sourbron2_FT',func2str(model.irf.model))
    gamma_no_scale=gamma_general*[0 0.0248    0.2826    0.0243    0.1031    0.5652 0]; % newly derived from ML with normalization with brain, arteries removal; - SUBMITTED PAPER
else
    gamma_no_scale=1e-1*[0  0.1 1 1 0.1 1]; % TH Garpebring - sampled aif (nove mene v Fp)
end

gamma=gamma_no_scale./rescale_parameters*length(t); % normalization and scaling

% prepare starting points
x_iter(1:no_of_curves,:,1)=repmat(x_start_rescaled,[no_of_curves/size(x_start_rescaled,1) 1 1]);

% initialization
ni=2*ones(no_of_curves,1);
rho=ones(no_of_curves,1);
SSD=zeros(no_of_curves,1);
func_eval=ones(no_of_curves,1);
no_of_par=nnz(~x_fixed);
ind_to_par=1:length(x_fixed);
ind_to_par=ind_to_par(~x_fixed);
mask_im=form_image((1:no_of_curves).',tissue(:,2));
mask_ind=mask_im(phantom_mask_crop);
c_ind=1:no_of_curves;
c_ind=c_ind.';
curve_step=1;
ub_r=ub;
lb_r=lb;

% reconstruct starting image and compute prior term
P_est_iter=form_image(x_iter(:,~x_fixed,1),data_separate(:,2));
% compute initial regularization term values - debug
% isotropic
for p=1:no_of_par
    L1_image=sqrt(sum(fgrad_1( P_est_iter(:,:,p),phantom_mask_crop).^2,3)); % compute gradient image
    L1_func(p)=gamma(ind_to_par(p))*sum(noise_std_map(phantom_mask_crop).^2.*L1_image(phantom_mask_crop));
    L1_func(p+no_of_par)=L1_func(p); % for no TV case
end
L1_func_iter(1,:)=L1_func;
lambda=lambda_base*ones(no_of_curves,1);
H=cell(no_of_curves,1);
% main proximal gradient iterations
for it=1:N_proximal%10
    iter=size(x_iter,3)
    x_iter(:,:,iter+1)=x_iter(:,:,iter); % if separate curves
    SSD(:,iter+1)=SSD(:,iter);
    if ~TV_regularization
        % exclude curve from fitting if allready in minimum
        c_ind(c_ind==0)=[];
        disp(['Curves in progress: ' num2str(length(c_ind))])
    end
    if isempty(c_ind) % stop if no curve rest
        break
    end
    for c_it=1:curve_step:length(c_ind) % go trough curves if set up
        c=c_ind(c_it,:);
        % correct ChP denoising if result is out of domain
        x_iter(c,~x_fixed,iter)=max(x_iter(c,~x_fixed,iter),lb_r);
        x_iter(c,~x_fixed,iter)=min(x_iter(c,~x_fixed,iter),ub_r); % correct ChP denoising if result is out of domain
        
        % set up criterial function for current curve
        crit_func=@(x)y_minus_Cfit_L2_2(x,t,[data_separate{c,1}],model,rescale_parameters,x_fixed);
        
        % evaluate crit.fcn for particular method
        [SSD_current dy_minus_Cfit_L2_2 H_est]=crit_func(x_iter(c,:,iter));
        H_est{1}=H_est{1}*noise_std(c).^-2; % for LM: I+J'WJ, i.e. noise inside
        
        dy_minus_Cfit_L2_2=dy_minus_Cfit_L2_2*noise_std(c).^-2;
        SSD_current=SSD_current*noise_std(c).^-2;
        
        n=iter;
        lambda(c,n)=lambda(c,max(n-1,1)); % just for if on line 197 in case of LM
        ni(c,n)=ni(c,max(n-1,1));
        y=x_iter(c,:,n);
        
        SSD(c,n)=SSD_current;
        func_eval(c,n)=1;
        
        bound_check=false(2*curve_step,no_of_par);
        
        % 1) perform step downwards
        if (rho(c,n)<0) % if wrong approximation - reduce lambda
            lambda(c,n)=lambda(c,max(n-1,1)); % presunuto z 137 kvuli zaporne lambde v BB
            lambda(c,n)=lambda(c,n)/ni(c,n);
            ni(c,n)=ni(c,n)*2; % this LM-modification [Madsen] is not necessary
            lambda(c,n)=max(lambda(c,n),lambda_min);
        else % if good approx. - keep/rise lambda
            lambda(c,n)=lambda(c,max(n-1,1))/max(1/3,1-(2*rho(c,n)-1)^3); % Levenberg-Marquardt scheme [Madsen]
            ni(c,n)=ni(c,1);
        end
        % truncate lambda to limits
        lambda(c,n)=min(lambda(c,n),lambda_max);
        lambda(c,n)=max(lambda(c,n),lambda_min);
        
        is_minimizing=false;
        
        lambda(c,n)=lambda(c,n)*2; % *2 due to first run of while
        while ~is_minimizing
            lambda(c,n)=lambda(c,n)/2;
            d=-dy_minus_Cfit_L2_2;
            % check and inverse the Hessian
            for nH=1:curve_step
                if any(isnan(H_est{nH}(:)))
                    H_est{nH}=eye(no_of_par);
                    warning('Hessian contains NaNs, gradient is used.')
                    if any(isnan(dy_minus_Cfit_L2_2(nH,~x_fixed)))
                        dy_minus_Cfit_L2_2(isnan(dy_minus_Cfit_L2_2))=0;
                        warning('Gradient contains NaNs, no step performed.')
                    end
                end
                H{(c(nH))}=eye(no_of_par)/lambda(c,n)+H_est{nH}; % for LM: I+J'WJ, i.e. noise (W) inside
                while rcond(H{(c(nH))})<1e-6
                    lambda(c,n)=lambda(c,n)/2;
                    H{(c(nH))}=eye(no_of_par)/lambda(c,n)+H_est{nH}; % for LM: I+J'WJ, i.e. noise (W) inside
                end
                d(nH,~x_fixed)=H{(c(nH))}\-dy_minus_Cfit_L2_2(nH,~x_fixed).'; % tohle je lepsi, ale potrebuji iH
            end
            
            % perform the step
            x_iter(c,:,n+1)=y+d; %% pozor pri GD
            
            % corrrect if current solution is outside region of definition
            x_iter(c,~x_fixed,n+1)=max(x_iter(c,~x_fixed,n+1),lb_r);
            x_iter(c,~x_fixed,n+1)=min(x_iter(c,~x_fixed,n+1),ub_r);
            
            % compute new SSD value
            [SSD_current grad]=crit_func(x_iter(c,:,n+1));
            grad=grad*noise_std(c).^-2;
            SSD_current=SSD_current*noise_std(c).^-2;
            func_eval(c,n)=func_eval(c,n)+1;
            
            % compute quality of approximation rho [Gavin](15)
            rho(c,n+1)=(SSD(c,n)-SSD_current)/(lambda(c,n)*d(1,~x_fixed)*(d(1,~x_fixed)-dy_minus_Cfit_L2_2(1,~x_fixed)).'); % Maruska - spravne
            
            % check if new solution is better
            if ((SSD_current-SSD(c,n))<=0)||(lambda(c,n)<=lambda_min)%<1e-8)%eps)%
                dy_minus_Cfit_L2_2=grad;
                is_minimizing=true; % breaks the while-loop
            end
        end
        % stop iterating when minimum found
        if (abs(SSD_current-SSD(c,n))<2*eps)%1e-15)%1e-12) % %  % stop iterating when minimum found
            x_iter(c,:,(n+1):(iter+1))=repmat(x_iter(c,:,(n+1)),1,iter+1-n);
            SSD(c,(n+1):(iter+1))=SSD_current;
            lambda(c,(n+1):(iter+1-1))=0;%lambda_base;%lambda{c}(n);
            if ~TV_regularization
                % exclude curve from further fitting
                c_ind(c_it)=0;
            end
            continue
        end
    end
    
    % 2) TV denoising
    % 2a) reconstruct image: parameter vector -> image
    P_est_iter=form_image(x_iter(:,~x_fixed,end),data_separate(:,2));
    P_before=P_est_iter;
    
    % compute prior-term value before ChP denoising
    % isotropic
    for p=1:no_of_par
        L1_image=sqrt(sum(fgrad_1( P_est_iter(:,:,p),phantom_mask_crop).^2,3));
        L1_func(p)=gamma(ind_to_par(p))*sum(L1_image(phantom_mask_crop));
        L1_func(p+no_of_par)=L1_func(p); % for no TV case
    end
    %     Proximal step
    if TV_regularization
        % Chambolle-Pock denoising
        [P_pom reg_term]=PD_denoising_H_TV(P_est_iter,H,phantom_mask_crop,gamma(~x_fixed));
        
        L1_func(no_of_par+1:end)=reg_term(2,:);
        P_est_iter=P_pom;
        
        % reformat estimated parameters: image->vector
        x_iter(:,ind_to_par,n+1)=form_vector(P_est_iter,tissue(:,2));
    end
    L1_func_iter(n+1,:)=L1_func;
    
    % Plot performance
    if mod(it,20)==0
        %%
        figure(1050)
        subplot(2,2,1)
        semilogy(1:iter,(sum(SSD(:,1:iter),1)))
        title('\sigma^-^2\Sigma||y-C(n,p)||_2^2')
        subplot(2,2,2)
        l_mean=mean(lambda,1);
        l_min=min(lambda,[],1);
        l_max=max(lambda,[],1);
        semilogy([l_min;l_mean;l_max].')
        title('\lambda')
        %         legend show
        figure(1061)
        which_par=no_of_aif_par+4;%ind_to_par;
        which_curves=unique(ceil(no_of_curves*rand(1,50)));
        subplot(1,2,1)
        plot(squeeze(x_iter(which_curves,which_par,1:iter)).')
        title(x_label(which_par))
        subplot(1,2,2)
        x_iter_diff=diff(x_iter(:,~x_fixed,1:iter),1,3);
        x_iter_diff=squeeze(sum(abs(x_iter_diff),1));
        semilogy([x_iter_diff.',sum(x_iter_diff).'])
        figure(1070)
        for p=1:no_of_par
            subplot(2,no_of_par,p)
            pom_rescale=rescale_parameters(~x_fixed);
            c_lim=quantile(reshape(P_before(:,:,p)./pom_rescale(p),1,[]),[0.1 0.95]);
            c_lim(2)=c_lim(2)+eps;
            imagesc(P_before(:,:,p)./pom_rescale(p),c_lim);
            axis image;colorbar('location','SouthOutside')
            
            subplot(2,no_of_par,no_of_par+p)
            pom_rescale=rescale_parameters(~x_fixed);
            imagesc(P_est_iter(:,:,p)./pom_rescale(p),c_lim);
            axis image;colorbar('location','SouthOutside')
        end
        figure(1050)
        subplot(2,2,3)
        cla
        semilogy(0:iter,(L1_func_iter(1:iter+1,1:no_of_par)))
        hold all
        semilogy(0:iter,(L1_func_iter(1:iter+1,(1+no_of_par):end)))
        title('\gamma\Sigma|\nablap|')
        subplot(2,2,4)
        semilogy((sum(L1_func_iter(1:iter,(no_of_par+1):end),2).'+sum(SSD(:,1:iter),1)))
        title('\sigma^-^2\Sigma||y-C(n,p)||_2^2 + \gamma\Sigma|\nablap|')
        drawnow
    end
end
time_spent=toc
disp(['Total function calls: ' num2str(sum(func_eval(:)))])
% profile viewer
