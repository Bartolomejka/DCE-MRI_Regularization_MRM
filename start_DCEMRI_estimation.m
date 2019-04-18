% Starting script of the spatially regularized perfusion parameter
% estimation
% Author: Michal Bartos - UTIA CAS 2019
%% select input data data by uncomenting
% RatPhantom with different noise realizations
folder='Data/';
% Synthetic rat phantom
% file_template='SyntheticEchoesOlivier_v04_v10_sec01_inp_con_tis_190416.mat';
% Preclinical dataset
% file_template='P13_sec01_con_tis_190417.mat';
% Clinical dataset
file_template='0009_KL_AN3_sec02_SVAL_CLANEK_all_SNRtr0_inp_con_tis_190417.mat';
% -------------------------------------------------------------------------

filenames=dir([folder,file_template]);
for file_number=1:length(filenames) % batch processing
    pause(1)
    close all
    disp(['File number: ' num2str(file_number) '/' num2str(length(filenames))]);
    % data_name=[filenames(file_number).folder '/' filenames(file_number).name];
    data_name=[folder filenames(file_number).name]; % old Matlab
    
    load(data_name);
    [data_path, data_name,~] = fileparts(data_name);
    
    % % % data reduction
%     tissue=tissue(1:round((end/6)),:);
    
    % mask
    [row,col]=size(data2{1});
    phantom_mask_crop=false(row,col);
    no_of_curves=size(tissue,1);
    for n=1:no_of_curves
        phantom_mask_crop(tissue{n,2})=true;
    end
    
    % visualize
    figure
    imagesc(medfilt2(data2{1},[3 3]))
    axis image;colormap gray
    figure
    imagesc(phantom_mask_crop)
    axis image;colormap gray
    figure
    plot(tissue{1,1})
    %% load aif from the input file
    aif=info.aif.aif{1};
    
    x_aif=0;
    info.aif.x_aif=x_aif;
    info.aif.model='sampled_AIF';
    
    % % visualize
    figure
    plot(aif)
    %% setup starting values
    % time
    Ts=info.acq.Ts/60; % sampling interval
    no_of_timesamples_or=length(tissue{1,1});
    t=(0:no_of_timesamples_or-1)*Ts; % time axis
    
    % normalization of signal by effective mean aif sample area
    norm_coeff=mean(aif);
    aif_norm_coeff=median(aif);
    % data normalizatin coefficient
    trf_norm_coeff=double(cellfun(@median,tissue(:,1)));
    trf_norm_coeff=quantile(trf_norm_coeff,0.75);
    
    % estimate noise std
    std_cellfun=@(x)function_stdEst2D(x); % MAD estimator of the noise
    noise_std_est=cellfun(std_cellfun,tissue(:,1));
    SNR_cellfun=@(x,y)20*log10(abs(mean(x))/y);
    SNR_est=cellfun(SNR_cellfun,tissue(:,1),mat2cell(noise_std_est,ones(size(noise_std_est))));
    
    SNR_map=form_image(SNR_est,tissue(:,2));
    noise_std_map=form_image(noise_std_est,tissue(:,2));
    
    figure;imagesc(SNR_map);colorbar('vert');axis image
    title('SNR [dB]')
    figure;imagesc(noise_std_map);colorbar('vert');axis image
    noise_std=mean(noise_std_est)
    title(['\sigma_N=' num2str(noise_std)])
    noise_std=double(noise_std_est);
    
    % remove curves with unsuficient SNR
    SNR_limit=0;%-200;%-5;%
    tissue(or(SNR_est<SNR_limit,isnan(SNR_est)),:)=[];
    if length(noise_std)>1
        noise_std(or(SNR_est<SNR_limit,isnan(SNR_est)),:)=[];
    end
    % update mask etc.
    % mask
    [row,col]=size(data2{1});
    phantom_mask_crop=false(row,col);
    no_of_curves=size(tissue,1);
    for n=1:no_of_curves
        phantom_mask_crop(tissue{n,2})=true;
    end
    
    % normalize noise std
    noise_std=noise_std/trf_norm_coeff;
    noise_std_map=form_image(noise_std,tissue(:,2));
    
    %% IRF model
    
    % TH - Sourbron - new parametrization
    model_irf=@TH_Sourbron2_FT;
    no_of_irf_par=5;
    irf_lb=[1e-3 Ts Ts 1e-4 -0.5]; % lower bound
    irf_ub=[100 3 100 3 1]; % upper bound
    % generate starting points
    % % x=[Fp Tc Te alpha tau]
    x0_irf=[1 0.1 2.5 0.4 0]; % starting point
    x0_irf=repmat(x0_irf,no_of_curves,1); % equal for all voxels
    x_fixed_irf=logical([0 0 0 0 0]); % full TH model
    x_irf_rescale_coeff=[1 10 1 10 10]; % scaling of the parameters -> convergence
    
    x_irf_label={'F_p' 'T_c' 'T_e' '\alpha' '\tau'};
    % /TH - Sourbron2
    % TH - Sourbron2, truncated - speedup, see Bartos - IUPESM Conference Prague
    no_of_irf_par=6;
    
    irf_lb=[irf_lb 0];
    irf_ub=[irf_ub Inf];
    x0_irf=[x0_irf repmat(t(end),no_of_curves,1)];
    x_fixed_irf=[x_fixed_irf true];
    x_irf_rescale_coeff=[x_irf_rescale_coeff 1];
    x_irf_label=[x_irf_label {'t_cut'}];
    x_irf_norm_coef=[trf_norm_coeff/aif_norm_coeff 1 1 1 1 1];
    % /TH - Sourbron2, truncated
    
    %% AIF model
    % use non-blind for sampled AIF
    x0_aif=x_aif;
    x_aif_rescale_coeff=1; %
    aif_lb=0;
    irf_lb(~strcmp(x_irf_label,'tau')==0)=0; % because of sampled aif
    aif_ub=0;
    x_aif_norm_coef=1;
    no_of_aif_par=1;
    x_aif_label={'fake'};
    
    x_fixed_aif=true(1);
    global aif_database
    % aif_database = struct('model',[0 aif],'t_cut',t(end),'Ts',Ts,'x_aif',[],'index',zeros(1,2^12,'uint16'),'data',[]);
    aif_database = struct('model',aif/aif_norm_coeff,'t_cut',t(end),'Ts',Ts,'x_aif',[],'index',zeros(1,2^12,'uint16'),'data',[]);
    model_aif=@aif_nonblind_real;
    %% general preparations
    x_label=[x_aif_label x_irf_label];
    
    rescale_parameters=[x_aif_rescale_coeff x_irf_rescale_coeff];
    x_fixed=[x_fixed_aif x_fixed_irf];
    x_norm_coef=[x_aif_norm_coef x_irf_norm_coef];
    
    pos_val={1:size(x0_irf,1),1:size(x0_aif,1)};
    x0_comb=generate_starting_points(pos_val);
    x_start=[x0_aif(x0_comb(:,2),:) x0_irf(x0_comb(:,1),:)];
    
    x_start_rescaled=x_start.*repmat(rescale_parameters,size(x_start,1),1);
    
    lb=[aif_lb irf_lb];
    ub=[aif_ub irf_ub];
    
    lb=lb.*rescale_parameters;
    ub=ub.*rescale_parameters;
    
    % no_of_starting_points=size(x_start,1);
    
    model.aif=struct('model',model_aif,'no_of_par',no_of_aif_par);
    model.irf=struct('model',model_irf,'no_of_par',no_of_irf_par);
    
    x0=x_start_rescaled(:,~x_fixed);
    
    lb=lb(~x_fixed);
    ub=ub(~x_fixed);
    
    x_aif=x0_aif;
    % x_irf=x0_irf;
    
    data_separate = tissue;
    data_separate(:,1)=cellfun(@(x)x.'/trf_norm_coeff,data_separate(:,1),'UniformOutput',false);
    
    %% model check
    [Cfit_pom,dCfit,aif_irf]=C_fit_FT(unique(x_start_rescaled,'rows'),t,model,rescale_parameters,x_fixed);
    
    figure
    plot(t,Cfit_pom)
    hold on
    plot(t(1:length(data_separate{1,1})),data_separate{1,1})
    figure
    plot(t,dCfit(:,:,1))
    legend(x_label(~x_fixed));
    
    figure
    subplot(1,2,1)
    plot(t,aif_irf.aif)
    title('aif')
    subplot(1,2,2)
    plot(t,aif_irf.irf)
    title('irf')
    
    %% run minimization
    proximalLM_ChP
    %% final computations and visualization
    figure
    for n=1:size(P_est_iter,3)
        subplot(1,no_of_par,n)
        pom_rescale=rescale_parameters(~x_fixed);
        imagesc(P_est_iter(:,:,n)./pom_rescale(n));
        axis image;colorbar('location','SouthOutside')
        title(x_label(ind_to_par(n)))
    end
    figure
    plot(squeeze(x_iter(1,~x_fixed,:)).'./repmat(rescale_parameters(~x_fixed),iter+1,1))
    legend(x_label(~x_fixed))
    figure
    subplot(2,2,1)
    semilogy(0:(iter-1),(sum(SSD(:,1:iter),1)))
    title('\Sigma_i\sigma_i^-^2||y_i-C(\cdot,p_i)||_2^2')
    subplot(2,2,2)
    l_mean=mean(lambda,1);
    l_min=min(lambda,[],1);
    l_max=max(lambda,[],1);
    semilogy([l_min;l_mean;l_max].')
    title('\lambda')
    subplot(2,2,3)
    plot(log10(L1_func_iter))
    title('log_1_0(\sigma^2*\gamma*\Sigma|\nablap|)')
    subplot(2,2,4)
    plot(log10(sum(L1_func_iter(:,(1+no_of_par):end),2).'+sum(SSD(:,1:size(L1_func_iter,1)),1)))
    title('log_1_0(SSD + \sigma^2*\gamma*\Sigma|\nablap|)')
    
    legend show
    
    [Cfit_start deriv_start aif_irf_start]=C_fit_FT(x_norm_coef.*x_iter(1,:,1),t,model,rescale_parameters,x_fixed);
    x_est=zeros(size(x_iter,1),size(x_iter,2));
    x_est(:,:)=x_iter(:,:,end-1);
    % return normalization by aif mean
    x_est=repmat(x_norm_coef,[size(x_est,1),1]).*x_est;
    COVB=cell(no_of_curves,1);
    clear aif_est Cfit_est irf_est
    for c=1:no_of_curves
        [Cfit_est(:,c) deriv_est aif_irf_est]=C_fit_FT(x_est(c,:),t,model,rescale_parameters,x_fixed);
        % unnormalize all
        aif_est(:,c)=single(aif_irf_est.aif.*aif_norm_coeff);
        irf_est(:,c)=single(aif_irf_est.irf);
        Cfit_est(:,c)=single(Cfit_est(:,c)*aif_norm_coeff);
        deriv_est=deriv_est*aif_norm_coeff;
        
        [Q,R]=qr(deriv_est,0);
        Rinv=R\eye(no_of_par);
        COVB{c}=(Rinv*Rinv')*sum((Cfit_est(:,c)-trf_norm_coeff*data_separate{c,1}).^2)/(length(t)-no_of_par);
        COVB{c}=diag(1./rescale_parameters(~x_fixed).^2)*COVB{c};
    end
    ind_to_see=randi(no_of_curves,1,5);
    figure
    hold on
    plot(t,trf_norm_coeff.*cat(2,data_separate{ind_to_see,1}))
    plot(t,Cfit_start*aif_norm_coeff,'g');
    plot(t,Cfit_est(:,ind_to_see),'k');
    title(['TRF: ' mat2str(ind_to_see)])
    
    Cfit_est=single(Cfit_est);
    
    % % % Visualize curves - not good for a lot of curves
    % figure
    % plot(t,aif_irf_start.aif,'g');
    % hold on
    % % plot(t,aif,'r');
    % plot(t,aif_est,'k');
    % title('aif')
    % figure
    % plot(t,aif_irf_start.irf,'g');
    % hold on
    % % plot(t,irf,'r');
    % plot(t,irf_est,'k');
    % title('irf')
    % figure
    % plot(t,cat(1,tissue{:,1}).',':r')
    % hold on
    % plot(t,Cfit_start,'g');
    % % plot(t,Cfit,'r');
    % plot(t,Cfit_est,'k');
    % title('c(t)')
    
    x_est=single(x_est./repmat(rescale_parameters,no_of_curves,1));
    %% convert result to perfan scheme
    clear x_est_std_full SSD_last
    if size(data_separate,2)>2
        data_separate(:,4:end)=[];
    end
    % conversion to all parameters
    if strcmp(func2str(model.irf.model),'TH_Sourbron2_FT')
        [x_est_full COVB_full x_irf_label_all]=TH_sourbron22complete_parameters(x_est(:,(no_of_aif_par+1):end-1),COVB,x_fixed(no_of_aif_par+1:end-1));
    else
        error('Unknow model.')
    end
    
    x_est_full=[x_est(:,1:no_of_aif_par) x_est_full];
    x_label_all=[x_aif_label x_irf_label_all];
    % std labels
    x_label_all=[x_label_all cellfun(@(x)['std(' x ')'],x_label_all,'Uniformoutput',false)];
    
    % std from covariances
    x_est_std_full=cell2mat(cellfun(@(x)sqrt(abs(diag(x))).',COVB_full,'UniformOutput',false));
    x_est_std_full=[zeros(size(x_est_std_full,1),no_of_aif_par) x_est_std_full];
    % estimated parameters and std maps
    data_separate(:,(1:(size(x_est_full,2)+size(x_est_std_full,2)))+3)=...
        num2cell([x_est_full x_est_std_full]);
    % data_separate(:,(1:size(x_est_full,2))+3)=num2cell(x_est_full); % without covariances
    % SSD, SNR
    SSD_last=sum((Cfit_est-trf_norm_coeff*[data_separate{:,1}]).^2,1);
    data_separate=[data_separate ...
        num2cell(SSD_last.') ...
        num2cell(10*log10(sum(Cfit_est.^2)./SSD_last).')];
    % estimated curves
    data_separate=[data_separate ...
        mat2cell(aif_est.',ones(1,size(Cfit_est,2))) ...
        mat2cell(irf_est.',ones(1,size(Cfit_est,2))) ...
        mat2cell(Cfit_est.',ones(1,size(Cfit_est,2)))];
    
    % saving
    info.perfan=struct;
    info.perfan.params=[x_aif_label x_irf_label];
    info.perfan.x_start=x_start;
    info.perfan.x_fixed=x_fixed;
    info.perfan.lb=lb;
    info.perfan.ub=ub;
    info.perfan.rescale_parameters=rescale_parameters;
    info.perfan.model=model;
    info.perfan.iterations_proximal=size(x_iter,3)-1;
    info.perfan.lambda_base=lambda_base;
    info.perfan.lambda_min=lambda_min;
    info.perfan.lambda_max=lambda_max;
    info.perfan.noise_std=noise_std; % po normalizaci!
    
    info.perfan.aif_norm_coeff=aif_norm_coeff;
    info.perfan.trf_norm_coeff=trf_norm_coeff;
    
    info.perfan.method='Levenberg-Marquardt';
    
    info.acq.Ts=Ts;
    
    info.tissue.descr=cat(2,{'Curve [-]' 'Mask' 'Mask coordinates'},x_label_all,{'SSD','SNR','AIF', 'IRF', 'Model','Estimate details'}); % without std estimates
    
    info.codewords = {'par'};
    if TV_regularization
        gamma_text=mat2str(gamma_no_scale(~x_fixed),2);
        notes_text=['TV_it' num2str(size(x_iter,3)-1) '_g' gamma_text(2:end-1)];
        info.perfan.TV_gamma=gamma_no_scale(~x_fixed);
    else
        notes_text=['noTV_it' num2str(size(x_iter,3)-1)];
    end
    notes_text=[notes_text 'start_' num2str(size(unique(x_start_rescaled,'rows'),1)) 'pt'];
    info.notes = {notes_text};
    info.perfan.iterations.x_iter=x_iter;
    info.perfan.iterations.SSD_iter=SSD;
    info.perfan.iterations.L1_iter=L1_func_iter;
    info.perfan.iterations.func_eval=func_eval;
    info.perfan.time_spent=time_spent;
    
    [filename info]=name_file(info);
    tissue = data_separate;
    
    % correct normalization
    tissue(:,1)=cellfun(@(x)x*trf_norm_coeff,data_separate(:,1),'UniformOutput',false);
    
    % save(filename,'info','tissue','data2','-v7.3'); % larger files
    save(filename,'info','tissue','data2','-v7');
    disp(['Data were saved as: ' filename])
    % //saving
end % // batch processing
