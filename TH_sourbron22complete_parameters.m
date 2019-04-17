function [param_complete,covariance_matrices,label_out]=TH_sourbron22complete_parameters(par_matrix,covariance_matrices,x_fixed)
% [parametry,K,label_out]=TH_sourbrno22complete_parameters(par_matrix)
% input:
%   par_matrix=[Fp Tc Te alpha tau] - original values or
%   par_matrix=[Fb Fp T Tc Te E PS vp ve Ktrans kep tau] - complete parameters
% output:
%   param_complete=[F_blood F_plasma T Tc Te E PS v_p v_e K_trans Kep tau]
%   param_complete=[Fp Tc Te alpha tau] - parameters as input to the related file,
%   i.e.: TH_Sourbron2_FT.m
%   K - transformation matrix for extension of original covariance
%       matrix, i.e. COVB_complete=K*COVB*K'
%   label_out={'F_b' 'F_p' 'T' 'T_c' 'T_e' 'E' 'PS' 'v_p' 'v_e' 'K^t^r^a^n^s' 'k_e_p' '\tau'}

switch size(par_matrix,2)
    case 5
        Hct=0.42;
        param_complete=zeros(size(par_matrix,1),12);
        
        if (nargout>1&&nargin>1)
            %             K=zeros([12,5,size(par_matrix,1)]);
            K=zeros(12,5);
            if ~iscell(covariance_matrices)
                covariance_matrices={covariance_matrices};
            end
            if (size(par_matrix,2)~=size(covariance_matrices{1},2))
                covariance_matrices=cellfun(@(x)covmtx_fixed(x,x_fixed),covariance_matrices,'UniformOutput',false);
            end
        end
        
        if nargout>2
            label_out={'F_b' 'F_p' 'T' 'T_c' 'T_e' 'E' 'PS' 'v_p' 'v_e' 'K^t^r^a^n^s' 'k_e_p' '\tau'};
        end
        
        for n=1:size(par_matrix,1)
            par=par_matrix(n,:);
            Fp=par(1);
            Tc=par(2);
            Te=par(3);
            alpha=par(4);
            tau=par(5);
            
            F_blood=Fp/(1-Hct);
            E=1-exp(-alpha);
            vp=Fp*Tc;
            T=alpha*Te+Tc;
            ve=Fp*alpha*Te; %ve=Fp*(T-Tc);
            Ktrans=Fp*E; % Ktrans=Fp(1-exp(-alpha))
            kep=Ktrans/ve; % kep=(1-exp(-alpha)/(alpha*Te)
            PS=alpha.*Fp;
            
            param_complete(n,:)=[F_blood Fp T Tc Te E PS vp ve Ktrans kep tau];
            
            if (nargout>1&&nargin>1)
                % par_matrix=[Fp T Tc Te tau] - original values
                % par_out=[F_blood F_plasma T Tc Te E PS v_p v_e K_trans Kep tau]
                K(1,1)=1/(1-Hct); %[dFb/dFp]
                K(2,1)=1; %[dFp/dFp]
                K(3,2:4)=[1 alpha Te]; %[dT/d(Tc dTe dalpha)]
                K(4,2)=1; %[dTc/dTc]
                K(5,3)=1; %[dTe dTe]
                K(6,4)=exp(-alpha); %[dE/dalpha]
                K(7,[1 4])=[alpha Fp];% [dPS/d(Fp alpha)]
                K(8,[1 3])=[Tc Fp];%[dvp/d(Fp Tc)]
                K(9,[1 3 4])=[alpha*Te Fp*alpha Fp*Te];% [dve/d(Fp Te alpha)]
                K(10,[1 4])=[E Fp*exp(-alpha)];%[dKtrans/d(Fp alpha)]
                K(11,3:4)=[(exp(-alpha)-1)/(alpha*Te^2) exp(-alpha)/(alpha*Te)];% [dkep/d(Te alpha)]
                K(12,5)=1; %[dtau/dtau]
                
                covariance_matrices{n}=K*covariance_matrices{n}*K.';
            end
        end
    case 12
        label_out={'F_p' 'Tc' 'Te' '\alpha' '\tau'};
        param_complete=par_matrix(:,[2 4 5 7 12]);
        param_complete(:,4)=param_complete(:,4)./param_complete(:,1); % alpha=PS/Fp
        
        if (nargout>1&&nargin>1) % tohle tady k nicemu neni
            K=zeros(12,5);
            K(2,1)=1;
            K(8,2)=1;
            K(9,3)=1;
            K(7,4)=1;
            K(12,5)=1;
        end
    otherwise
        error('Wrong number of the parameters in the input matrix.')
end
end
function [K]=covmtx_fixed(K,x_fixed)
I=zeros(length(x_fixed));
I(~x_fixed,~x_fixed)=K;
K=I;
end