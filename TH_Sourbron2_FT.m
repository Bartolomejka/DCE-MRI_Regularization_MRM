function [H,dH] = TH_Sourbron2_FT(par,w,deriv)
% parametrization [Fp Tc Te alpha tau]


w=w(:);
Fp=par(1);
Tc=par(2);
Te=par(3);
alpha=par(4);
tau=par(5);
t_cut=par(6);
T=Tc+Te.*alpha;
% tau=exp(-1j*par(5)*(0:(length(w)-1)));
phase_shift=exp(-1j*par(5)*w);

p=w*1j;
s=p;

% part to be truncated, i.e. shifted GK (Kety) model
kep=(1-exp(-alpha))./alpha./Te;
H_GK=1./(s+kep);
H_GK=H_GK.*(1-exp(-alpha)).*exp(-kep.*(t_cut-tau-Tc)); % reduce energy
H_GK=H_GK.*exp(-1j*(t_cut-tau)*w); % shift to t_cut

% % part to be truncated, i.e. shifted GK (Kety) model
% H_GK=1./(p/Ktrans+1/ve); % GK model
% H_GK=H_GK*exp(-Ktrans/ve*(t_cut-tau-Tc)); % reduce energy
% H_GK=H_GK.*exp(-1j*(t_cut-tau)*w); % shift to t_cut

H_TH=(T+s.*Tc.*Te).*(s.*Tc+alpha).*(1-exp(-(s.*Tc+alpha))) ./ ...
    (s.*(T+s.*Tc.*Te).*(s.*Tc+alpha)+alpha.*(1-exp(-(s.*Tc+alpha))));

% H_TH=(exp(-(a+Tc*p))-1).*(a+Tc*p).*(a*ve+vp*(p*c*Tc+a))./...
%     (a^2*(exp(-(a+Tc*p))-1)-(a+Tc*p).*p*Tc.*(c*(a+Tc*p)+a));
% % H=H*sqrt(length(w)); % jeste zkontrolovat, jak je to s temi energiemi

H=H_TH-H_GK; % truncation
H=H.*phase_shift; % time-shift
if nargout>1 % if derivatives needed
    if nargin>2
        ind=1:length(deriv);
        ind=ind(deriv);
        dH=zeros(length(w),nnz(deriv));
        dH_GK=zeros(length(w),nnz(deriv));
        for n=1:length(ind)
            switch ind(n)
                case 1
                    %                     dH/dFp
                    dH(:,n)=H./Fp;
                case 2
                    %                     dH/dTc
                    dH(:,n)=((exp(-alpha-Tc.*s)-1.0).*(Te.*s+1.0).*(alpha+Tc.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(s.*(exp(-alpha-Tc.*s)-1.0).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(exp(-alpha-Tc.*s)-1.0).*1.0./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s)).^2.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s).*(s.^2.*(Tc+Te.*alpha+Tc.*Te.*s)+s.*(Te.*s+1.0).*(alpha+Tc.*s)+alpha.*s.*exp(-alpha-Tc.*s))-(s.*exp(-alpha-Tc.*s).*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s));
                    dH_GK(:,n)=(exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*(exp(-alpha)-1.0).^2)./(Te.*alpha.*(s-(exp(-alpha)-1.0)./(Te.*alpha)));
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 3
                    %dH/dTe
                    dH(:,n)=((exp(-alpha-Tc.*s)-1.0).*(alpha+Tc.*s).^2)./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+s.*(exp(-alpha-Tc.*s)-1.0).*1.0./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s)).^2.*(alpha+Tc.*s).^3.*(Tc+Te.*alpha+Tc.*Te.*s);
                    dH_GK(:,n)=(1.0./Te.^2.*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*1.0./(s-(exp(-alpha)-1.0)./(Te.*alpha)).^2.*(exp(-alpha)-1.0).^2)./alpha-(1.0./Te.^2.*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*(exp(-alpha)-1.0).^2.*(Tc-t_cut+tau))./(alpha.*(s-(exp(-alpha)-1.0)./(Te.*alpha)));
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 4
                    %                     dH/dalpha
                    dH(:,n)=((exp(-alpha-Tc.*s)-1.0).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(Te.*(exp(-alpha-Tc.*s)-1.0).*(alpha+Tc.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))-(exp(-alpha-Tc.*s).*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(exp(-alpha-Tc.*s)-1.0).*1.0./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s)).^2.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s).*(-exp(-alpha-Tc.*s)+s.*(Tc+Te.*alpha+Tc.*Te.*s)+alpha.*exp(-alpha-Tc.*s)+Te.*s.*(alpha+Tc.*s)+1.0);
                    dH_GK(:,n)=(exp(-alpha).*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)))./(s-(exp(-alpha)-1.0)./(Te.*alpha))-(exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*(exp(-alpha)-1.0).*((1.0./alpha.^2.*(exp(-alpha)-1.0).*(Tc-t_cut+tau))./Te+(exp(-alpha).*(Tc-t_cut+tau))./(Te.*alpha)))./(s-(exp(-alpha)-1.0)./(Te.*alpha))+exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*1.0./(s-(exp(-alpha)-1.0)./(Te.*alpha)).^2.*(exp(-alpha)-1.0).*((1.0./alpha.^2.*(exp(-alpha)-1.0))./Te+exp(-alpha)./(Te.*alpha));
                    dH(:,n)=dH(:,n)-dH_GK(:,n).*exp(-1j*(t_cut-tau)*w); % shift to t_cut-tau and diff
                    dH(:,n)=dH(:,n).*phase_shift; % time-shift of derivatives
                case 5
                    %                     dH/dtau
                    dH(:,n)=H.*-1j.*w; % time derivative, neni zkonrolovano
                    dH_GK(:,n)=-(s.*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*exp(-s.*(t_cut-tau)).*(exp(-alpha)-1.0))./(s-(exp(-alpha)-1.0)./(Te.*alpha))+(exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*exp(-s.*(t_cut-tau)).*(exp(-alpha)-1.0).^2)./(Te.*alpha.*(s-(exp(-alpha)-1.0)./(Te.*alpha)));
                    dH(:,n)=dH(:,n)+dH_GK(:,n).*phase_shift;
            end
        end
    else
        dH=zeros(length(w),5);
        dH(:,1)=H./Fp; % together with dH_GK/dF
        %                     dH/dTc
        dH(:,2)=((exp(-alpha-Tc.*s)-1.0).*(Te.*s+1.0).*(alpha+Tc.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(s.*(exp(-alpha-Tc.*s)-1.0).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(exp(-alpha-Tc.*s)-1.0).*1.0./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s)).^2.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s).*(s.^2.*(Tc+Te.*alpha+Tc.*Te.*s)+s.*(Te.*s+1.0).*(alpha+Tc.*s)+alpha.*s.*exp(-alpha-Tc.*s))-(s.*exp(-alpha-Tc.*s).*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s));
        dH(:,3)=((exp(-alpha-Tc.*s)-1.0).*(alpha+Tc.*s).^2)./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+s.*(exp(-alpha-Tc.*s)-1.0).*1.0./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s)).^2.*(alpha+Tc.*s).^3.*(Tc+Te.*alpha+Tc.*Te.*s);
        dH(:,4)=((exp(alpha+Tc.*s)-1.0).*(Tc+alpha+Tc.*Te.*s))./(alpha.*(exp(alpha+Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+alpha+Tc.*Te.*s))+((exp(alpha+Tc.*s)-1.0).*(alpha+Tc.*s))./(alpha.*(exp(alpha+Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+alpha+Tc.*Te.*s))+(exp(alpha+Tc.*s).*(alpha+Tc.*s).*(Tc+alpha+Tc.*Te.*s))./(alpha.*(exp(alpha+Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+alpha+Tc.*Te.*s))+1.0./(alpha.*(exp(alpha+Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+alpha+Tc.*Te.*s)).^2.*(exp(alpha+Tc.*s)-1.0).*(alpha+Tc.*s).*(Tc+alpha+Tc.*Te.*s).*(-exp(alpha+Tc.*s)-alpha.*exp(alpha+Tc.*s)+s.*(Tc+alpha+Tc.*Te.*s)+s.*(alpha+Tc.*s)+1.0);
        
        
        dH(:,1)=0; % it is done different way
        dH_GK(:,2)=-(s.*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*exp(-s.*(t_cut-tau)).*(exp(-alpha)-1.0))./(s-(exp(-alpha)-1.0)./(Te.*alpha))+(exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*exp(-s.*(t_cut-tau)).*(exp(-alpha)-1.0).^2)./(Te.*alpha.*(s-(exp(-alpha)-1.0)./(Te.*alpha)));
        dH_GK(:,3)=(1.0./Te.^2.*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*1.0./(s-(exp(-alpha)-1.0)./(Te.*alpha)).^2.*(exp(-alpha)-1.0).^2)./alpha-(1.0./Te.^2.*exp(-((exp(-alpha)-1.0).*(Tc-t_cut+tau))./(Te.*alpha)).*(exp(-alpha)-1.0).^2.*(Tc-t_cut+tau))./(alpha.*(s-(exp(-alpha)-1.0)./(Te.*alpha)));
        dH_GK(:,4)=((exp(-alpha-Tc.*s)-1.0).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(Te.*(exp(-alpha-Tc.*s)-1.0).*(alpha+Tc.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))-(exp(-alpha-Tc.*s).*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s))+(exp(-alpha-Tc.*s)-1.0).*1.0./(alpha.*(exp(-alpha-Tc.*s)-1.0)-s.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s)).^2.*(alpha+Tc.*s).*(Tc+Te.*alpha+Tc.*Te.*s).*(-exp(-alpha-Tc.*s)+s.*(Tc+Te.*alpha+Tc.*Te.*s)+alpha.*exp(-alpha-Tc.*s)+Te.*s.*(alpha+Tc.*s)+1.0);
        
        dH(:,5)=H.*-1j.*w; % time derivative, neni zkonrolovano
        dH_GK(:,5)=(exp(-((exp(-alpha) - 1)*(Tc - t_cut + tau))/(Te*alpha))*exp(-s*(t_cut - tau))*(exp(-alpha) - 1)^2)/(Te*alpha*(s - (exp(-alpha) - 1)/(Te*alpha))) - (s*exp(-((exp(-alpha) - 1)*(Tc - t_cut + tau))/(Te*alpha))*exp(-s*(t_cut - tau))*(exp(-alpha) - 1))/(s - (exp(-alpha) - 1)/(Te*alpha));
        dH(:,5)=dH(:,5)+dH_GK(:,5).*phase_shift;
        
        dH=dH(:,2:4)-dH_GK(:,2:4).*repmat(exp(-1j*(t_cut-tau)*w),1,3); % truncation of derivatives
        
        dH=dH(:,2:4).*repmat(phase_shift,1,3); % time-shift of derivatives
        
    end
    dH=Fp.*dH;
end
H=Fp.*H;