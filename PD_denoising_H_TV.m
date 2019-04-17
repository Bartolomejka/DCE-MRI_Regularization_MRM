function [x reg_term]=PD_denoising_H_TV(y,H,varargin)
% Modified TV denoising for perfusion estimation based on paper:
% A. Chambolle and T. Pock, "A first-order primal-dual algorithm for convex
% problems with applications to imaging," J. Math. Imaging Vis., vol. 40,
% no. 1, pp. 120-145, 2011. Algorithm 1. Gradient and divergence operators
% are from Christian Bredies (TGV JPEG reconstruction);
%
% Michal Bartos, UTIA CAS 2017

% debug=true;
debug=false;

if nargin>2
    % remap operators if mask is present
    fgrad1=@(x)fgrad_1(x,varargin{1});
    bdiv1=@(x)bdiv_1(x,varargin{1});
else
    fgrad1=@(x)fgrad_1(x);
    bdiv1=@(x)bdiv_1(x);
end

if nargin<=3
    fgrad=@(x)fgrad_m(x);
    bdiv=@(x)bdiv_m(x);
elseif nargin==4
    reg=varargin{2};
    % remap operators if mask and regularization is present
    fgrad=@(x)fgrad_m_r(x);
    bdiv=@(x)bdiv_m_r(x);
end

norm_A=sqrt(max(reg)*8); % norm of the operator, 8 for forward differences operator

% Initialization of variables
u=y;
v=fgrad(u);
x=u;

% compute initial value of the regularization term
reg_term(1,:)=sum(sum(sqrt(sum(v.^2,4)),2),1);

theta=0.1;
tau=1/norm_A;
N_iter=200;
beta=5*nnz(varargin{1})/max(reg_term(1,:))*1e-1;
sigma=tau*beta;

% check convergence limit
disp(['\tau*\sigma*||A||^2<1: ' num2str(tau*sigma*norm_A.^2)])

v=0*v; % somehow better init.

aux_I=eye(size(H{1},1));
itaukappa_I=cellfun(@(x)(tau*x+aux_I)\aux_I,H,'UniformOutput',false);

for n=1:N_iter
    if debug
        L1_image=sqrt(sum(fgrad(x).^2,4));
        Data_image=0.5*spec_mult(H,x-y).*(x-y);
        crit_fce(n)=sum(Data_image(:))+sum(L1_image(:));
    end
    % Traditonall PD algorithm with fixed step-size
    %   1) projection step
        v_plus_grad_x=v+sigma*fgrad(x);
        v_max=max(1,sqrt(sum(v_plus_grad_x.^2,4))); % isotropic - projection onto ball
        v_new(:,:,:,1)=v_plus_grad_x(:,:,:,1)./v_max;
        v_new(:,:,:,2)=v_plus_grad_x(:,:,:,2)./v_max;
    %   2) quadratic step
        u_new=u-tau*(-bdiv(v_new))+spec_mult(H,tau*y); % spec_mult(kappa,tau*d) can be precomputed
        u_new=spec_mult(itaukappa_I,u_new);
    %   3) update variables
        x=u_new+theta*(u_new-u);
        u=u_new;
        v=v_new;
end

reg_term(2,:)=sum(sum(sqrt(sum(fgrad(x).^2,4)),2),1);

if debug
    figure(2002)
    hold all
    plot(crit_fce)
    h_legend=legend;
    h_legend.String{end}=['\theta=' num2str(theta) ', \tau=' num2str(tau(1))...
        ', \sigma=' num2str(beta*tau(1)) ', \beta=' num2str(beta)];
    %     legend(h_legend,'show')
    drawnow
end

    function fg=fgrad_m(x)
        size_x=size(x);
        fg=zeros([size_x 2]);
        for m=1:size(x,3)
            fg(:,:,m,:)=fgrad1(x(:,:,m));
        end
    end

    function fd=bdiv_m(x)
        size_x=size(x);
        fd=zeros(size_x(1:3));
        for m=1:size(x,3)
            fd(:,:,m)=bdiv1(x(:,:,m,:)); % should be bdiv1(squeeze(x(:,:,m,:))) but it does not matter and squeeze is slower
        end
    end

    function fg=fgrad_m_r(x)
        size_x=size(x);
        fg=zeros([size_x 2]);
        for m=1:size(x,3)
            fg(:,:,m,:)=reg(m)*fgrad1(x(:,:,m));
        end
    end

    function fd=bdiv_m_r(x)
        size_x=size(x);
        fd=zeros(size_x(1:3));
        for m=1:size(x,3)
            fd(:,:,m)=reg(m)*bdiv1(x(:,:,m,:)); % should be bdiv1(squeeze(x(:,:,m,:))) but it does not matter and squeeze is slower
        end
    end

    function y=spec_mult(mtx,x)
        mask_vec=repmat(varargin{1}(:),size(x,3) ,1);
        x2=reshape(x(mask_vec),length(mtx),[]).';
        y1=zeros(size(mtx{1},1),length(mtx));
        for m=1:length(mtx)
            y1(:,m)=mtx{m}*x2(:,m);
        end
        y2=reshape(y1.',1,[]);
        y=zeros(numel(x),1);
        y(mask_vec)=y2;
        y=reshape(y,size(x));
    end

    function y=solve_system(mtx,x)
        mask_vec=repmat(varargin{1}(:),size(x,3) ,1);
        x2=reshape(x(mask_vec),length(mtx),[]).';
        y1=zeros(size(mtx{1},1),length(mtx));
        for m=1:length(mtx)
            y1(:,m)=mtx{m}\x2(:,m);
        end
        y2=reshape(y1.',1,[]);
        y=zeros(numel(x),1);
        y(mask_vec)=y2;
        y=reshape(y,size(x));
    end       
end