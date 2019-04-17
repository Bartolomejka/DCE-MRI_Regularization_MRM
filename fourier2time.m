function x=fourier2time(X,N)
% x=fourier2time(X,N)
% It takes half-spectrum X, duplicates it and does Inverse FFT to return x
% in time domain, possibly truncated to N samples.

K=min(size(X));
if K==1
    X=X(:);
end
% % LP filter
% W=tukeywin(2*length(X)-1,1); % 0 - rect win, 1 - Hann win
% W=W(length(X):end);
% X=X.*repmat(W,1,K);

% % duplicating of spectrum - for odd number of samples in time domain
% X_full=zeros(2*length(X)-1,K,class(X));
% X_full(1:length(X),:)=X;
% X_full((length(X)+1):end,:)=conj(flipud(X(2:end,:)));

% duplicating of spectrum - for even number of samples in time domain
X_full=zeros(2*(length(X)-1),K,class(X));
X_full(1:length(X),:)=X;

% conjugate symetrization is not needed if ifft is used with 'symmetric' flag
% X_full((length(X)+1):end,:)=conj(X((end-1):-1:2,:)); % zrychleni

% conversion to time domain
x=(ifft(X_full,'symmetric'));

% possible truncating of signal
if nargin>1
    x=x(1:N,:);
end