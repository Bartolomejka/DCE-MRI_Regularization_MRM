%Function to calculate the gradient of a [n,m] matrix using forward
%differences. The boundary extension is repetition of the last element
%
%(grad_u)_[i,j]=(u_[i+1,j]-u_[i,j],u_[i,j+1]-u_[i,j])

function [grad] = fgrad_1(u,varargin)

[n,m]=size(u);

grad = zeros(n,m,2);

%Set (1,0) derivatives
%grad(n,:,1) = 0;
% grad(1:n-1,:,1) = u(2:n,:) - u(1:n-1,:);
grad(1:n-1,:,1) = diff(u,1,1); % is faster

%Set (0,1) derivatives
%grad(:,m,2) = 0;
% grad(:,1:m-1,2) = u(:,2:m) - u(:,1:m-1);
grad(:,1:m-1,2) = diff(u,1,2); % is faster

if ~isempty(varargin)
    % correction on mask boundary
    phantom_mask=varargin{1};
    % rows
    mask=[~logical(diff(phantom_mask));true(1,size(phantom_mask,2))];
    grad(:,:,1)=grad(:,:,1).*(mask&phantom_mask);
    % cols
    mask=[~logical(diff(phantom_mask,1,2)),true(size(phantom_mask,1),1)];
    grad(:,:,2)=grad(:,:,2).*(mask&phantom_mask);
end