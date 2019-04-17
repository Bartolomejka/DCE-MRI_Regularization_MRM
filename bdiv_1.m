%Function to calculate the divergence of a [n,m]-matrix using backward
%differences. The boundary extension mode is 0 for the 0' and n/m'th
%row/colomn
%
%!Remember! grad*=-div i.e. -div is dual operator of gradient
%
%div(v)_[i,j]=v_[i,j]-v_[i-1,j] + v_[i,j]-v_[i,j-1]

function [div_v] = bdiv_1(v,varargin)

n = size(v,1);
m = size(v,2);

if ~isempty(varargin)
    % correction on mask boundary
    phantom_mask=varargin{1};
    div_v = zeros(n,m,2);
    % rows
%     div_v(1,:,1)=v(1,:,1);
%     div_v(n,:,1) = -v(n-1,:,1);
%     div_v(2:n-1,:,1) = v(2:n-1,:,1) - v(1:n-2,:,1);
    div_v(2:n,:,1) = v(2:n,:,1) - v(1:n-1,:,1);
    mask_diff=diff([false(1,size(phantom_mask,2));phantom_mask;false(1,size(phantom_mask,2))]);
    mask1=mask_diff(1:end-1,:);
    mask1=cat(3,mask1,zeros(n,m));
    div_v(mask1>0)=v(mask1>0);
    mask2=mask_diff(2:end,:);
    mask2=cat(3,mask2,zeros(n,m));
    div_v(mask2<0)=div_v(mask2<0)-v(mask2<0);
    % cols
    div_v(:,2:m,2) = v(:,2:m,2) - v(:,1:m-1,2);
    mask_diff=diff([false(size(phantom_mask,1),1),phantom_mask,false(size(phantom_mask,1),1)],1,2);
    mask1=mask_diff(:,1:end-1);
    mask1=cat(3,zeros(n,m),mask1);
    div_v(mask1>0)=v(mask1>0);
    mask2=mask_diff(:,2:end);
    mask2=cat(3,zeros(n,m),mask2);
    div_v(mask2<0)=div_v(mask2<0)-v(mask2<0);
    
    div_v=sum(div_v,3);
    div_v(~phantom_mask)=0;
    
else
    % original Bredies
    div_v = zeros(n,m);
    
    div_v(1,:) = v(1,:,1);
    div_v(n,:) = -v(n-1,:,1);
    div_v(2:n-1,:) = v(2:n-1,:,1) - v(1:n-2,:,1);
    
    
    div_v(:,1) = div_v(:,1) + v(:,1,2);
    div_v(:,m) = div_v(:,m) - v(:,m-1,2);
    div_v(:,2:m-1) = div_v(:,2:m-1) + v(:,2:m-1,2) - v(:,1:m-2,2);
%     div_v(:,2:m-1) = div_v(:,2:m-1) + diff(v(:,1:m-1,2),1,2);
end