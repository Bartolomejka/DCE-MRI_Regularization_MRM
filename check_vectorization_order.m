function reformat=check_vectorization_order(tissue)
% Checks, whether the tissue curves are stored in a correct order (column
% first), otherwise returns 'false' for further processing.
% Inputs: tissue - cellarray with defined format; data2 - cellarray with
% background image

[grid_x grid_y]=meshgrid(1:size(tissue{1,2},2),1:size(tissue{1,2},1));
x=cellfun(@(x)grid_x(x),tissue(:,2));
y=cellfun(@(x)grid_y(x),tissue(:,2));

if ~any(diff(x)<0)
    % column first
    reformat=false;
elseif ~any(diff(y)<0)
    % row first
    reformat=true;
    warning('Reordering of the data (horizontal to vertical) need to be performed.')
else
    error('Random ordering of curves - regularization not possible')
end