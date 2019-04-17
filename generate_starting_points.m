function [start_pts]=generate_starting_points(parameters_matrix)
% generate matrix of starting points
% for matrix(double) (every row is original starting point)
%   from matrix of posible values: parameters_matrix=[a1 a2 a3;b1 b2 b3] ->
%   start_pts=[a1 b1;a1 b2;a1 b3; a2 b1; a2 b2; a2 b3;a3 b1;a3 b2;a3 b3]
% for vector(cell) from matrix of posible values: parameters_matrix={[a1 a2
%   a3] [b1 b2]) ->
%   start_pts=[a1 b1;a1 b2; a2 b1; a2 b2; a3 b1;a3 b2]

if iscell(parameters_matrix)
    no_of_parameters=length(parameters_matrix);
    no_of_combination=ones(1,no_of_parameters+1);
    for n=no_of_parameters:-1:1
        no_of_values=length(parameters_matrix{n});
        no_of_combination(n)=no_of_combination(n+1)*no_of_values;
    end
    start_pts=zeros(no_of_combination(1),no_of_parameters);
    for n=1:no_of_parameters
        p=parameters_matrix{n};
        p=repmat(p',1,no_of_combination(n+1))';
        p=p(:);
        start_pts(:,n)=repmat(p,no_of_combination(1)/length(p),1);
    end

else
    [M N]=size(parameters_matrix);
    no_of_combination=N^M;
    start_pts=zeros(no_of_combination,M);
    for n=1:M
        p=parameters_matrix(n,:);
        p=repmat(p',1,no_of_combination/N^n)';
        start_pts(:,n)=repmat(p(:),N^(n-1),1);
    end
end
end