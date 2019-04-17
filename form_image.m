function image=form_image(x,masks)
% reconstruct parameter images from vector of parameters based on mask
no_of_par=size(x,2);
no_of_curves=size(x,1);
[row col]=size(masks{1});
P_pom=zeros(row,col);
image=zeros(row,col,no_of_par);

for p=1:no_of_par
    for c=1:no_of_curves
        P_pom(masks{c})=x(c,p);
    end
    image(:,:,p)=P_pom;
end



%     for p=1:no_of_par
%         P_pom=zeros(row,col);
%         P_pom(phantom_mask_crop)=x_iter(:,ind_to_par(p),end);
%         P_est_iter(:,:,p)=P_pom;
%     end