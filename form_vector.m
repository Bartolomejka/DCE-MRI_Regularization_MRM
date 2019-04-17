function v=form_vector(image,masks)
% form a vector of parameters from parameter maps and masks
% if mask is larger then one pixel, mean of parameters under mask is used

no_of_par=size(image,3);
no_of_curves=size(masks,1);
v=zeros(no_of_curves,no_of_par);

    for c=1:no_of_curves
        for p=1:no_of_par
            im=image(:,:,p);
            v(c,p)=mean(im(masks{c}));
        end
%         [row,col] = find(masks{c});
%         v(c,:)=mean(image(row,col,:));      
    end