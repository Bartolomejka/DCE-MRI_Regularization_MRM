% Function for getting parametres maps for perfusion parameters 
function  [curves_name,parameters_vizualization,parameters_tissue,fullmask,curves,concentration] = GetParamMaps(tissue, parameters_tissue)
%%
if strcmp(parameters_tissue{end},'Estimate details')
    parameters_tissue=parameters_tissue(1:end-1);
    tissue=tissue(:,1:length(parameters_tissue));
end;
for j=1:length(parameters_tissue)-3
    if length(tissue{1,j+3})==1
       parameters_vizualization{j} = (zeros(size(tissue{1,2})));
    else
       parameters_vizualization{j} = (zeros(1,1));    
    end;
end;
num_rois = size(tissue,1);
fullmask = zeros(size(tissue{1,2}));
curves_name=cell(1,1);
curves=cell(1,1);
curves{1,1}=0;
for f=1:num_rois
    curves_column=1;
    mask = full(tissue{f,2});
    mask(mask>0)=1;
    for g=1:length(parameters_tissue)-3
        if length(tissue{1,g+3})==1
          parameters_vizualization{g}(mask>0)=0;
          pom=tissue{f,g+3};
          pom(isnan(pom))=0; % Bartolomej - to remove NaNs possibly present aftre J'J matrix inversion
          parameters_vizualization{g} = parameters_vizualization{g} + mask*pom;
        else
          curves{f,curves_column}=tissue{f,g+3};    
          curves_name{1,curves_column}=parameters_tissue{1,g+3};
          curves_column=curves_column+1;
        end;
    
    end;
    fullmask(mask>0) = f;
end;

concentration=tissue(1:end,1);
