function [aif daif] = aif_nonblind_real(x_aif,w,t_cut,aif_ind)
% AIF for non-blind case - real sampled
% It is stored as a database of allready computed AIFs
% It needs to have prepared global variable 'aif_database'
% in the form:
%     aif_database.model=@aif_model_function
%     aif_database.index=[];
%     aif_database.data=[];
%     aif_database.t_cut=t_cut;
%     aif_database.x_aif=x_aif;
%     aif_database.Ts=Ts;

global aif_database



if any(aif_ind)
    error('Sampled AIF is not possible in blind estimation.')
    [aif daif] = aif_database.model(x_aif,w,t_cut,aif_ind);
else
    % check database correctness
    if (t_cut~=aif_database.t_cut)||any(aif_database.x_aif~=x_aif)||(aif_database.Ts~=pi/w(end))
        aif_database.index=uint16([]);
        aif_database.data=[];
        aif_database.t_cut=t_cut;
        aif_database.x_aif=x_aif;
        aif_database.Ts=pi/w(end);
    end
    
    N_w=length(w);
    N_index=length(aif_database.index);
    if N_index<N_w
        to_add=zeros(1,N_w-N_index,class(aif_database.index));
        aif_database.index=[aif_database.index to_add];
        
        %         aif = aif_database.model(x_aif,w,t_cut,aif_ind);
        aif = AIF(aif_database.model,N_w);
        aif_database.data=[aif_database.data {aif}];
        aif_database.index(N_w)=length(aif_database.data);
    elseif aif_database.index(N_w)==0
        aif = AIF(aif_database.model,N_w);
        aif_database.data=[aif_database.data {aif}];
        aif_database.index(N_w)=length(aif_database.data);
    else
        aif = aif_database.data{aif_database.index(N_w)};
    end
    daif=[];
end

    function AIF=AIF(aif,N_w)
        AIF = fft(aif,2*(N_w-1));
        AIF = AIF(1:N_w);
    end
end