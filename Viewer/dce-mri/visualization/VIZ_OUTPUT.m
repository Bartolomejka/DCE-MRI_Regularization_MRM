%Function for visualisation
function VIZ_OUTPUT(parametr,parameter_name, frame, mask,curves,curves_name,concetration)
%%
global var
global data
%%
%Figure properities
hlavni=figure('Name',parameter_name,'Resize','off','MenuBar','none','Toolbar','figure','units','pixels','Color',[0.6 0.6 0.6]);
set(hlavni,'Position',[50 50 450 650]);   % size of figure (xposition,yposition,width,height)
set(hlavni,'CloseRequestFcn',@hlavni_zavri);
%Control if some of output figures are open
try existjizobrazene=var.priznak; catch end;
if (exist('existjizobrazene'))
    var.priznak=var.priznak+1;
    var.zavirani{var.priznak}=struct('hlavni',hlavni);
else
    var.priznak=1;
    var.zavirani{var.priznak}=struct('hlavni',hlavni);
end;
slidertrans = uicontrol('Style','slider','Position',[50 90 350 30],...
    'Max',100,'Callback',@slidertrans_Callback);
timetext = uicontrol('Style','text','String',['Transparency: 0%' ],...
    'Position',[165,70,105,20],'HorizontalAlignment','center','BackGroundColor',[0.6 0.6 0.6]);
krivky_button=uicontrol('Style','togglebutton','Position',[40 50 150 25],...
    'Callback',@krivky_button_Callback,'String','Show curves');
contrast_button=uicontrol('Style','togglebutton','Position',[40 20 150 25],...
    'Callback',@contrast_button_Callback,'String','Parameter range');
edit_contrast1=uicontrol('Position',[340 20 60 25],'Style','edit','Callback',@contrast_button_Callback);
maxtext = uicontrol('Style','text','String',['max' ],...
    'Position',[300,20,30,20],'HorizontalAlignment','center','BackGroundColor',[0.6 0.6 0.6]);
edit_contrast2=uicontrol('Position',[240 20 60 25],'Style','edit','Callback',@contrast_button_Callback);
mintext = uicontrol('Style','text','String',['min' ],...
    'Position',[200,20,30,20],'HorizontalAlignment','center','BackGroundColor',[0.6 0.6 0.6]);
plocha = axes('Units','Pixels','Position',[60 150 350 500]);
%%
%Drawing into figure
transp = ones(size(frame)) * 0;
transp(mask==0) = 1;
parametr_max=max(max(parametr));
parametr_min=min(min(parametr));
set(edit_contrast1,'String',parametr_max);
set(edit_contrast2,'String',parametr_min);
transp = flipud(transp);
% frame=double(frame);
% frame=uint8(frame);
g = grs2rgb(flipud(frame), gray(256));
handles.plocha=imagesc(flipud(real(parametr)));colormap(jet)
colorbar;
hold on
h = imshow(g,'InitialMag', 'fit');
hold off
set(h, 'AlphaData', transp)
title(parameter_name)
axis xy
var.mask_length=length(mask);

impixelinfo2(hlavni,handles.plocha); % Bartos - impixelinfo2;

%%
% FUNCTIONS callbacks
    function slidertrans_Callback(source,eventdata)
        % on modification of time slider
        %         impixelinfo;
        %         pomoc=parametr;
        %         aktual_contrast_max=(get(edit_contrast1,'String'));
        %         aktual_contrast_min=(get(edit_contrast2,'String'));
        %         aktual_contrast_max=str2num(aktual_contrast_max);
        %         aktual_contrast_min=str2num(aktual_contrast_min);
        %         if length(aktual_contrast_max)==0 || length(aktual_contrast_min)==0
        %         elseif aktual_contrast_max<aktual_contrast_min
        %             helpdlg('Contrast max has to be higher than contrast min.','Info');
        %         else
        ts=round(get(slidertrans,'Value'));
        set(timetext,'String',['Transparency: ' num2str(ts) '%']);
        transp = ones(size(frame)) * ts/100;
        transp(mask==0) = 1;
        transp = flipud(transp);
        %             g = grs2rgb(flipud(frame), gray);
        %
        %             handles.plocha=imagesc(flipud(parametr),[aktual_contrast_min aktual_contrast_max]);colormap(jet)
        % %             impixelinfo2; deleted by Bartos
        %             colorbar;
        %             hold on
        %             h = imshow(g,'InitialMag', 'fit');
        %             hold off
        set(h, 'AlphaData', transp)
        %             title(parameter_name)
        %             axis xy
        %             parametr=pomoc;
        
        %         end;
        
    end

%%
    function krivky_button_Callback(source,eventdata)
        
        map_size=size(data{end}.data2{1});
        [xovas,yovas]=ginput2(1);
                impixelinfo2;
        xovas=round(xovas);
        %         yovas=length(mask)-round(yovas)+1;
        yovas=size(mask,1)-round(yovas)+1;    % Bartos- oprava v indexaci v masce
        if var.pocitadlo==var.nmr_param_not_curves_counter
            for zaviraci_parametr=1:var.nmr_param_not_curves_counter
                if isempty(var.zavirani_krivek{zaviraci_parametr})
                else
                    %                     delete(var.zavirani_krivek{zaviraci_parametr}.krivkyfig);
                end;
            end;
            var.pocitadlo=0;
        elseif var.pocitadlo==1
            if var.zavirani_krivek{1}.krivkyfig==0
            else
                delete(var.zavirani_krivek{1}.krivkyfig);
            end;
            var.pocitadlo=0;
        else
            var.pocitadlo=0;
        end;
        %          if xovas<0 | yovas<0 | xovas>size(mask,2) | yovas>length(mask) % divne - opravit
        if xovas<0 | yovas<0 | xovas>size(mask,2) | yovas>size(mask,1) % Bartos - oprava v indexaci maskou
            helpdlg('No curves in this area','No curves');
            hodnota=0;
        else
            hodnota=mask(yovas,xovas);
            counter=1;
            namesforsee={};
            valuesforsee={};
            for nmr_parameter=1:length(var.mapy_parametru)
                %                 if var.mapy_parametru{nmr_parameter}==0 % Bartos - proc ?
                if numel(var.mapy_parametru{nmr_parameter})==1 % Bartos - proc ? -> odstraneno
                else
                    namesforsee{counter}=var.nazvy_parametru{nmr_parameter+3};
                    valuesforsee{counter}=var.mapy_parametru{nmr_parameter}(yovas,xovas);
                    counter=counter+1;
                end;
            end;
        end;
        if hodnota==0
            helpdlg('No curves in this area.','No curves');
        else
            if curves{1,1}==0
                nmr_curves=0;
                model_crv_to_graph_concentration=0;
            else
                for nmr_curves=1:length(curves(1,1:end))
                    is_Model_big_curve=strcmp(curves_name{nmr_curves},'Model');
                    is_model_curve=strcmp(curves_name{nmr_curves},'model');
                    if is_Model_big_curve==1 | is_model_curve==1
                        model_crv_to_graph_concentration=1;
                        curve_number_showing=nmr_curves;
                    else
                        try model_crv_to_graph_concentration=model_crv_to_graph_concentration; catch end;
                        if (exist('model_crv_to_graph_concentration'))
                            if model_crv_to_graph_concentration==1
                                krivkyfig=figure(100+nmr_curves-1);
                                clf
                                set(krivkyfig,'CloseRequestFcn',@close_krivky);
                                subplot(2,1,1);
                                axis_actual=0:(data{1}.info.acq.Ts):(length(curves{hodnota,nmr_curves})-1)*(data{1}.info.acq.Ts);
                                plot(axis_actual,curves{hodnota,nmr_curves});xlabel('time [s]');ylabel('amplitude');
                                title(curves_name{1,nmr_curves})
                                %                                 subplot(2,1,2);
                                % Bartos - uitable used instead of text
                                tab_pos=get(h_sub,'Position');
                                tab_pos(2)=tab_pos(2)-1.15*tab_pos(4);
                                
                                namesforsee=[namesforsee(1:(length(namesforsee)/2-1)) namesforsee(end-1:end)];
                                valuesforsee=[[valuesforsee(1:(length(valuesforsee)/2-1)) valuesforsee(end-1:end)];...
                                    [valuesforsee(length(valuesforsee)/2:end-2) 0 0]]';
                                h_tab=uitable('Data',valuesforsee,'RowName',namesforsee,...
                                    'ColumnName',column_names(1:size(valuesforsee,2)),'Units','normalized','Position',tab_pos);
                                text_pos=tab_pos;
                                text_pos(2)=text_pos(2)-0.05;text_pos(4)=0.05;
%                                 position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(map_size(1)-yovas+1)];
                                position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(yovas)];
                                uicontrol('Parent',gcf,'Style','text','String',position_x_y,'Units','normalized','Position',text_pos);
                                %/Bartos
                                %                    text (0.1,0.8,namesforsee);
                                %                    text (0.7,0.8,valuesforsee);
                                %                                 position_x_y=['Position: y=' num2str(yovas) '  x=' num2str(xovas)];
                                %                                 text(0.1,0.4,position_x_y);
                                %                                 axis off;
                                var.pocitadlo=var.pocitadlo+1;
                                var.zavirani_krivek{var.pocitadlo}=struct('krivkyfig',krivkyfig);
                            else
                                krivkyfig=figure(100+nmr_curves);
                                clf
                                set(krivkyfig,'CloseRequestFcn',@close_krivky);
                                subplot(2,1,1);
                                axis_actual=0:(data{1}.info.acq.Ts):(length(curves{hodnota,nmr_curves})-1)*(data{1}.info.acq.Ts);
                                plot(axis_actual,curves{hodnota,nmr_curves});xlabel('time [s]');ylabel('amplitude');
                                title(curves_name{1,nmr_curves})
                                % Bartos - uitable used instead of text
                                tab_pos=get(h_sub,'Position');
                                tab_pos(2)=tab_pos(2)-1.15*tab_pos(4);
                                %                                 namesforsee=[namesforsee(1:(length(namesforsee)/2-1)) namesforsee(end-1:end)];
                                %                                 valuesforsee=[[valuesforsee(1:(length(valuesforsee)/2-1)) valuesforsee(end-1:end)];...
                                %                                     [valuesforsee(length(valuesforsee)/2:end-2) 0 0]]';
                                h_tab=uitable('Data',valuesforsee,'RowName',namesforsee,...
                                    'ColumnName',column_names(1:size(valuesforsee,2)),'Units','normalized','Position',tab_pos);
                                text_pos=tab_pos;
                                text_pos(2)=text_pos(2)-0.05;text_pos(4)=0.05;
                                %                                 position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(map_size(1)-yovas+1)];
                                position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(yovas)];
                                uicontrol('Parent',gcf,'Style','text','String',position_x_y,'Units','normalized','Position',text_pos);
                                %/Bartos
                                %                                 subplot(2,1,2);
                                %                                 text (0.1,0.8,namesforsee);
                                %                                 text (0.7,0.8,valuesforsee);
                                %                                 position_x_y=['Position: y=' num2str(yovas) '  x=' num2str(xovas)];
                                %                                 text(0.1,0.4,position_x_y);
                                %                                 axis off;
                                %                                 model_crv_to_graph_concentration=0;
                                var.pocitadlo=var.pocitadlo+1;
                                var.zavirani_krivek{var.pocitadlo}=struct('krivkyfig',krivkyfig);
                            end;
                        else
                            krivkyfig=figure(100+nmr_curves);
                            clf
                            set(krivkyfig,'CloseRequestFcn',@close_krivky);
                            h_sub=subplot(2,1,1);
                            axis_actual=0:(data{1}.info.acq.Ts):(length(curves{hodnota,nmr_curves})-1)*(data{1}.info.acq.Ts);
                            plot(axis_actual,curves{hodnota,nmr_curves});xlabel('time [s]');ylabel('amplitude');
                            title(curves_name{1,nmr_curves})
                            % Bartos - uitable used instead of text
                            tab_pos=get(h_sub,'Position');
                            tab_pos(2)=tab_pos(2)-1.15*tab_pos(4);
                            ind_std=strncmp('std(',namesforsee,4);
                            namesforsee=namesforsee(~ind_std);
                            valuesforsee=[valuesforsee(~ind_std);...
                                [valuesforsee(ind_std) cell(1,sum(~ind_std)-sum(ind_std))]]';
                            column_names={'\mu','\sigma'};
                            h_tab=uitable('Data',valuesforsee,'RowName',namesforsee,...
                                'ColumnName',column_names(1:size(valuesforsee,2)),'Units','normalized','Position',tab_pos);
                            text_pos=tab_pos;
                            text_pos(2)=text_pos(2)-0.05;text_pos(4)=0.05;
%                             position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(map_size(1)-yovas+1)];
                            position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(yovas)];
                            uicontrol('Parent',gcf,'Style','text','String',position_x_y,'Units','normalized','Position',text_pos);
                            %/Bartos
                            %                    subplot(2,1,2);
                            %                    text (0.1,0.8,namesforsee);
                            %                    text (0.7,0.8,valuesforsee);
                            %                             position_x_y=['Position: y=' num2str(yovas) '  x=' num2str(xovas)];
                            %                             text(0.1,0.4,position_x_y);
                            %                             axis off;
                            model_crv_to_graph_concentration=0;
                            var.pocitadlo=var.pocitadlo+1;
                            var.zavirani_krivek{var.pocitadlo}=struct('krivkyfig',krivkyfig);
                        end;
                    end;
                end;
            end;
        end;
        if hodnota==0
        else
            if model_crv_to_graph_concentration==1
                krivkyfig=figure(100+nmr_curves);
            else
                krivkyfig=figure(100+nmr_curves+1);
            end;
            set(krivkyfig,'CloseRequestFcn',@close_krivky);
            clf
            subplot(2,1,1);
            if model_crv_to_graph_concentration==1
                model_crv_to_graph_concentration=0;
                hold on
                axis_actual_conc=0:(data{1}.info.acq.Ts):(length(concetration{hodnota})-1)*(data{1}.info.acq.Ts);
                plot(axis_actual_conc,concetration{hodnota});xlabel('time[s]');ylabel('amplitude');
                plot(axis_actual_conc,curves{hodnota,curve_number_showing},'r');xlabel('time[s]');ylabel('amplitude');
                hold off
                legend('Concentration','Model');
            else
                
                axis_actual_conc=0:(data{1}.info.acq.Ts):(length(concetration{hodnota})-1)*(data{1}.info.acq.Ts);
                plot(axis_actual_conc,concetration{hodnota});xlabel('time[s]');ylabel('amplitude');
            end;
            title('Concentration')
            % Bartos - uitable used instead of text
            tab_pos=get(h_sub,'Position');
            tab_pos(2)=tab_pos(2)-1.15*tab_pos(4);
            %             namesforsee=[namesforsee(1:(length(namesforsee)/2-1)) namesforsee(end-1:end)];
            %             valuesforsee=[[valuesforsee(1:(length(valuesforsee)/2-1)) valuesforsee(end-1:end)];...
            %                 [valuesforsee(length(valuesforsee)/2:end-2) 0 0]]';
            h_tab=uitable('Data',valuesforsee,'RowName',namesforsee,...
                'ColumnName',column_names(1:size(valuesforsee,2)),'Units','normalized','Position',tab_pos);
            text_pos=tab_pos;
            text_pos(2)=text_pos(2)-0.05;text_pos(4)=0.05;
%             position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(map_size(1)-yovas+1)];
            position_x_y=['Position: x=' num2str(xovas) ' y=' num2str(yovas)];
            uicontrol('Parent',gcf,'Style','text','String',position_x_y,'Units','normalized','Position',text_pos);
            %/Bartos
            %             subplot(2,1,2);
            %             text (0.1,0.8,namesforsee);
            %             text (0.7,0.8,valuesforsee);
            %             position_x_y=['Position: y=' num2str(yovas) '  x=' num2str(xovas)];
            %             text(0.1,0.4,position_x_y);
            %             axis off;
            var.pocitadlo=var.pocitadlo+1;
            var.zavirani_krivek{var.pocitadlo}=struct('krivkyfig',krivkyfig);
        end;
        
    end

%%
    function hlavni_zavri(hObject,eventdata)
        ind=hObject.Number;
        delete(var.zavirani{ind}.hlavni);
        var.zavirani{ind}.hlavni=0;
    end
%%
    function close_krivky(hObject,eventdata)
        if exist('var','var')
            ind=hObject.Number-100;
            delete(var.zavirani_krivek{ind}.krivkyfig);
            var.zavirani_krivek{ind}=[];
        end
    end
%%
    function contrast_button_Callback(source,eventdata)
        %         impixelinfo;
        %         pomoc=parametr;
        aktual_contrast_max=(get(edit_contrast1,'String'));
        aktual_contrast_min=(get(edit_contrast2,'String'));
        
        aktual_contrast_max=str2num(aktual_contrast_max);
        aktual_contrast_min=str2num(aktual_contrast_min);
        if isempty(aktual_contrast_max) || isempty(aktual_contrast_min)
            helpdlg('Enter numbers please','Info');
        elseif aktual_contrast_max<=aktual_contrast_min
            helpdlg('Contrast max has to be higher than contrast min.','Info');
        else
            %            ts=round(get(slidertrans,'Value'));
            %            transp = ones(size(frame)) * ts/100;
            %            transp(mask==0) = 1;
            %            transp = flipud(transp);
            %            g = grs2rgb(flipud(frame), gray);
            %            handles.plocha=imagesc(flipud(parametr),[aktual_contrast_min aktual_contrast_max]);colormap(jet)
            % %            impixelinfo2; Bartos
            %            colorbar;
            %            hold on
            %            h = imshow(g,'InitialMag', 'fit');
            %            hold off
            %            set(h, 'AlphaData', transp)
            %            title(parameter_name)
            %            axis xy
            %            parametr=pomoc;
            set(get(handles.plocha,'Parent'),'CLim',[aktual_contrast_min aktual_contrast_max])
        end;
        
    end
end


