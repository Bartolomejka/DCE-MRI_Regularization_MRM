%Script for GUI visualization of perfusion parameters, prepared for main
%GUI
%which is above this module. Script needs functions, that are in our SVN
%structure
function varargout = VISUALIZATION(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @VISUALIZATION_OpeningFcn, ...
                   'gui_OutputFcn',  @VISUALIZATION_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end

function VISUALIZATION_OpeningFcn(hObject, eventdata, handles, varargin)
%%
handles.output = hObject;
handles.priznak=0;
% Update handles structure
guidata(hObject, handles);
addpath('../../dce-mri/visualization');%adding path for functions that this script needs
addpath('../../dce-mri/visualization/ginput_impixelinfo');
end

function varargout = VISUALIZATION_OutputFcn(hObject, eventdata, handles)
%%
% Get default command line output from handles structure
varargout{1} = handles.output;
set(handles.figure1,'CloseRequestFcn',@hlavni_close);
end

function opening_button_Callback(hObject, eventdata, handles)
%%
global data
global var
% Path to files start
global path_data;
% var.pathname=path_data;
% Path to files end
data={}; 
try zaznam_parametru=data{end}.info.perfan.params;catch end;
try zaznam_zdroje=data{end}.info.source;catch end;
try zaznam_tissue=data{end}.tissue;catch end;
if size(data,1)>0 & (exist('zaznam_parametru')) & (exist('zaznam_zdroje')) & (exist('zaznam_tissue'))                     
   b1=questdlg('There are some loaded data. Do you want clear current data and load new dataset?','Warning','Yes','No','No');
else
   b1='Yes';
end;
switch b1
  case 'Yes'
      data={};  
      set(handles.text1,'String','Source');
%       pathfiles = '';
%                 if exist('var','var')
%                     if isfield(var,'pathname')
                        pathfiles = path_data;
%                     end
%                 end
      [filename_p, pathname] = uigetfile( ...
      {'*.mat', 'Mat files';}, ...
      'Pick a file',pathfiles);
      if filename_p==0
         handles.priznak=0;
         guidata(hObject, handles);
         helpdlg('No file selected.','Info');
      else
         path=[pathname,filename_p];
         load(path);
         set(handles.text3,'String', 'Frame'); 
         try correctfile=info.codewords;catch end;
         if (exist ('correctfile'))
         for ii=1:length(info.codewords) has_par=strcmp('par',info.codewords(ii)); if has_par==1 break; end;end;
         end;
         if (exist ('correctfile')) & has_par==1             
           parameters_tissue=info.tissue.descr;
           data{1}.info=info;
           var.mapy_parametru={};
           var.nazvy_parametru={};
           [curves_name,parameters_vizualization,parameters_tissue,fullmask,curves,concentration] = GetParamMaps(tissue,parameters_tissue);
           var.mapy_parametru=parameters_vizualization;
           var.nazvy_parametru=parameters_tissue;
           handles.curves=curves;
           data{1}.tissue=tissue;
           handles.concentration=concentration;
           handles.curves_name=curves_name;
           handles.mask=fullmask;
           guidata(hObject, handles);
           %This try and if is new - because of new file format - one picture in
           %data1 and data2 after perfusion analysis
           try newfilestructurewithonedata1picture=data2;catch end;
           if (exist ('newfilestructurewithonedata1picture'))
%                data{1}.data1=data1;
               data{1}.data2=data2;
               dohromady=num2str(length(data2));
               dohromady=['Frame <1-' dohromady '>'];
               set(handles.text3,'String', dohromady);
               handles.priznak=1;
               guidata(hObject, handles);
               helpdlg('Loading complete.','Info');
               set(handles.text1,'String','Loaded');
               path_data=pathname;
           else
               %This is old
           for i=1:length(info.source)
              waytoload=[pathname info.source{1,i}];
              try data2=load(info.source{1,i}, 'data2');catch end;
              try data2=load(waytoload, 'data2');catch end;
              try researchifcorrect=data2.data2;catch end;
              if (exist('researchifcorrect'))
              else
                  clear data2;
              end;
              if (exist('data2'))              
                 data2=data2.data2; 
%                  try data{1}.data1=load(info.source{1,i}, 'data1');catch end;
%                  try data{1}.data1=load(waytoload, 'data1');catch end;
%                  data{1}.data1=data{1}.data1.data1;                                          
                 dohromady=num2str(length(data2));
                 dohromady=['Frame <1-' dohromady '>'];
                 set(handles.text3,'String', dohromady); 
                 data{1}.data2=data2;
                 guidata(hObject, handles);
                 helpdlg('Loading complete.','Info');
                 set(handles.text1,'String','Loaded');
                 handles.priznak=1;
                 guidata(hObject, handles);
                 path_data=pathname;
                 break;              
               end;
            end;%end of old
            end;%end of new
            try data2=data2; catch end;
            if (exist('data2'))             
              return;
            else  
              backwards_add=pwd;
              cd(pathname);
              [filename_p, pathname] = uigetfile( ...
              {'*.mat', 'Mat files';}, ...
              'Choose file with first or second harmonic in mat format.',pathfiles);
              cd(backwards_add);
              if filename_p==0
                 handles.priznak=0;
                 guidata(hObject, handles);
                 clear global data;               
                 helpdlg('No file selected.','Info');
              else
                 path=[pathname,filename_p];
                 load(path, 'data2');              
                 try data2=data2; catch end;
                 if (exist('data2'))
%                     load(path, 'data1');
%                     data{1}.data1=data1;    
                    dohromady=num2str(length(data2));
                    dohromady=['Frame <1-' dohromady '>'];
                    set(handles.text3,'String', dohromady); 
                    data{1}.data2=data2;
                    handles.priznak=1;
                    guidata(hObject, handles);                
                    helpdlg('Loading complete.','Info');
                    set(handles.text1,'String','Loaded');
                 else
                    handles.priznak=0;
                    guidata(hObject, handles);
                    clear global data;
                    helpdlg('Wrong file selected.','Info');
                 end;
              end;
            end;
         else 
            handles.priznak=0;
            guidata(hObject, handles);
            clear global data;
            helpdlg('Wrong file selected.','Info');
         end;      
      end;
  case 'No'
      parameters_tissue=data{end}.info.perfan.params;
      tissue=data{end}.tissue;
      data{end}.tissue=tissue;
      info=data{end}.info;
      data{end}.info=info;
      var.mapy_parametru={};
      var.nazvy_parametru={};
      [curves_name,parameters_vizualization,parameters_tissue,fullmask,curves,concentration] = GetParamMaps(tissue,parameters_tissue);
      var.mapy_parametru=parameters_vizualization;
      var.nazvy_parametru=parameters_tissue;
      handles.curves=curves;
      handles.concentration=concentration;
      handles.curves_name=curves_name;
      handles.mask=fullmask;
      guidata(hObject, handles);
      data2=data{end}.data2;
      data{end}.data2=data2;
%       data1=data{end}.data1;
%       data{end}.data1=data1;
      dohromady=num2str(length(data2));
      dohromady=['Frame <1-' dohromady '>'];
      set(handles.text3,'String', dohromady);
      handles.priznak=1;
      guidata(hObject, handles);
end;   
end   
 
function Showing_button_Callback(hObject, eventdata, handles)
%%
global data
global var
if handles.priznak==1
   frame=get(handles.edit1);
   frame=frame.String;
   try frame_num=str2num(frame); catch end;
   if (exist('frame_num')) & frame_num<=(length(data{1}.data2)) & frame_num>0
      parameters_vizualization=var.mapy_parametru;
      parameters_tissue=var.nazvy_parametru;
      mask=handles.mask;
      curves=handles.curves;
      curves_name=handles.curves_name;
      concentration=handles.concentration;
%       vyber_dat=get(handles.radiobutton_data1,'Value');
%       frame=data{1}.data2{1,frame_num};
      frame=mean(cat(3,data{1}.data2{1,:}),3);
%       if vyber_dat==1
%         frame=data{1}.data1{1,frame_num};
%       else
%         frame=data{1}.data2{1,frame_num};    
%       end;
      krivkovy_index=1;
      parametr_name_krivky=cell(1,1);
      
      nmr_param_not_curves_counter=0;
      for nmr_param=1:length(parameters_vizualization)
        if parameters_vizualization{nmr_param}==0
           nmr_param_not_curves_counter=nmr_param_not_curves_counter+1;
        end;        
      end;
    nmr_param_not_curves=length(parameters_vizualization)-nmr_param_not_curves_counter;
    var.nmr_param_not_curves=nmr_param_not_curves;
    var.nmr_param_not_curves_counter=nmr_param_not_curves_counter;
    were_closed=0;
      for pocet_parametru=1:length(parameters_vizualization)
        parametr=parameters_vizualization{pocet_parametru};
        parameter_name=parameters_tissue{pocet_parametru+3};
        if parametr==0            
        else
            try existjizobrazene=var.priznak; catch end;
            if (exist('existjizobrazene'))
                try priznakovani=var.priznak_before; catch end;
                if (exist('priznakovani'))
                    if were_closed==0
                      if var.priznak==var.priznak_before
                        were_closed=1;
                        for zaviraci_parametr=1:var.priznak_before
                          if var.zavirani{zaviraci_parametr}.hlavni==0
                          else
                             delete(var.zavirani{zaviraci_parametr}.hlavni);                
                          end;
                        end;
                        var.priznak=0;            
                      end;  
                    end;
                    
                 end;
                
                VIZ_OUTPUT(parametr,parameter_name, frame, mask,curves,curves_name,concentration);
                try existjifigykrivek=var.pocitadlo; catch end;
                if (exist('existjifigykrivek'))
                   if var.pocitadlo==nmr_param_not_curves_counter | var.pocitadlo==nmr_param_not_curves_counter+1
                       if var.pocitadlo==nmr_param_not_curves_counter
                         for zaviraci_parametr=1:nmr_param_not_curves_counter
                           if isempty(var.zavirani_krivek{zaviraci_parametr})
                           else
                              delete(var.zavirani_krivek{zaviraci_parametr}.krivkyfig);                
                           end;
                         end;
                       else
                          for zaviraci_parametr=1:nmr_param_not_curves_counter+1
                           if var.zavirani_krivek{zaviraci_parametr}.krivkyfig==0
                           else
                              delete(var.zavirani_krivek{zaviraci_parametr}.krivkyfig);                
                           end;
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
                else
                   var.pocitadlo=0;
                end;                                
            else             
                try existjifigykrivek=var.pocitadlo; catch end;
                if (exist('existjifigykrivek'))
                  if var.pocitadlo==nmr_param_not_curves_counter | var.pocitadlo==nmr_param_not_curves_counter+1
                    for zaviraci_parametr=1:nmr_param_not_curves_counter
                      if var.zavirani_krivek{zaviraci_parametr}.krivkyfig==0
                      else
                        delete(var.zavirani_krivek{zaviraci_parametr}.krivkyfig);                
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
                else
                    var.pocitadlo=0;
                end;                                                                                
                VIZ_OUTPUT(parametr,parameter_name, frame, mask,curves,curves_name,concentration); 
            end;
        end;
      end;
      var.priznak_before=var.priznak;
   else
      helpdlg('You have to enter a number in given range.','Info');
   end;
else
   helpdlg('Nothing to show.','Info');
end;
end

function edit1_Callback(hObject, eventdata, handles)
%%
end

function edit1_CreateFcn(hObject, eventdata, handles)
%%
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
%%
function hlavni_close(hObject,eventdata)
global var
% global path_data
try existjiparam_not_curve=var.nmr_param_not_curves_counter; catch end;
if (exist('existjiparam_not_curve'))
    nmr_param_not_curves_counter=var.nmr_param_not_curves_counter;
    nmr_param_not_curves=var.nmr_param_not_curves;
end;
        delete(gcf);
        try existjizobrazene=var.priznak; catch end;
        if (exist('existjizobrazene'))
                if var.priznak==nmr_param_not_curves
                    for zaviraci_parametr=1:nmr_param_not_curves
                       if var.zavirani{zaviraci_parametr}.hlavni==0
                       else
                          delete(var.zavirani{zaviraci_parametr}.hlavni);                
                       end;
                    end;
                end;
        end;
        try existjifigykrivek=var.pocitadlo; catch end;
        if (exist('existjifigykrivek'))
                if var.pocitadlo==nmr_param_not_curves_counter
                    for zaviraci_parametr=1:nmr_param_not_curves_counter
                       if isempty(var.zavirani_krivek{zaviraci_parametr})
                       else
                          delete(var.zavirani_krivek{zaviraci_parametr}.krivkyfig);                
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
        else
                var.pocitadlo=0;
        end;
        % Path to files start
%         path_data=var.pathname;
        % Path to files end
        clear global data;    
        clear global handles;
        clear global var;       % global variables delete, only global variable data is kept
        rmpath('../../dce-mri/visualization');
        rmpath('../../dce-mri/visualization/ginput_impixelinfo');
        main_mri
        rmpath('../visualization/');
end

function checkbox_data_Callback(hObject, eventdata, handles)
%%
end

function radiobutton_data1_Callback(hObject, eventdata, handles)
%%
end

function radiobutton_data2_Callback(hObject, eventdata, handles)
%%
end


% --- Executes on button press in whole_sequence.
function whole_sequence_Callback(hObject, eventdata, handles)
global data
global path_data;
if handles.priznak==1
         for i=1:length(data{1}.info.source)
              waytoload=[path_data data{1}.info.source{1,i}];
              try data2=load(data{1}.info.source{1,i}, 'data2');catch end;
              try data2=load(waytoload, 'data2');catch end;
              try researchifcorrect=data2.data2;catch end;
              if (exist('researchifcorrect'))
              else
                  clear data1;
              end;
              if (exist('data2'))              
                 data2=data2.data2; 
%                  try data{1}.data1=load(info.source{1,i}, 'data1');catch end;
%                  try data{1}.data1=load(waytoload, 'data1');catch end;
%                  data{1}.data1=data{1}.data1.data1;                                          
                 dohromady=num2str(length(data2));
                 dohromady=['Frame <1-' dohromady '>'];
                 set(handles.text3,'String', dohromady); 
                 data{1}.data2=data2;
                 guidata(hObject, handles);
                 helpdlg('Loading complete.','Info');
                 set(handles.text1,'String','Loaded');
                 
                 guidata(hObject, handles);
                 
                 break;              
               end;
          end;  
       
          try data2=data2; catch end;
            if (exist('data2'))             
              return;
            else  
              backwards_add=pwd;
              cd(path_data);
              [filename_p, pathname2] = uigetfile( ...
              {'*.mat', 'Mat files';}, ...
              'Choose file with first or second harmonic in mat format.',path_data);
              cd(backwards_add);
              if filename_p==0                                                
                 helpdlg('No file selected.','Info');
              else
                 path=[pathname2,filename_p];
                 load(path, 'data2');              
                 try data2=data2; catch end;
                 if (exist('data2'))
%                     load(path, 'data1');
%                     data{1}.data1=data1;    
                    dohromady=num2str(length(data2));
                    dohromady=['Frame <1-' dohromady '>'];
                    set(handles.text3,'String', dohromady); 
                    data{1}.data2=data2;
                    
                    guidata(hObject, handles);                
                    helpdlg('Loading complete.','Info');
                    set(handles.text1,'String','Loaded');
                 else
                    
                    
                    helpdlg('Wrong file selected.','Info');
                 end;
              end;
            end;
else
   helpdlg('No data loaded.','Info');
end;




end
