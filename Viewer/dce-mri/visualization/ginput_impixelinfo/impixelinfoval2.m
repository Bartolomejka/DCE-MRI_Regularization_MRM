function huic = impixelinfoval2(parent,himage)


iptcheckhandle(parent,{'figure','uipanel','uicontainer'},mfilename, ...
   'HPARENT',1);
% checkImageHandleArray(himage,mfilename);

[himage,axesHandles,hFig] = imhandles(himage);

huic = uicontrol('Parent',parent,...
   'Style','text',...
   'Units','pixels',...
   'HorizontalAlignment','left',...
   'BusyAction','queue',...
   'enable', 'inactive',...
   'Interruptible','off',...
   'Visible','on',...
   'DeleteFcn',@deleteFcn,...
   'Tag','pixelinfo text');

setChildColorToMatchParent(huic,parent);

% declare so they will have scope in nested functions
imageModels = getimagemodel(himage);
xystring = '(X, Y)';
twoSpaces = '  ';
oneSpace = ' ';
stringForMultipleImages = 'Pixel Value';

is_rset = strcmp(get(himage(1),'tag'),'rset overview');

% create Context menu
cmenu = uicontextmenu('Parent', hFig);
set(himage, 'UIContextMenu', cmenu);         
menuItemHandle = uimenu(cmenu,...
   'Label','Copy pixel info', ...
   'Callback',@copyPixInfoToClipboard);

   %---------------------------------------
   function copyPixInfoToClipboard(obj,evt) %#ok<INUSD>
      clipboard('copy', pixelString);
   end

callbackID = iptaddcallback2(hFig,'WindowButtonMotionFcn', @displayPixelInfo);
reactToImageChangesInFig(himage,huic,@reactDeleteFcn,@reactRefreshFcn);
registerModularToolWithManager(huic,himage);


   %-------------------------------
   function reactDeleteFcn(obj,evt) %#ok<INUSD>
      if ishghandle(huic)
         delete(huic);
      end
   end


   %-------------------------------
   function reactRefreshFcn(obj,evt) %#ok<INUSD>
       imageModels = getimagemodel(himage);
   end


   %--------------------------
   function deleteFcn(obj,evt) %#ok<INUSD>

      iptremovecallback(hFig,'WindowButtonMotionFcn',callbackID);
      if ishghandle(cmenu)
         delete(cmenu);
      end
      if ishghandle(menuItemHandle)
         delete(menuItemHandle);
      end

   end


   %---------------------------------
   function displayPixelInfo(obj,evt) %#ok<INUSD>
       global var
       maska=var.mask_length;
       [index,x,y] = findAxesThatTheCursorIsOver(axesHandles,maska);
       him = himage(index);

       if isempty(him) || isCursorOutsideVisibleImage(him,x,y)
           % if we are not over an image...
           displayDefaultString;
       
       elseif strcmp(get(him,'tag'),'rset overview')
           % if we are over an rset overview image...
           createXString = @(x) sprintf('%1d', round(x));
           createYString = createXString;
           pixelString = sprintf('(%s,%s%s)', ...
               createXString(x), ...
               oneSpace, ...
               createYString(y));
           set(huic,'String',pixelString);
           
       else
           % else display pixel info string as usual
           [nrows,ncols,ndims] = size(get(him,'Cdata')); %#ok<NASGU>
           if nrows == 0 || ncols == 0
               displayDefaultString;
             
           else
               [isXDataDefault,isYDataDefault] = isDefaultXYData(him);
               xdata = get(him,'XData');
               ydata = get(him,'YData');
               
               if isXDataDefault && isYDataDefault
                   createXString = @(x) sprintf('%1d', round(x));
                   createYString = createXString;

               else
                   createXString = @(x) sprintf('%1.2f', x);
                   createYString = createXString;

                   isXDataTheSame = xdata(2) == xdata(1);
                   isYDataTheSame = ydata(2) == ydata(1);
                   if isXDataTheSame
                     createXString = @(x) sprintf('%1.2f', xdata(1));
                   end
                   
                   if isYDataTheSame
                     createYString = @(y) sprintf('%1.2f', ydata(1));
                   end
               end

               % Construct pixel value string
               imModel = imageModels(index);
               rp = axes2pix(nrows,ydata,y);
               cp = axes2pix(ncols,xdata,x);
               r = min(nrows,max(1,round(rp)));
               c = min(ncols,max(1,round(cp)));

%                valueString = getPixelInfoString(imModel,r,c);

%                 xova=round(x);
%                 yova=round(y)+1;
                xova=c;
                yova=r;
                
                current_bild= parent;
                my_name2=get(current_bild,'Parent');
                my_name=get(my_name2,'Name');
                maps_param=var.nazvy_parametru;
                current_name=0;
                for pointer=1:length(maps_param)-3
                    issame=strcmp(maps_param{pointer+3},my_name);
                    if issame==0
                        current_name=current_name+1;
                    else
                        current_name=current_name+1;
                        break
                    end;
                end;
                
               valueString=num2str(var.mapy_parametru{current_name}(yova,xova));
               
               valueString = fixENotation(valueString);

               pixelString = sprintf('(%s,%s%s)%s%s', ...
                                     createXString(x), ...
                                     oneSpace, ...
                                     createYString(y), ...
                                     twoSpaces, ...
                                     valueString);
               set(huic,'String',pixelString);
           end
       end

       %----------------------------
       function displayDefaultString

           if numel(himage) ~= 1
               defaultPixelInfoString = stringForMultipleImages;
           elseif is_rset
               defaultPixelInfoString = '';
           else
               defaultPixelInfoString = getDefaultPixelInfoString(imageModels);
           end
           pixelString = sprintf('%s%s%s',xystring,twoSpaces, ...
               defaultPixelInfoString);
           set(huic,'String',pixelString);
       end

   end % displayPixelInfo

   % Set position of string so that it could fit the largest pixel value.
   % Turn pixel reporting on by default, after resetting hPixelInfoValue's
   % string to ''.
   set(huic,'String','(0000,0000) [0.00E+00 0.00E+00 0.00E+00]');
   pixInfoExtent = get(huic,'Extent');
   fudge = 2;
   posPixInfoValue = [1 1 pixInfoExtent(3)+3*fudge pixInfoExtent(4)];
   set(huic,'Position',posPixInfoValue);
   pixelString = sprintf('%s%s%s',xystring,twoSpaces, stringForMultipleImages);
   set(huic,'String',pixelString);
   set(huic,'Visible','on');

end % impixelinfoval


%--------------------------------------------------------
function isOutside = isCursorOutsideVisibleImage(him,x,y)

hSp = imshared.getimscrollpanel(him);

if isempty(hSp)
    isOutside = false;

else
    api = iptgetapi(hSp);
    rect = api.getVisibleImageRect();
    xmin = rect(1);
    ymin = rect(2);
    width = rect(3);
    height = rect(4);
    outsideRect = x < xmin || x > xmin + width || y < ymin || y > ymin ...
        + height;
    if outsideRect
        isOutside = true;
    else
        isOutside = false;
    end
end

end % isCursorOutsideVisibleImage
