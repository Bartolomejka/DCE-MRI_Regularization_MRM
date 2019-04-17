function hpanel = impixelinfo2(varargin)


[h,parent] = parseInputs(varargin{:});

if strcmp(get(parent,'Type'),'figure')
    parentIsFigure = true;
else
    parentIsFigure = false;
end

imageHandles = imhandles(h);
hFig = ancestor(h,'Figure');

% hFig will be a cell array if imageHandles is an array, even though
% imageHandles may belong to the same figure.
if iscell(hFig) 
  hFig = hFig{1};
end

if isempty(imageHandles)
    eid = sprintf('Images:%s:noImageInFigure',mfilename);
    msg = 'The figure must contain at least one image.';
    error(eid,'%s',msg);
end

hPixInfoPanel = createPanel;

reactToImageChangesInFig(imageHandles,hPixInfoPanel,@reactDeleteFcn,[]);
registerModularToolWithManager(hPixInfoPanel,imageHandles);

if isequal(parent,ancestor(imageHandles,'figure')) && ...
        strcmp(get(parent,'Visible'),'on')
    figure(parent);
end

if nargout > 0
    hpanel = hPixInfoPanel;
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function reactDeleteFcn(obj,evt) %#ok<INUSD>
        if ishghandle(hPixInfoPanel)
            delete(hPixInfoPanel);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function hPixInfoPanel = createPanel

        units = 'Pixels';
        posPanel = [1 1 300 20];
        visibility = 'on';

        if parentIsFigure
            backgrndColor = get(parent,'Color');
        else
            backgrndColor = get(parent,'BackgroundColor');
        end

        fudge = 2;

        hPixInfoPanel = uipanel('Parent',parent,...
            'Units',units,... 
            'Position',posPanel,...
            'Tag','pixelinfo panel',...
            'Visible',visibility,...
            'Bordertype','none',...
            'BackgroundColor', backgrndColor);

        hPixelInfoLabel = uicontrol('Parent',hPixInfoPanel,...
            'Style','text',...
            'String','Pixel info:', ...
            'Tag','pixelinfo label',...
            'Units','pixels',...
            'Visible',visibility,...
            'BackgroundColor',backgrndColor);

        labelExtent = get(hPixelInfoLabel,'Extent');
        posLabel = [posPanel(1) posPanel(2) labelExtent(3) labelExtent(4)];
        set(hPixelInfoLabel,'Position',posLabel);

        % initialize uicontrol that will contain the pixel info values.
        hPixelInfoValue = impixelinfoval2(hPixInfoPanel,imageHandles);
        posPixInfoValue = get(hPixelInfoValue,'Position');
        set(hPixelInfoValue,'Position',[posLabel(1)+posLabel(3) posPanel(2) ...
            posPixInfoValue(3) posPixInfoValue(4)]);
        posPixInfoValue = get(hPixelInfoValue,'Position');

        % link visibility of hPixInfoPanel and its children 
        hlink = linkprop([hPixInfoPanel hPixelInfoLabel hPixelInfoValue],...
            'Visible');
        setappdata(hPixInfoPanel,'linkToChildren',hlink);

        newPanelWidth = posPixInfoValue(1)+posPixInfoValue(3)+fudge;
        newPanelHeight = max([posLabel(4) posPixInfoValue(4)]) + 2*fudge;
        set(hPixInfoPanel,'Position', [posPanel(1) posPanel(2) newPanelWidth ...
            newPanelHeight]);

        set(hPixInfoPanel,'Visible','on');  % must be in this function.

    end

end  %main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,parent] = parseInputs(varargin)

iptchecknargin(0,2,nargin,mfilename);

switch nargin
    case 0
        %IMPIXELINFO
        h = get(0, 'CurrentFigure');
        if isempty(h)
            eid = sprintf('Images:%s:noImageInFigure',mfilename);
            msg = sprintf('%s expects a current figure containing an image.', ...
                upper(mfilename));
            error(eid,'%s',msg);
        end
        parent = h;

    case 1
        h = varargin{1};
        if ~ishghandle(h)
            eid = sprintf('Images:%s:invalidGraphicsHandle',mfilename);
            msg = 'H is an invalid graphics handle.';
            error(eid,'%s',msg);
        end
        parent = ancestor(h,'Figure');

    case 2
        parent = varargin{1};
        if ishghandle(parent)
            type = get(parent,'type');
            if ~strcmp(type,'uipanel') && ~strcmp(type,'uicontainer') && ...
                    ~strcmp(type,'figure')
                eid = sprintf('Images:%s:invalidParent',mfilename);
                msg = 'HPARENT must be a handle to a UIPANEL or a FIGURE.';
                error(eid,'%s%s',msg);
            end
        else
            eid = sprintf('Images:%s:invalidGraphicsHandle',mfilename);
            msg = 'HPARENT is an invalid graphics handle.';
            error(eid,'%s',msg);
        end

        h = varargin{2};
        checkImageHandleArray(h,mfilename);
end

end
