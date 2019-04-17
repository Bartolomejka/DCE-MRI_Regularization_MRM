function reactToImageChangesInFig(target_image,h_caller,deleteFcn,refreshFcn)
%reactToImageChangesInFig sets up listeners to react to image changes.
%   reactToImageChangesInFig(TARGET_IMAGE,H_CALLER,DELETE_FCN) calls
%   DELETE_FCN if TARGET_IMAGE is deleted or if the CData property of the
%   target image TARGET_IMAGE changes.  DELETE_FCN is a function handle
%   specified by the caller that deletes the tool, H_CALLER. TARGET_IMAGE
%   can be array of handles to graphics image objects.
%
%      DELETE_FCN is called when:
%      ==========================
%      * TARGET_IMAGE is deleted
%      * the CData property of TARGET_IMAGE is modified
%
%   reactToImageChangesInFig(TARGET_IMAGE,H_CALLER,DELETE_FCN,REFRESH_FCN)
%   calls only the DELETE_FCN if TARGET_IMAGE is deleted and calls the
%   REFRESH_FCN if the CData property of TARGET_IMAGE is modified.  The
%   REFRESH_FCN is a function handle specified by the caller that would
%   cause H_CALLER or to refresh itself.
%
%      DELETE_FCN is called when:
%      ==========================
%      * TARGET_IMAGE is deleted
%
%      REFRESH_FCN is called when:
%      ===========================
%      * the CData property of TARGET_IMAGE is modified
%
%   See also IMPIXELINFO,IMPIXELINFOVAL.

%   Copyright 2004-2008 The MathWorks, Inc.
%   $Revision: 1.1.8.7 $  $Date: 2008/12/22 23:48:03 $

% validate input
iptchecknargin(3,4,nargin,mfilename);
checkImageHandleArray(target_image,mfilename);

% refreshFcn function is optional
if nargin == 4
    refresh_fcn_provided = true;
else
    refresh_fcn_provided = false;
end

% call the deleteFcn if image is destroyed
objectDestroyedListener = iptui.iptaddlistener(target_image,...
    'ObjectBeingDestroyed',deleteFcn);
storeListener(h_caller,objectDestroyedListener);

% call appropriate function if image cdata changes
if refresh_fcn_provided
    imageCDataChangedListener = iptui.iptaddlistener(target_image,...
        'CData','PostSet',refreshFcn);
else
    imageCDataChangedListener = iptui.iptaddlistener(target_image,...
        'CData','PostSet',deleteFcn);
end
storeListener(h_caller,imageCDataChangedListener);

% create and install the listener manager the modular tool.
if refresh_fcn_provided
%     listener_manager = iptui.listenerManager(h_caller,refreshFcn);
    listener_manager = iptui.iptaddlistener(h_caller,'ObjectChanged',refreshFcn);

else
    listener_manager = iptui.listenerManager(h_caller,[]);
end
listener_manager.installManager(h_caller);


%-----------------------------------------------------
function storeListener(h_caller,listener)
% this function stores the listeners in the appdata of the h_caller object
% using the makelist function.

% grab current list of listeners from the caller's appdata
listenerList = getappdata(h_caller,'imageChangeListeners');
if isempty(listenerList)
    listenerList = makeList;
end
listenerList.appendItem(listener);
setappdata(h_caller,'imageChangeListeners',listenerList);

