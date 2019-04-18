classdef xjViewHandler < handle
    %XJVIEWHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xjViewFigHandle;%the figure handle of an xjView
        maskName;%name of T or F map file
    end
    
    methods
        function obj=xjViewHandler(mask,xjViewFigHandle,selection,minSize)
            if(~exist('selection','var'))
                selection='';
            end
            if(~exist('minSize','var'))
                minSize=0;
            end
            if(~exist('xjViewFigHandle', 'var'))
                xjViewFigHandle=[];
            end
            if(~exist('mask','var'))
                mask=[];
            end
            obj.maskName=mask;
            if(isempty(xjViewFigHandle))
                newInstance=true;
                xjViewFigHandle=xjViewHandler.get_xjView(newInstance);
            end
            guiData=guidata(xjViewFigHandle);
            if(~isempty(mask))
                if(~xjViewHandler.IsOpen(xjViewFigHandle,mask))
                    guiData.loadImagePush.Callback(xjViewFigHandle,[],mask);
                    if(~isempty(selection))
                        if(strcmp(selection,'+')||strcmp(selection,'-'))
                            guiData.allIntensityRadio.Callback(guiData.allIntensityRadio,[],selection);
                        end
                    end
                    if(~isempty(minSize))
                        set(guiData.clusterSizeThresholdEdit,'String',num2str(minSize));
                        guiData.clusterSizeThresholdEdit.Callback(guiData.clusterSizeThresholdEdit,[]);
                    end
                    obj.maskName=mask;
                    spm_orthviews('reposition');
                end
            end
            obj.xjViewFigHandle=xjViewFigHandle;
        end
        function selectPositiveOnly(obj)
            handles=guidata(obj.xjViewFigHandle);
            handles.allIntensityRadio.Callback(handles.allIntensityRadio,[],'+');
        end
        function h=getXjviewHandles(obj)
            h=guidata(obj.xjViewFigHandle);
        end
        function h=getXjviewFigHandle(obj)
            h=obj.xjViewFigHandle;
        end
    end
    methods (Static = true)
        function h=get_xjView(newInstance)
            if(~exist('newInstance','var'))
                newInstance=false;
            end
            if(newInstance)
                h=xjview;
            else
                h=[];
                hs=get(0, 'children');
                for i=1:length(hs)
                    name=hs(i).Name;
                    if(length(name)<6)
                        continue;
                    end
                    if(strcmp(name(1:6), 'xjView'))
                        h=hs(i);
                        break
                    end
                end
                if(isempty(h))
                    return;
                end
                FM.addFigure(h,'xjView');
            end
            p=getpixelposition(h);
            w=1500;
            ratio=w/p(3);
            p=get(h,'Position');
            p=p.*[1 1 ratio ratio];
            set(h,'position', p);
            figure(h);
        end
        function IS=IsOpen(xjViewFigHandle,img)
            %testing whether img is already opened by the xjView
            IS=false;
            guiData=guidata(xjViewFigHandle);
            if(~isfield(guiData,'imageFileName'))
                return;
            end
            IS=strcmp(img,guiData.imageFileName{1});
        end
        function IS=isvalidhandle(h)
            IS=false;
            if(isempty(h))
                return;
            end
            if(~ishandle(h))
                return;
            end
            if(~isvalid(h))
                return;
            end
            name=get(h,'Name');
            if(length(name)<6)
                return;
            end
            IS=strcmp(name(1:6),'xjView');
        end
    end    
end

