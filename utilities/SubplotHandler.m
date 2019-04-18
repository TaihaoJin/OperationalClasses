classdef SubplotHandler < handle
    %SUBPLOTHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        QPGuiHandler;
        hAxs;
        hFig;
        position;
        ButtonDownFcn;
        rows;
        cols;
        rank;%rank = 1, 2, 3, .... rows*cols
        lx;
        ly;
    end
    
    methods
        function obj = SubplotHandler()
            obj.QPGuiHandler=[];
            obj.hAxs=[];
            obj.hFig=[];
            obj.position=[];
            obj.ButtonDownFcn=[];
            obj.rows=[];
            obj.cols=[];
            obj.rank=[];
            obj.lx=[];
            obj.ly=[];
        end
        
        function sub = clone(obj)
            sub=SubplotHandler();
            sub.QPGuiHandler=obj.QPGuiHandler;
            sub.hAxs=obj.hAxs;
            sub.hFig=obj.hFig;
            sub.position=obj.position;
            sub.ButtonDownFcn=obj.ButtonDownFcn;
            sub.rows=obj.rows;
            sub.cols=obj.rows;
            sub.rank=obj.rank;
        end
        
        function updateQPGuiHandler(obj, hObject)
            if(isempty(obj.QPGuiHandler))
                obj.QPGuiHandler=quickPlotterGuiDataHandler();
            end
            obj.QPGuiHandler.update(hObject);
        end
        
        function updateData(obj, hObject)
            obj.QPGuiHandler.updateData(hObject);
        end
        
        function updatePlottingOption(obj, hObject)
            obj.QPGuiHandler.updatePlottingOption(hObject);
        end
        
        function updatePlottingType(obj, hObject)
            obj.QPGuiHandler.updatePlottingType(hObject);
        end
        
        function plot(obj)
            obj.QPGuiHandler.plot(obj)
        end
        
        function plot0(obj, xData, yData, xInds, yInds)   
            
            CommonMethods.purgeOverlappingAxes(obj.hFig, obj.hAxs)
            figure(obj.hFig);
            if(isempty(obj.hAxs)||~ishandle(obj.hAxs))
                obj.hAxs=axes('position', obj.position, 'ButtonDownFcn', obj.ButtonDownFcn);                
            end
            
            axes(obj.hAxs);
            
            XL=obj.xLabel;
            YL=obj.yLabel;
            if(~exist('xData','var'))
                xData=obj.xData;
                yData=obj.yData;
            end
            if(~exist('xInds', 'var'))
                xInds=obj.xInds;
                yInds=obj.yInds;
            end            
            
            Style=obj.plottingStyle;
            option=obj.plottingOption;
            if(strcmpi(Style, 'Line'))
                X=xData(:,xInds);
                Y=yData(:,yInds);
                plot(X,Y);
            elseif(strcmpi(Style, 'Scatter'))
                X=xData(:,xInds);
                Y=yData(:,yInds);
                scatter(X,Y);
            elseif(strcmpi(Style, 'Histogram'))
                Y=yData(:,yInds);
                hist(Y);
                XL=YL;
                YL='Counts';
                if(strcmp(option,'Y Hist_Group'))
                    YL='Counts\_Group';
                end
            elseif(strcmpi(Style, 'Plotyy'))
                X1=xData(:,xInds{1});
                Y1=yData(:,yInds{1});
                X2=xData(:,xInds{2});
                Y2=yData(:,yInds{2});       
                [Ax,~,~]=plotyy(X1,Y1,X2,Y2);
                set(Ax(1), 'ButtonDownFcn', obj.ButtonDownFcn);
                set(Ax(2), 'ButtonDownFcn', obj.ButtonDownFcn);
                ylabel(Ax(2),YL{2},'fontsize', 14);
                YL=YL{1};
                obj.hAxs=Ax;
            end
            xlabel(XL,'fontsize',14);
            ylabel(YL,'fontsize',14); 
            set(gca, 'ButtonDownFcn', obj.ButtonDownFcn);
            obj.hAxs=gca;
        end
        
        function des = getDescription(obj)
            row=floor((obj.rank-1)/obj.cols)+1;
            col=mod((obj.rank-1),obj.cols)+1;
            des=[num2str(row) 'r_' num2str(col) 'c'];
        end
        
        function markCurrentPoint(obj,cp)
            axes(obj.hAxs(1));
            
            if(~isempty(obj.lx))
                try
                    delete(obj.lx);
                catch 
                end
            end
            
            if(~isempty(obj.ly))
                try
                    delete(obj.ly);
                catch
                end
            end
            
            r=0.05;
            xr0=xlim;
            xr=r*(xr0(2)-xr0(1));
            
            yr0=ylim;
            yr=r*(yr0(2)-yr0(1));
            
            x=cp(1,1);
            y=cp(1,2);
            
            xr=[max(x-xr,xr0(1));min(x+xr,xr0(2))];
            xr=xlim;
%            yr=[max(y-yr,yr0(1));min(y+yr,yr0(2))];
            yr=ylim;
            
            obj.lx=line(xr,[y;y],'color', 'green');
            obj.ly=line([x;x],yr,'color', 'green');
            
            xlim(xr0);
            ylim(yr0);
        end
        
        function row = getCurrentRow(obj, cp)
            row=[];
            if(~isempty(obj.QPGuiHandler))
                row=obj.QPGuiHandler.getCurrentRow(cp);
            end
        end
    end
    
end

