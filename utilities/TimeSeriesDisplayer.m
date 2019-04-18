classdef TimeSeriesDisplayer
    %TIMESERIESDISPLAYER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        figTitle;
        titles;%titles of subplots
        files;
        tSeries;
        tSeries_R;%time series of processed
        yy=[];
        xx;
        xRange;%range of x values
        noVols;
        num;%no of 4D images
        rows;%row and cols is the arrangement of subplots
        cols;
        dims;
        noSess;
        fh=[-1 -1 -1];%figure hadle
        Axes;%Axes of subplots
        Lines;%Line objects of subplots
        Betas;
    end
    
    properties (Dependent = true)
        Residuals;
    end
        
    methods
        function obj = TimeSeriesDisplayer(files,yy,statsDir,conName)
%            obj.figTitle=figTitle;
            obj.files=files;
            obj.yy=yy;       
            len=length(obj.files);
            obj.num=len;
            ts=cell(1,len);
            xx=cell(1,len);
            dims=cell(1,len);
            xRanges=zeros(len,2);
            for i=1:len
                cprintf('blue','%s\n',['reading files: ' obj.files{i}]);
                ts{i}=spm_read_vols(spm_vol(obj.files{i}));
                titles{i}=CommonMethods.getFileName(obj.files{i});              
                xx{i}=(1:size(ts{i},4))';
                xRanges(i,1)=min(xx{i});
                xRanges(i,2)=max(xx{i});
                dims{i}=size(ts{i});
            end
            obj.xx=xx;
            obj.fh(1) = figure();
            cols=floor(sqrt(len));
            rows=cols;
            if(cols*rows<len)
                cols=cols+1;
            end
            if(cols*rows<len)
                rows=rows+1;
            end
            obj.rows=rows;
            obj.cols=cols;
            obj.tSeries=ts;
            obj.titles=titles;
            obj.tSeries_R=ts;
            obj.dims=dims;
            obj.noVols=cellfun(@(x) x(4),dims);
            obj.Axes=zeros(1,obj.num);
            obj.Lines=zeros(1,obj.num);
            obj.xRange=[min(xRanges(:,1)) max(xRanges(:,2))];
            
            if(exist('statsDir','var'))               
                %build corrected time series
                load(fullfile(statsDir,'SPM.mat'));
                mask=fullfile(statsDir,'mask.img');
                maskVol=spm_read_vols(spm_vol(mask));
                
                %the design matrix
                desMats=SPM.xX.X;
                %column names of the design matrix
                desNames=SPM.xX.name;
                noSess=length(SPM.Sess);
                obj.noSess=noSess;
                noVols=obj.noVols(1);
                dim=dims{1};
                dimt=[dim(4) dim(1)*dim(2)*dim(3)];
                noPars=length(desNames);
                Betas=arrayfun(@(p)spm_read_vols(spm_vol(fullfile(statsDir,['beta_' sprintf('%04d',p) '.img']))), 1:noPars,'UniformOutput',false);                
                for s=1:noSess
                    ts=obj.tSeries{s};
                    pos=CommonMethods.findstr_StrInCellArray(desNames,['Sn(' num2str(s) ')']);
                    columns=find(cellfun(@(x) ~isempty(x), pos));
                    parNames=desNames(columns);
                    desMat=desMats((noVols*(s-1)+1):s*noVols, columns);
                    betas=arrayfun(@(x)Betas{x},columns,'UniformOutput',false);
                    %exacted desMat and beta of the run s
                    
                    conds=CommonMethods.strsplit_CellArray({'+' '-'}, conName);     
                    pos=cellfun(@(x) CommonMethods.findstr_CellArrayInStr(x,conds), parNames, 'UniformOutput', false);
                    condColumns=find(cellfun(@(x) ~isemptyr(x), pos));
                    obj.yy{s}=desMat(:,condColumns);
                    
                    regressingColumns=find(~ismember(1:length(columns), condColumns));
                    %regressingColumns: the column numbers of the
                    %regressor. Numbering is within the columns in desMat
                    desMat=desMat(:,regressingColumns);
                    betas=betas(regressingColumns);            
                    %excluded the clolumns containing the condition name
                    %from the design matrix and the betas
                    betasr=betas{1};
                    for i=2:length(betas)
                        betasr=cat(4,betasr,betas{i});
                    end
                    betasr=reshape(betasr,[dimt(2) size(betasr,4)])';
                    
                    tsr=ts;
                    tsr=reshape(tsr, [dimt(2) dimt(1)])';
                    residuals=tsr-desMat*betasr;     
                    obj.tSeries_R{s}=reshape(residuals',dim);
                end
             obj.fh(2) = figure();
             obj.fh(3) = figure();
           end
        end
        
        function display(obj,v)
            [ranges, rangesr,rangesyy] = obj.getTSRanges(v);
            rmy=max(cellfun(@(x) x.range, ranges));
            rmyr=max(cellfun(@(x) x.range, rangesr));
            rmyy=max(cellfun(@(x) x.range, rangesyy));
            
            figure(obj.fh(1));            
            for i=1:obj.num
                y=squeeze(obj.tSeries{i}(v(1),v(2),v(3),:));
                yy=obj.yy{i};
                subplot(obj.rows,obj.cols,i);
                title(obj.titles{i});
                ylimt=CommonMethods.getDispRange(y,rmy);
                if(isempty(yy))
                    plot(obj.xx{i},y);
%                    set(gca,'ylim',ylimt,'ytickmode', 'auto');
                else
                    yylimt=CommonMethods.getDispRange(yy,rmyy);                    
                    ax=plotyy(obj.xx{i},y,obj.xx{i},yy);
%                    set(ax(1),'ylim',ylimt,'ytickmode', 'auto');
%                    set(ax(2),'ylim',yylimt,'ytickmode', 'auto');
                end
            end
            
            figure(obj.fh(2));            
            for i=1:obj.num
                y=squeeze(obj.tSeries_R{i}(v(1),v(2),v(3),:));
                yy=obj.yy{i};
                subplot(obj.rows,obj.cols,i);
                title(obj.titles{i});
                ylimt=CommonMethods.getDispRange(y,rmyr);
                if(isempty(yy))
                    plot(obj.xx{i},y);
 %                   set(gca,'ylim',ylimt,'ytickmode', 'auto');
                else
                    yylimt=CommonMethods.getDispRange(yy,rmyy);
                    ax=plotyy(obj.xx{i},y,obj.xx{i},yy);
 %                   set(ax(1),'ylim',ylimt,'ytickmode', 'auto');
%                    set(ax(2),'ylim',yylimt,'ytickmode', 'auto');
                end
            end
            
            figure(obj.fh(3));            
            for i=1:obj.num
                y=squeeze(obj.tSeries{i}(v(1),v(2),v(3),:))-squeeze(obj.tSeries_R{i}(v(1),v(2),v(3),:));
                yy=obj.yy{i};
                subplot(obj.rows,obj.cols,i);
                title(obj.titles{i});
                ylimt=CommonMethods.getDispRange(y,rmyr);
                if(isempty(yy))
                    plot(obj.xx{i},y);
 %                   set(gca,'ylim',ylimt,'ytickmode', 'auto');
                else
                    yylimt=CommonMethods.getDispRange(yy,rmyy);
                    ax=plotyy(obj.xx{i},y,obj.xx{i},yy);
 %                   set(ax(1),'ylim',ylimt,'ytickmode', 'auto');
%                    set(ax(2),'ylim',yylimt,'ytickmode', 'auto');
                end
            end
            
        end
        
        function [ranges, rangesr, rangesyy] = getTSRanges(obj,v)
            %returns the data range at the voxel v
            ranges=cell(1,obj.num);
            rangesr=cell(1,obj.num);            
            rangesyy={};         
            
            for i=1:obj.num
                y=squeeze(obj.tSeries{i}(v(1),v(2),v(3),:));
                yr=squeeze(obj.tSeries_R{i}(v(1),v(2),v(3),:));
                ranges{i}=CommonMethods.getRange(y);
                rangesr{i}=CommonMethods.getRange(yr);
                yy=obj.yy{i};
                if(~isempty(yy))
                    rangesyy{end+1}=CommonMethods.getRange(obj.yy{i});
                end
            end
        end
    end
    
end

