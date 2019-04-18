classdef spmTabDataHandler < handle
    %SPMTABDATAHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    %This class is to handle the tabData of SPM volume search
    properties
        tabData
    end
    
    methods
        function obj=spmTabDataHandler(tabData, V)
            if(~exist(fullfile(CommonEnvMethods.getMatlabProjectPath(),'MRI', 'Scripts', 'Modified Scripts', 'saved_objects','emptySPMTabData.mat'),'file'))
                save(fullfile(CommonEnvMethods.getMatlabProjectPath(),'MRI', 'Scripts', 'Modified Scripts', 'saved_objects','emptySPMTabData.mat'), 'tabData');     
            end
            if(isempty(tabData))
                load(fullfile(projectPath,'MRI', 'Scripts', 'Modified Scripts', 'saved_objects','emptySPMTabData.mat'));
            end
            obj.tabData=tabData;
            obj.retrieveClusterInfo();
            tabData=obj.tabData;
            if(exist('V','var'))
                tabData.mat=V.mat;
            else
                obj.tabData.mat=[];
                return;
            end
            if(~isempty(tabData.dat))
            coors=cellfun(@(x)x',tabData.dat(:,11),'UniformOutput',false);
            coors=cell2mat(coors);
            volIndexes=CommonMethods.coor2ind_mat(coors,tabData.mat);
            volInd1=sub2ind(V.dim, volIndexes(:,1), volIndexes(:,2), volIndexes(:,3));
            else
                coors=[];
                volIndexes=[];
                volInd1=[];
            end
            tabData.volIndexes=volIndexes;
            tabData.volInd1=volInd1;
            tabData.coors=coors;
            tabData.vars=arrayfun(@(x)[tabData.hdr{2,x} '_' tabData.hdr{1,x}], 1:10, 'UniformOutput', false);
            tabData.vars{end+1}=tabData.hdr{2,11};
            tabData.vars{end+1}='clusterInd';
            obj.tabData=tabData;
        end
        function retrieveClusterInfo(obj)
            %assign cluster index for each local maximum
            tabData=obj.tabData;
            if(isempty(tabData.dat))
                tabData.p_cluster=[];
                tabData.clusterInds=[];
                tabData.clusterPeak=[];
            else
                p_cluster=cellfun(@(x)x,tabData.dat(:,3),'UniformOutput',false);
                inds=find(cellfun(@(x)~isempty(x),p_cluster));%index of the first voxel (cluster peak) of each cluster
                p_cluster=arrayfun(@(x)tabData.dat{x,3},inds);%cluster level p value of each cluster
                clusterSize=arrayfun(@(x)tabData.dat{x,4},inds);%cluster level p value of each cluster
                tabData.p_cluster=p_cluster;
                len=size(tabData.dat, 1);
                clusterInds=zeros(1,len);
                clusterPeak=false(1,len);
                clusterPeak(inds)=true;
                lent=length(inds);%number of clusters
                cols=size(tabData.dat,2);
                for i=1:lent
                    iI=inds(i);
                    if(i<lent)
                        iF=inds(i+1);
                    else
                        iF=len;
                    end
                    clusterInds(iI:iF)=i;
                    tabData.dat(iI:iF,3)={p_cluster(i)};
                    tabData.dat(iI:iF,4)={clusterSize(i)};
                    tabData.dat(iI:iF,cols+1)={i};
                end
                tabData.clusterInds=clusterInds;
                tabData.clusterPeak=clusterPeak;
            end
            obj.tabData=tabData;
        end
        function [vars, values]=getLocalMaximum(obj,ind)            
            vars=obj.tabData.vars;
            values=obj.tabData.dat(ind,:);
        end
        function clusterInd=getClusterInd(obj,volInds1)
            clusterInd=nan;
            %this function return the cluster index in tabData 
            %volInds1 is 1-d voxel indexes of a cluster
            if(~isfield(obj.tabData, 'volInd1'))
                return;
            end
            nnzs=arrayfun(@(x)nnz(volInds1==x),obj.tabData.volInd1);
            idx=find(nnzs);
            if(isempty(idx))
                return;
            end
            clusterInd=obj.tabData.clusterInds(idx(1));
        end
        function p=getClusterP(obj,ind)
            %ind is the cluster index
            p=obj.tabData.p_cluster(ind);
        end
        
        function table=getSummaryTable(obj, cutoff_p_cluster)
 
%     'set'    'set'    'cluster'    'cluster'    'cluster'    'voxel'         'voxel'         'voxel'    'voxel'     'voxel'     ''          
%     'c'      'p'      'p(cor)'     'equivk'     'p(unc)'     'p(FWE-cor)'    'p(FDR-cor)'    'T'        'equivZ'    'p(unc)'    'x,y,z {mm}'
            tab=sprintf('%s\t','');
            columns=[11 4 3 6 7 8];
            cols=length(columns);
            obj.renewSelection(cutoff_p_cluster);
            table={};
            for c=1:cols
                col=obj.getColumnText(columns(c));
                table=CommonMethods.insertColumnToTable(table,col, -1, tab);
            end
        end
        
        function renewSelection(obj,cutoff_p_cluster)
            if(isempty(obj.tabData.dat))
                selection=[];
            else
                selection_c=cellfun(@(x) x, obj.tabData.dat(:,3))<=cutoff_p_cluster;%selecting the local peaks in the selected clusters only
                selection_cPeak=obj.tabData.clusterPeak';
                %extending the cluster peak selection to include up to two additional
                %local peaks from each cluster
                idx=find(selection_cPeak);
                len=length(idx);
                for i=1:len
                    iI=idx(i);
                    if(i<len)
                        iF=min(iI+2,idx(i+1)-1);
                    else
                        iF=min(iI+2,len);
                    end
                    selection_cPeak(iI:iF)=true;
                end
                selection_cPeak=(selection_c.*selection_cPeak)>0;%selecting the cluster peaks of the selected clusters only
                selection_v=cellfun(@(x) x, obj.tabData.dat(:,6))<0.05;%selecting the local peaks meets voxel level FRW threshold
                selection=(selection_v.*selection_c+selection_cPeak)>0;%in selected clusters, selecting local peaks either meet the voxel level criterion or being the cluster peaks.
            end
            obj.tabData.selection=selection;
        end
        
        function coors=getPeakCoors(obj)
            if(isempty(obj.tabData.dat))
                coors=[];
            else
                coors=obj.tabData.dat(obj.tabData.selection,11);
                coors=[cellfun(@(x)x(1), coors) cellfun(@(x)x(2), coors) cellfun(@(x)x(3), coors)];
            end
        end
        
        function col=getColumnText(obj,c)
            col=obj.tabData.hdr(:,c);
            rows=size(obj.tabData.dat,1);
            if(isfield(obj.tabData, 'selection'))
                selection=obj.tabData.selection;
            else
                selection=true(1,rows);
            end
            fmt=obj.tabData.fmt;
            fmt([1 3 5 6 7 10])={'%0.12f'};
            for r=1:rows
                if(~selection(r))
                    continue;
                end
                v=obj.tabData.dat{r,c};
                len=length(v);
                strT='';
%                 for t=1:len
%                     cmd=['sprintf( ''' obj.tabData.fmt{c} ''', v)'];
%                     str=eval(cmd);
%                     strT=[strT ' ' str];
%                 end
                cmd=['sprintf( ''' fmt{c} ''', v)'];
                strT=eval(cmd);
                if(isempty(strT))
                    strT=' ';
                end
                col{end+1}=strT;
            end
        end
        
        function coor=getClusterPeakCoor(obj,cInd)
            inds=find(obj.tabData.clusterPeak);
            ind=inds(cInd);
            coor=obj.tabData.coors(ind,:);
        end
    end
    methods (Static = true)
        function tabData=getEmptyTabData()
            load(fullfile(CommonEnvMethods.getMatlabProjectPath(),'MRI', 'Scripts', 'Modified Scripts', 'saved_objects','emptySPMTabData.mat'));
            tabData.dat={};
        end
        
    end
end

