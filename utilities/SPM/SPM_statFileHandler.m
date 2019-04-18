classdef SPM_statFileHandler < handle
    %SPM_STATFILEHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fname;%a spm T/F contrast image
        volIndexes;%voxel indexes of the voxels; N*3 matrix
        coors;%voxel coordinates; Nx3 matrix
        intensity;%voxel values
        volIndexes1D;%1 D (single-column) voxel indexes;
        clusterInds;%the cluster index of each voxel, will be rearranged according to the clusters in tabData
        M;% the rotaion/translation matrix
        DIM;% dimension
        TF; %t test or f test? 'T','F' or []
        df;% degree of freedome for T/F test
        tdaHandler;%an instance of TDAtlasHandler
        cutoff_c;%the cluster level significance cutoff
        cutoff_v;%voxel level cutoff p value
        cutoff_FWE;%cutoff p value for local maxima output
        indexes_c;%one column vovel indexes
        %indexes_c is a cellarray of one-column voxel indexes. each element is a vector of indexes.
        %the length of the vector is the cluster size, and each element in
        %the vector is an index pointing to a position (row) in cor or
        %cor_singlecolumn.
        atlasIndexes;%the atlas indexes of each voxel. It is a N*6 integer matrix
        %atlasIndexes1D;%single-column atlas indexes of each voxel;
        sideness;%the sideness of each voxel: 1 left, 2 center, 3 right
        defined;%whether the voxel is inside any defined structures in the atlas
        N; %the number of selected voxels.
        xjViewHandles;%an xjView handles
        intensityCutoff;%voxel level intensity cutoff, function of pValue and df
        summaryStruct;%the summary tables and default file names
        tabData;%the tab separated data read in by spm containing all stat information, keeping it to use lagecy functions
        tbdHandler;%the tabData handler
        SPMmat;%the SPM.mat file
    end
    
    methods
        function obj = SPM_statFileHandler(xjViewFigHandle, mask, SPMmat, cutoff_c, cutoff_v, cutoff_FWE)
            %mask: An spm T or F map, defined by TF
            %SPMmat: file name of the "SPM.mat" file containing matching
            %mask
            %cutoff_c: cutoff for FWE corrected cluster level p value
            %cutoffv: cutoff for uncorrected voxel level p value
            
            %             if(~exist('TF', 'var'))
            %                 [~,name,~]=fileparts(maskg wiht spmT or spmF is expected for the class constructure SPM_statFileHandler');
            %                 end
            %             end
            
            if(~exist('mask','var'))
                mask=[];
            end
            if(~exist(mask,'file'))
                mask=fullfile('/home2/backup/home1/data/IRB13278/data/DTI/group/FA/twowayANOVA_flexible_repeatedMeasure/mask_spmICV','spmT_0001.nii');
            end
            if(~exist('SPMmat','var'))
                SPMmat=[];
            end
            if(isempty(SPMmat))
                folder=fileparts(mask);
                SPMmat=fullfile(folder,'SPM.mat');
            end
            if(~exist('xjViewFigHandle', 'var'))
                xjViewFigHandle=[];
            end
            obj.SPMmat=SPMmat;
            if(~xjViewHandler.isvalidhandle(xjViewFigHandle))
                xjViewFigHandle=xjview;
            end
             objt=xjViewHandler(mask, xjViewFigHandle);
             objt.selectPositiveOnly();
             xjViewHandles=objt.getXjviewHandles();
            
            if(~exist('cutoff_c','var'))
                cutoff_c=0.05;
            end
            if(~exist('cutoff_v','var'))
                cutoff_v=0.001;
            end
            if(~exist('cutoff_FWE','var'))
                cutoff_FWE=0.05;
            end
            obj.fname=mask;
            obj.xjViewHandles=xjViewHandles;
            [cor, intensity, cor_singlecolumn,M,DIM,TF,df] = obj.mask2coord(mask);
            %            obj.intensityCutoff=xjViewHandles.intensityThreshold{1};
            obj.intensityCutoff=CommonMethods.p2t(cutoff_v,df,TF);
            selection=(intensity>obj.intensityCutoff) + (intensity < -obj.intensityCutoff);
            selection=selection>0;
            obj.volIndexes=cor(selection,:);
            obj.intensity=intensity(selection);
            obj.volIndexes1D=cor_singlecolumn(selection);
            obj.coors=CommonMethods.ind2coor_mat(obj.volIndexes, M);
            obj.M=M;
            obj.DIM=DIM;
            obj.TF=TF;
            obj.df=df;
            obj.cutoff_c=cutoff_c;
            obj.cutoff_v=cutoff_v;
            obj.cutoff_FWE=cutoff_FWE;
            obj.buildClusters();
            obj.apply_P_cFWE_Cutoff();
            obj.tdaHandler=TDAtlasHandler();
            [atlasIndexes,sideness,defined]=obj.tdaHandler.getAtlasIndexes(obj.coors);
            obj.atlasIndexes=atlasIndexes;
            obj.sideness=sideness;
            obj.defined=defined;
        end
        
        function tbl=getOutputTableText(obj)
            %CommonMethods.applyingLineIndentation(line, level, dilimiter)
            tbl=obj.tbdHandler.getSummaryTable(obj.cutoff_c);
            coors=obj.tbdHandler.getPeakCoors();
            [structures, atlasIndexes]=obj.tdaHandler.findStructures(coors);
            tableNames=arrayfun(@(x)obj.tableRoiName(structures(x,:)),1:size(structures,1),'UniformOutput',false);
            col=horzcat({' ' 'Region'}, tableNames);
            tab=sprintf('%s\t','');
            tbl=CommonMethods.insertColumnToTable(tbl, col,-1,tab);
        end
        function exportOutputTable(obj, fname)
            fid=fopen(fname,'wt');
            tbl=obj.getOutputTableText();
            for i=1:length(tbl)
                fprintf(fid,'%s\n',tbl{i});
            end
            fclose(fid);
        end
        function [cor, intensity, cor_singlecolumn,M,DIM,TF,df] = mask2coord(obj,mask)
            %Modified from xjView
            % [cor, intensity, cor_singlecolumn] = mask2coord(mask, checkdimension)
            %
            % This is to retrieve the coordinate of a mask file, or a matrix of 3-D
            %
            % mask: an image file or a matrix (3-D), with some of the elements are
            % non-zeros
            % checkdimension: check if the dimension is checkdimension, if not, return empty
            % matrix
            % cor: a N*3 matrix, which each row a coordinate in matrix
            % intensity: a N*1 matrix, which encodes the intensity of each voxel.
            % cor_singlecolumn: a N*1 matrix, each row is the index in the matrix
            % M: rotation matrix
            % DIM: dimension
            % TF: t test or f test? 'T','F' or []
            % df: degree of freedome for T/F test
            %
            % Example:
            %   mask = zeros(4,3,2);
            %   mask(1,2,1) = 1;
            %   mask(3,2,2) = 1;
            %   mask2coord(mask)
            %
            %   mask2coord('spmT_0002.img')
            %   mask2coord('spmT_0002.img',[41 48 35])
            %
            % Xu Cui
            % 2004-9-20
            % last revised: 2005-04-30
            
            if nargin < 2
                checkdimension = 0;
            end
            
            if ischar(mask)
                V = spm_vol(mask);
                mask = spm_read_vols(V);
                M = V.mat;
                DIM = V.dim;
                TF = 'T';
                T_start = strfind(V.descrip,'SPM{T_[')+length('SPM{T_[');
                if isempty(T_start); T_start = strfind(V.descrip,'SPM{F_[')+length('SPM{F_['); TF='F'; end
                if isempty(T_start)
                    TF=[]; df=[];
                else
                    T_end = strfind(V.descrip,']}')-1;
                    df = str2num(V.descrip(T_start:T_end));
                end
            else
                M = [];
                TF = [];
                df = [];
            end
            
            dim = size(mask);
            %             if length(checkdimension)==3
            %                 if dim(1)~= checkdimension(1) | dim(2)~= checkdimension(2) | dim(3)~= checkdimension(3)
            %                     y = [];
            %                     disp('dimension doesn''t match')
            %                     return
            %                 end
            %             end
            
            %Taihao Jin 11/10/2016
            if(V.dt(1)==16)
                pos = find(~isnan(mask));
            else
                pos = find(mask~=0);
            end
            
            cor_singlecolumn = pos;
            intensity = mask(pos);
            
            y = zeros(length(pos),3);
            
            y(:,3) = ceil(pos/(dim(1)*dim(2)));
            pos = pos - (y(:,3)-1)*(dim(1)*dim(2));
            y(:,2) = ceil(pos/dim(1));
            pos = pos - (y(:,2)-1)*(dim(1));
            y(:,1) = pos;
            
            cor = y;
            DIM = dim;
        end
        
        function apply_P_cFWE_Cutoff(obj)
            %this function select the clusters based on the FWE adjusted
            %cluster level signification.
            p_cluster_c=obj.cutoff_c;  
            p_cluster=arrayfun(@(x)obj.tbdHandler.getClusterP(x),obj.clusterInds);
            selected=p_cluster<=p_cluster_c;
            obj.indexes_c=obj.indexes_c(selected);            
            obj.clusterInds=obj.clusterInds(selected);            
        end
        
        function buildClusters(obj)
            if(isempty(obj.intensity))
                tabData=spmTabDataHandler.getEmptyTabData();
            else
                tabData=obj.xjViewHandles.volumePush.Callback(obj.xjViewHandles.volumePush,'export50',obj.SPMmat);
            end
            
            obj.tabData=tabData;
            obj.tbdHandler=spmTabDataHandler(tabData,spm_vol(obj.fname));
                       
            volIndexes=obj.volIndexes;%3xlen 3D indexes
            %the one D index
            clusterIndexes=spm_clusters(volIndexes',26);
            %the index indicating the cluster of the voxels belong to
            
            clusters=unique(clusterIndexes);
            lenc0=max(clusters);%the number of clusters
            
            %counting the cluster sizes
            sizes0=zeros(1,lenc0);
            indexes0=cell(1,lenc0);
            for i=1:length(clusters);
                indc=clusters(i);
                inds=find(clusterIndexes==indc);
                sizes0(i)=length(inds);
                indexes0{i}=inds;
            end
            sizeCutoff=str2num(obj.xjViewHandles.clusterSizeThresholdEdit.String);
            
            indexes0=indexes0(sizes0>sizeCutoff);
            
            lenc=length(indexes0);
            
            indexesc=cell(1,lenc);
            obj.clusterInds=nan(1,lenc);
            for i=1:lenc
                indexes=indexes0{i};
                indc=obj.tbdHandler.getClusterInd(obj.volIndexes1D(indexes));
                indexesc{i}=indexes;
                obj.clusterInds(i)=indc;
            end
            selection=~isnan(obj.clusterInds);
            obj.indexes_c=indexesc(selection);
            obj.clusterInds=obj.clusterInds(selection);
        end
        
        function exportClusterAnatomicalMasks(obj, level)
            %this function export the intersecting masks between the
            %clusters and the anatomical regions in TD atlas
            %level: the level in the atlas - 1 for the hemisphere level , 2 for the lobe level,
            %3 for the gyrus level, 4 for the tissue level (GM/WM), 5 for the cell
            %level (BA, nuclei), 6 for aal.
            obj.fname=obj.xjViewHandles.imageFileName{1};
            [folder,name,ext]=fileparts(obj.fname);
            folder=fullfile(folder,'maskImgs',name);
            if(~exist(folder,'file'))
                mkdir(folder);
            end
            
            if(~exist('level','var'))
                level=3;
            end
            len=length(obj.indexes_c);
            V=spm_vol(obj.fname);
            for c=1:len
                inds=obj.indexes_c{c};
                idx=obj.atlasIndexes(inds,level);
                uas=unique(idx);%unique atlas indexes for the level
                for u=1:length(uas)
                    aid=uas(u);%the atlas index
                    if(aid==0)
                        continue;
                    end
                    selected=idx==aid;
                    inds_selected=inds(selected);
                    vInds1D=obj.volIndexes1D(inds_selected);%the 1D voxel indexes of the voxels in the cluster belong to the anatomical region
                    sname=obj.tdaHandler.getStructureName(level,aid);
                    vols=zeros(obj.DIM);
                    vols(vInds1D)=1;
                    fname=fullfile(folder,[name '_' sname '_' num2str(length(vInds1D)) 'voxels.nii']);
                    V.fname=fname;
                    spm_write_vol(V,vols);
                end
            end
        end
        
        function exportClusterMasks(obj)
            %this function export the the clusters mask
            fname=obj.xjViewHandles.imageFileName{1};
            [folder,name,ext]=fileparts(fname);
            folder=fullfile(folder,'maskImgs');
            if(~exist(folder,'file'))
                mkdir(folder);
            end
            len=length(obj.indexes_c);
            V=spm_vol(fname);
            V.dt(1)=2;
            vols=zeros(V.dim);
            for c=1:len
                %                cInd=obj.clusterInds(c);
                %               coor=obj.tbdHandler.getClusterPeakCoor(cInd);
                %               cStr=CommonMethods.coor2str(coor);
                inds=obj.indexes_c{c};
                inds=obj.volIndexes1D(inds);
                vols(inds)=1;
            end
            V.fname=fullfile(folder,[name '_clusters' ext]);
            spm_write_vol(V,vols);
        end
        
        function exportSummaryTables(obj)
            %exports four files about the informations of the local maxima
            %in the same directory of mask
            
            tabData=obj.tabData;
            if(~exist('hObject','var'))
                hObject=[];
            end
            
            rows=size(tabData.dat,1);
            cols=size(tabData.dat,2);
            p_cluster_c=obj.cutoff_c;
            p_FWE_c=obj.cutoff_FWE;
            coors=cellfun(@(x)x',tabData.dat(:,11),'UniformOutput',false);
            coors=cell2mat(coors);
            p_FWE=arrayfun(@(x)tabData.dat{x,6},1:rows);
            p_FDR=arrayfun(@(x)tabData.dat{x,7},1:rows);
            Ts=arrayfun(@(x)tabData.dat{x,8},1:rows);
            p_cluster=cellfun(@(x)x,tabData.dat(:,3),'UniformOutput',false);
            
            inds=find(cellfun(@(x)~isempty(x),p_cluster));%index of the first voxel of each cluster
            p_cluster=arrayfun(@(x)tabData.dat{x,3},inds);
            p_cluster_u=arrayfun(@(x)tabData.dat{x,5},inds);%uncorrected cluster p values
            clusterSizes=arrayfun(@(x)tabData.dat{x,4},inds);
            
            clusterIndexes=zeros(rows,1);
            for i=1:length(inds)
                iI=inds(i);
                if(i<length(inds))
                    iF=inds(i+1)-1;
                else
                    iF=rows;
                end
                for it=iI:iF
                    clusterIndexes(it)=i;
                end
            end
            selected_c=find(p_cluster<p_cluster_c);%selectedClusters
            %   rows=size(tabData.dat,1)+12;%12 extra lines are needed for the header info.
            minClusterSize=min(clusterSizes(selected_c));
            TDH=TDAtlasHandler();
            [AtlasIndexes,sideness,defined]=TDH.getAtlasIndexes(coors);
            selected=zeros(rows,1);
            for c=1:length(selected_c)
                cind=selected_c(c);
                cur_inds=find(clusterIndexes==cind.*defined);
                selected(cur_inds(1))=1;
                pvs=p_FWE(cur_inds);
                cur_selectedInds=cur_inds(pvs<p_FWE_c);
                selected(cur_selectedInds)=1;
            end
            coors_selected=coors(find(selected),:);
            TDH.updateSelection(coors_selected);
            
            extraCols=length(size(AtlasIndexes,2));
            cols=cols+extraCols;
            
            curImg=obj.fname;
            
            pst=['p le ' num2str(obj.cutoff_v)];
            kst=['k ge ' num2str(minClusterSize)];
            [cur_dir,name,ext]=fileparts(curImg);
            cur_name=['Significant Local Maxima_' name '_' pst '_' kst]
            %    fname=uigetfile(fullfile(cur_dir,'*.tsv'),'select file');
            %             [fname, path]=uiputfile(fullfile(cur_dir,'*.tsv'),'Specifya tab separated file to export significant local maxima', cur_fname);
            
            ext='.tsv';
            
            path=cur_dir;
            
            fname_all=[cur_name '_all' ext];
            fname_cluster=[cur_name '_cluster' ext];
            fname_table=[cur_name '_table' ext];
            fname_unique=[cur_name '_unique' ext];% saving only one local maximum from each anatomic region
            
            fname_all=fullfile(path,fname_all);
            fname_table=fullfile(path, fname_table);
            
            fid_all=fopen(fname_all,'wt');%export all local maxima
            fid_unique=fopen(fname_unique,'wt');%export each anatomic region once for each cluster
            fid_cluster=fopen(fname_cluster,'wt');%export only the cluster peack
            fid_table=fopen(fname_table, 'wt');%export as the table format
            
            tab=sprintf('%s\t','');
            
            line=tabData.tit;
            for i=2:cols
                line=[line tab ''];
            end
            fprintf(fid_all,'%s\n',line);
            fprintf(fid_unique,'%s\n',line);
            fprintf(fid_cluster,'%s\n',line);
            
            line=tabData.str;
            for i=2:cols
                line=[line tab ''];
            end
            fprintf(fid_all,'%s\n',line);
            fprintf(fid_unique,'%s\n',line);
            fprintf(fid_cluster,'%s\n',line);
            
            cellStr=tabData.ftr;
            rowst=size(cellStr,1);
            colst=size(cellStr,2);
            for r=1:rowst
                line=cellStr{r,1};
                for c=2:colst
                    line=[line tab cellStr{r,c}];
                end
                for c=colst+1:cols
                    line=[line tab ''];
                end
                fprintf(fid_all,'%s\n',line);
                fprintf(fid_unique,'%s\n',line);
                fprintf(fid_cluster,'%s\n',line);
            end
            
            cellStr=tabData.hdr;
            
            %the second row of the header is in wrong order and now is correcting
            %it
            st=cellStr{2,1};
            cellStr{2,1}=cellStr{2,2};
            cellStr{2,2}=st;
            %end of correction
            
            rowst=size(cellStr,1);
            colst=size(cellStr,2);
            for r=1:rowst
                line=cellStr{r,1};
                for c=2:colst
                    line=[line tab cellStr{r,c}];
                end
                for c=colst+1:cols
                    line=[line tab ''];
                end
                fprintf(fid_all,'%s\n',line);
                fprintf(fid_unique,'%s\n',line);
                fprintf(fid_cluster,'%s\n',line);
            end
            
            %exporting table header
            line=['MNI coordinates (x y z)' tab 'Cluster extent (k)' tab 'Cluster-level' tab 'Voxel-level' tab 'Voxel-level' tab 'T' tab 'Region description (selected local maxima)' tab 'peak counts'];
            fprintf(fid_table, '%s\n', line);
            table_col_inds=[4 3 6 7 8];
            
            line=['' tab '' tab 'P_corrected' tab 'P_FWE-corrected' tab 'P_FDR-corrected' tab '' tab ''];
            fprintf(fid_table, '%s\n', line);
            
            %writing the data
            cellStr=tabData.dat;
            rowst=size(cellStr,1);
            colst=size(cellStr,2);
            
            pcutoff=1.05;
            
            onelinestructures={};
            tableROIShortNames={};
            peakCounts=[];
            tableLines={};
            tableROIIndexes=[];
            for r=1:rowst
                uq=false;%whether a unique anatomic rigion in the cluster
                if(r==1)
                    setStr=cellStr(r,1:2);
                else
                    cellStr(r,1:2)=setStr;
                end
                
                coor=cellStr{r,colst};
                cellarraystructure = obj.tdaHandler.findStructures(coor');
                %       onelinestructure=onelinestructure{1};
                
                %        tableRoiName=tableRoiShortName(cellarraystructure);
                tableRoiName=obj.tableRoiName(cellarraystructure);
                unspecified=false;
                lastInCluster=false;
                if(~isempty(strfind(tableRoiName,'undefined')))
                    unspecified=true;
                    if(r<rowst)
                        if(isempty(cellStr{r+1,3})&&~isempty(cellStr{r,3}))
                            cellStr(r+1,3:5)=cellStr(r,3:5);
                            lastInCluster=false;
                        else
                            lastInCluster=true;
                        end
                    end
                end
                
                onelinestructure='';
                for i=1:length(cellarraystructure)-1
                    onelinestructure=[onelinestructure cellarraystructure{i}];
                end
                
                if(~isempty(cellStr{r,3}))%first local maximum of the cluster
                    clusterStr=cellStr(r,3:5);
                    onelinestructures={onelinestructure};
                    uniqueStructure=true;
                    firstInCluster=true;
                    uq=true;
                else
                    firstInCluster=false;
                    cellStr(r,3:5)=clusterStr;
                    if(ismember(onelinestructure,onelinestructures))
                        uniqueStructure=false;
                    else
                        uniqueStructure=true;
                        onelinestructures{end+1}=onelinestructure;
                    end
                end
                
                line=num2str(cellStr{r,1});
                for c=2:colst-1
                    line=[line tab num2str(cellStr{r,c})];
                end
                
                st_coor='';
                for t=1:3
                    st_coor=[st_coor ' ' num2str(coor(t))];
                end
                line=[line tab st_coor];
                
                for c=1:extraCols
                    line=[line tab cellarraystructure{c}];
                end
                fprintf(fid_all,'%s\n',line);
                
                %         pfwec=str2double(cellStr{r,6});
                %         pfdrc=str2double(cellStr{r,7});
                %         pcluster=str2double(cellStr{r,3});
                pfwec=cellStr{r,6};
                pfdrc=cellStr{r,7};
                pcluster=cellStr{r,3};
                if(min(pfwec,pfdrc)<pcutoff)
                    if(uniqueStructure)
                        uq=true;
                    end
                elseif(firstInCluster&&pcluster<pcutoff)
                    uq=true;
                end
                
                if(uq)
                    fprintf(fid_unique,'%s\n', line);
                end
                if(firstInCluster&&pcluster<pcutoff)
                    fprintf(fid_cluster,'%s\n', line);
                end
                
                %exporting to the table
                if(unspecified)
                    if(~lastInCluster)
                        continue;
                    end
                end
                
                line_table=st_coor;
                frmt={'%6.0f' '%7.4f' '%7.4f' '%7.4f' '%7.2f'};
                for i=1:length( table_col_inds)
                    cmd=['sprintf(''' frmt{i} ''', cellStr{r,table_col_inds(i)})'];
                    line_table=[line_table tab eval(cmd)];
                end
                
                %exporting to the table file
                roiIndex=0;
                if((~iscell(tableRoiName)&&~ischar(tableRoiName))||isempty(tableRoiName))
                    uniqueStructure=false;
                elseif(isempty(tableROIShortNames))
                    uniqueStructure=true;
                    tableROIShortNames{end+1}=tableRoiName;
                    peakCounts(end+1)=1;
                    roiIndex=length(peakCounts);
                else
                    [~, pos]=ismember(tableRoiName, tableROIShortNames);
                    if(any(pos))
                        uniqueStructure=false;
                        peakCounts(pos(1))=peakCounts(pos(1))+1;
                        roiIndex=pos(1);
                    else
                        uniqueStructure=true;
                        tableROIShortNames{end+1}=tableRoiName;
                        peakCounts(end+1)=1;
                        roiIndex=length(peakCounts);
                    end
                end
                line_table=[line_table tab tableRoiName];
                export=false;
                %        if(min(pfwec,pfdrc)<pcutoff) %there are too many false discovery
                %        corrected p values lower than the cutoff
                %        if(pfdrc<pcutoff)
                if(pfwec<pcutoff)
                    if(uniqueStructure)
                        export=true;
                    end
                elseif(firstInCluster&&pcluster<pcutoff)
                    export=true;
                end
                
                if(export)
                    %            fprintf(fid_table,'%s\n', line_table);
                    tableLines{end+1}=line_table;
                    tableROIIndexes(length(tableLines))=roiIndex;
                end
            end
            
            fclose(fid_all);
            fclose(fid_unique);
            fclose(fid_cluster);
            
            for l=1:length(tableLines)
                roiIndex=tableROIIndexes(l);
                tableLines{l}=[tableLines{l} tab num2str(peakCounts(roiIndex))];
            end
            
            names=arrayfun(@(x)tableROIShortNames{x},tableROIIndexes,'UniformOutput',false);
            shorterNames=cellfun(@(x)strrep(x,'L ',''),names,'UniformOutput',false);
            shorterNames=cellfun(@(x)strrep(x,'R ',''),shorterNames,'UniformOutput',false);
            
            len=length(tableLines);
            sNodes=TDH.getSelectedStructureNodes(1);
            
            for l=1:len-2
                sname=shorterNames{l};
                [~,pos]=ismember(sname,shorterNames(l+1:len));
                lent=length(pos);
                if(lent>0)
                    p=pos(1)+l;
                    if(p>l+1)
                        tmp1=shorterNames{p};
                        tmp2=tableLines{p};
                        
                        for pt=p:-1:l+2
                            shorterNames{pt}=shorterNames{pt-1};
                            tableLines{pt}=tableLines{pt-1};
                        end
                        
                        shorterNames{l+1}=tmp1;
                        tableLines{l+1}=tmp2;
                    end
                end
            end
            for l=1:length(tableLines)
                roiIndex=tableROIIndexes(l);
                fprintf(fid_table,'%s\n', tableLines{l});
            end
            fclose(fid_table);
        end
        
        function exportAnatomicalVoxelCounts(obj)
            %Part of externalFunction(hObject, action,varargin) that I have
            %added into xjview. Have not yet fulled adapted into this
            %class.
            %20171014 
            msg='the function exportAnatomicalVoxelCounts in SPM_statFileHandler has not fulled implementd.';
            fprintf(1,'%s\n',msg);
            return;
            
            X=varargin{1};
            if(X.erode)
                DB=X.DB_sorted_eroded;
            else
                DB=X.DB_sorted;
            end
            cur_dir=varargin{2};
            cur_fname=varargin{3};
            mat=st.vols{1}.mat;
            mni=handles.currentDisplayMNI{1};
            coors=CommonMethods.mni2ind_cols(mni',mat)';
            [coors,keep]=CommonMethods.removeDuplicatedRows(coors);
            mni=mni(find(keep),:);
            xcoors=mni(:,1);
            left=xcoors<=-0.5;
            right=xcoors>=0.5;
            center=~(left+right);
            coorsl=coors(left,:);
            coorsr=coors(right,:);
            coorsc=coors(center,:);
            noVoxelsTobal=length(mni);
            
            len=length(mni);
            [fname, path]=uiputfile(fullfile(cur_dir,'*.tsv'),'Specifya tab separated file to export anatomical region voxel counts', cur_fname);
            fname=fullfile(path,fname);
            
            len1=length(DB);
            VoxelCounts=cell(1,len1);
            rows=0;
            for d=1:len1
                db=DB{d};
                anatomys=db.anatomy;
                anatomys=cellfun(@(x)strrep(x,' ','_'), anatomys, 'UniformOutput', false);
                counts.type=db.type;
                counts.anatomy=db.anatomy;
                counts.noVoxels=zeros(length(anatomys),4);
                typesl=arrayfun(@(x)db.mnilist(coorsl(x,1),coorsl(x,2),coorsl(x,3)),1:size(coorsl,1));
                typesr=arrayfun(@(x)db.mnilist(coorsr(x,1),coorsr(x,2),coorsr(x,3)),1:size(coorsr,1));
                typesc=arrayfun(@(x)db.mnilist(coorsc(x,1),coorsc(x,2),coorsc(x,3)),1:size(coorsc,1));
                for i=1:length(anatomys)
                    anatomy=anatomys{i};
                    counts.noVoxels(i,1)=nnz(typesl==i);
                    counts.noVoxels(i,2)=nnz(typesr==i);
                    counts.noVoxels(i,3)=nnz(typesc==i);
                    if(i==42)
                        i=i;
                    end
                    counts.noVoxels(i,4)=nnz(db.mnilist==i);
                end
                VoxelCounts{d}=counts;
                rows=max(rows,length(anatomys));
            end
            
            fid=fopen(fname,'wt');
            cols=len1*6;
            line=['Total number of voxels: ' num2str(size(coors,1))];
            for i=2:cols
                line=[line tab ''];
            end
            fprintf(fid,'%s\n',line);
            
            line='';
            for i=1:len1
                line=[line tab VoxelCounts{i}.type tab ' left' tab ' right' tab ' center' tab ' size' tab ' porsion'];
            end
            fprintf(fid,'%s\n',line);
            
            for r=1:rows
                line='';
                for i=1:len1
                    counts=VoxelCounts{i};
                    anatomys=counts.anatomy;
                    lent=length(anatomys);
                    if(r<=lent)
                        linet=anatomys{r};
                        for t=1:4
                            linet=[linet tab num2str(counts.noVoxels(r,t))];
                        end
                        linet=[linet tab num2str(sum(counts.noVoxels(r,1:3))/counts.noVoxels(r,4))];
                    else
                        linet='';
                        for t=1:5
                            linet=[linet tab ''];
                        end
                    end
                    line=[line tab linet];
                end
                fprintf(fid,'%s\n',line);
            end
            
            fclose(fid);
        end
        function makeSummaryTables_xjView(obj, hObject)
            %ollder version of of exportSummaryTables.
            %hObject is a spmExtension GUI handle
            %    exportDB(guiData.DB);
            %             chandler=SPM_statFileHandler(guiData.xjViewStruct.handles,guiData.currentImgFileName);
            %             chandler.exportClusterAnatomicalMasks(3);
            %
            tabData=obj.tabData;
            if(~exist('hObject','var'))
                hObject=[];
            end
            
            if(~isempty(hObject))
                guiData=guidata(hObject);
            else
                guiData=obj.xjViewHandles;
            end
            %    [onelinestructure, cellarraystructure] = cuixuFindStructure(tabData.dat{1,11}', guiData.X.DB);
            
            rows=size(tabData.dat,1);
            cols=size(tabData.dat,2);
            p_cluster_c=0.05;
            p_FWE_c=0.01;
            coors=cellfun(@(x)x',tabData.dat(:,11),'UniformOutput',false);
            coors=cell2mat(coors);
            p_FWE=arrayfun(@(x)tabData.dat{x,6},1:rows);
            p_FDR=arrayfun(@(x)tabData.dat{x,7},1:rows);
            Ts=arrayfun(@(x)tabData.dat{x,8},1:rows);
            p_cluster=cellfun(@(x)x,tabData.dat(:,3),'UniformOutput',false);
            inds=find(cellfun(@(x)~isempty(x),p_cluster));%index of the first voxel of each cluster
            p_cluster=arrayfun(@(x)tabData.dat{x,3},inds);
            p_cluster_u=arrayfun(@(x)tabData.dat{x,5},inds);%uncorrected cluster p values
            clusterSizes=arrayfun(@(x)tabData.dat{x,4},inds);
            clusterIndexes=zeros(rows,1);
            for i=1:length(inds)
                iI=inds(i);
                if(i<length(inds))
                    iF=inds(i+1)-1;
                else
                    iF=rows;
                end
                for it=iI:iF
                    clusterIndexes(it)=i;
                end
            end
            selected_c=find(p_cluster<p_cluster_c);%selectedClusters
            %   rows=size(tabData.dat,1)+12;%12 extra lines are needed for the header info.
            TDH=TDAtlasHandler();
            [AtlasIndexes,sideness,defined]=TDH.getAtlasIndexes(coors);
            selected=zeros(rows,1);
            for c=1:length(selected_c)
                cind=selected_c(c);
                cur_inds=find(clusterIndexes==cind.*defined);
                selected(cur_inds(1))=1;
                pvs=p_FWE(cur_inds);
                cur_selectedInds=cur_inds(pvs<p_FWE_c);
                selected(cur_selectedInds)=1;
            end
            coors_selected=coors(find(selected),:);
            TDH.updateSelection(coors_selected);
            
            extraCols=length(size(AtlasIndexes,2));
            cols=cols+extraCols;
            if(isfield(guiData,'currentImgileName'))
                curImg=guiData.currentImgFileName;
            else
                curImg=obj.fname;
            end
            pst=['p le ' num2str(guiData.xjViewStruct.handles.pValue)];
            kst=['k ge ' num2str(guiData.xjViewStruct.handles.clusterSizeThreshold)];
            [cur_dir,name,ext]=fileparts(curImg);
            cur_fname=['Significant Local Maxima_' name '_' pst '_' kst '.tsv'];
            %    fname=uigetfile(fullfile(cur_dir,'*.tsv'),'select file');
            [fname, path]=uiputfile(fullfile(cur_dir,'*.tsv'),'Specifya tab separated file to export significant local maxima', cur_fname);
            [~, name, ext]=fileparts(fname);
            fname_all=[name '_all' ext];
            fname_cluster=[name '_cluster' ext];
            fname_table=[name '_table' ext];
            
            fname_all=fullfile(path,fname_all);
            fname_unique=fullfile(path, fname);% saving only one local maximum from each anatomic region
            fname_table=fullfile(path, fname_table);
            
            fid_all=fopen(fname_all,'wt');%export all local maxima
            fid_unique=fopen(fname_unique,'wt');%export each anatomic region once for each cluster
            fid_cluster=fopen(fname_cluster,'wt');%export only the cluster peack
            fid_table=fopen(fname_table, 'wt');%export as the table format
            
            tab=sprintf('%s\t','');
            
            line=tabData.tit;
            for i=2:cols
                line=[line tab ''];
            end
            fprintf(fid_all,'%s\n',line);
            fprintf(fid_unique,'%s\n',line);
            fprintf(fid_cluster,'%s\n',line);
            
            line=tabData.str;
            for i=2:cols
                line=[line tab ''];
            end
            fprintf(fid_all,'%s\n',line);
            fprintf(fid_unique,'%s\n',line);
            fprintf(fid_cluster,'%s\n',line);
            
            cellStr=tabData.ftr;
            rowst=size(cellStr,1);
            colst=size(cellStr,2);
            for r=1:rowst
                line=cellStr{r,1};
                for c=2:colst
                    line=[line tab cellStr{r,c}];
                end
                for c=colst+1:cols
                    line=[line tab ''];
                end
                fprintf(fid_all,'%s\n',line);
                fprintf(fid_unique,'%s\n',line);
                fprintf(fid_cluster,'%s\n',line);
            end
            
            cellStr=tabData.hdr;
            
            %the second row of the header is in wrong order and now is correcting
            %it
            st=cellStr{2,1};
            cellStr{2,1}=cellStr{2,2};
            cellStr{2,2}=st;
            %end of correction
            
            rowst=size(cellStr,1);
            colst=size(cellStr,2);
            for r=1:rowst
                line=cellStr{r,1};
                for c=2:colst
                    line=[line tab cellStr{r,c}];
                end
                for c=colst+1:cols
                    line=[line tab ''];
                end
                fprintf(fid_all,'%s\n',line);
                fprintf(fid_unique,'%s\n',line);
                fprintf(fid_cluster,'%s\n',line);
            end
            
            %exporting table header
            line=['MNI coordinates (x y z)' tab 'Cluster extent (k)' tab 'Cluster-level' tab 'Voxel-level' tab 'Voxel-level' tab 'T' tab 'Region description (selected local maxima)' tab 'peak counts'];
            fprintf(fid_table, '%s\n', line);
            table_col_inds=[4 3 6 7 8];
            
            line=['' tab '' tab 'P_corrected' tab 'P_FWE-corrected' tab 'P_FDR-corrected' tab '' tab ''];
            fprintf(fid_table, '%s\n', line);
            
            %writing the data
            cellStr=tabData.dat;
            rowst=size(cellStr,1);
            colst=size(cellStr,2);
            
            pcutoff=1.05;
            
            onelinestructures={};
            tableROIShortNames={};
            peakCounts=[];
            tableLines={};
            tableROIIndexes=[];
            for r=1:rowst
                uq=false;%whether a unique anatomic rigion in the cluster
                if(r==1)
                    setStr=cellStr(r,1:2);
                else
                    cellStr(r,1:2)=setStr;
                end
                
                coor=cellStr{r,colst};
                [onelinestructure, cellarraystructure] = cuixuFindStructure(coor', guiData.X.DB);
                %       onelinestructure=onelinestructure{1};
                
                %        tableRoiName=tableRoiShortName(cellarraystructure);
                tableRoiName=tableRoiName(cellarraystructure);
                unspecified=false;
                lastInCluster=false;
                if(~isempty(strfind(tableRoiName,'undefined')))
                    unspecified=true;
                    if(r<rowst)
                        if(isempty(cellStr{r+1,3})&&~isempty(cellStr{r,3}))
                            cellStr(r+1,3:5)=cellStr(r,3:5);
                            lastInCluster=false;
                        else
                            lastInCluster=true;
                        end
                    end
                end
                
                onelinestructure='';
                for i=1:length(cellarraystructure)-1
                    onelinestructure=[onelinestructure cellarraystructure{i}];
                end
                
                if(~isempty(cellStr{r,3}))%first local maximum of the cluster
                    clusterStr=cellStr(r,3:5);
                    onelinestructures={onelinestructure};
                    uniqueStructure=true;
                    firstInCluster=true;
                    uq=true;
                else
                    firstInCluster=false;
                    cellStr(r,3:5)=clusterStr;
                    if(ismember(onelinestructure,onelinestructures))
                        uniqueStructure=false;
                    else
                        uniqueStructure=true;
                        onelinestructures{end+1}=onelinestructure;
                    end
                end
                
                line=num2str(cellStr{r,1});
                for c=2:colst-1
                    line=[line tab num2str(cellStr{r,c})];
                end
                
                st_coor='';
                for t=1:3
                    st_coor=[st_coor ' ' num2str(coor(t))];
                end
                line=[line tab st_coor];
                
                for c=1:extraCols
                    line=[line tab cellarraystructure{c}];
                end
                fprintf(fid_all,'%s\n',line);
                
                %         pfwec=str2double(cellStr{r,6});
                %         pfdrc=str2double(cellStr{r,7});
                %         pcluster=str2double(cellStr{r,3});
                pfwec=cellStr{r,6};
                pfdrc=cellStr{r,7};
                pcluster=cellStr{r,3};
                if(min(pfwec,pfdrc)<pcutoff)
                    if(uniqueStructure)
                        uq=true;
                    end
                elseif(firstInCluster&&pcluster<pcutoff)
                    uq=true;
                end
                
                if(uq)
                    fprintf(fid_unique,'%s\n', line);
                end
                if(firstInCluster&&pcluster<pcutoff)
                    fprintf(fid_cluster,'%s\n', line);
                end
                
                %exporting to the table
                if(unspecified)
                    if(~lastInCluster)
                        continue;
                    end
                end
                
                line_table=st_coor;
                frmt={'%6.0f' '%7.4f' '%7.4f' '%7.4f' '%7.2f'};
                for i=1:length( table_col_inds)
                    cmd=['sprintf(''' frmt{i} ''', cellStr{r,table_col_inds(i)})'];
                    line_table=[line_table tab eval(cmd)];
                end
                
                %exporting to the table file
                roiIndex=0;
                if((~iscell(tableRoiName)&&~ischar(tableRoiName))||isempty(tableRoiName))
                    uniqueStructure=false;
                elseif(isempty(tableROIShortNames))
                    uniqueStructure=true;
                    tableROIShortNames{end+1}=tableRoiName;
                    peakCounts(end+1)=1;
                    roiIndex=length(peakCounts);
                else
                    [~, pos]=ismember(tableRoiName, tableROIShortNames);
                    if(any(pos))
                        uniqueStructure=false;
                        peakCounts(pos(1))=peakCounts(pos(1))+1;
                        roiIndex=pos(1);
                    else
                        uniqueStructure=true;
                        tableROIShortNames{end+1}=tableRoiName;
                        peakCounts(end+1)=1;
                        roiIndex=length(peakCounts);
                    end
                end
                line_table=[line_table tab tableRoiName];
                export=false;
                %        if(min(pfwec,pfdrc)<pcutoff) %there are too many false discovery
                %        corrected p values lower than the cutoff
                %        if(pfdrc<pcutoff)
                if(pfwec<pcutoff)
                    if(uniqueStructure)
                        export=true;
                    end
                elseif(firstInCluster&&pcluster<pcutoff)
                    export=true;
                end
                
                if(export)
                    %            fprintf(fid_table,'%s\n', line_table);
                    tableLines{end+1}=line_table;
                    tableROIIndexes(length(tableLines))=roiIndex;
                end
            end
            
            fclose(fid_all);
            fclose(fid_unique);
            fclose(fid_cluster);
            
            for l=1:length(tableLines)
                roiIndex=tableROIIndexes(l);
                tableLines{l}=[tableLines{l} tab num2str(peakCounts(roiIndex))];
            end
            
            names=arrayfun(@(x)tableROIShortNames{x},tableROIIndexes,'UniformOutput',false);
            shorterNames=cellfun(@(x)strrep(x,'L ',''),names,'UniformOutput',false);
            shorterNames=cellfun(@(x)strrep(x,'R ',''),shorterNames,'UniformOutput',false);
            
            len=length(tableLines);
            sNodes=TDH.getSelectedStructureNodes(1);
            
            for l=1:len-2
                sname=shorterNames{l};
                [~,pos]=ismember(sname,shorterNames(l+1:len));
                lent=length(pos);
                if(lent>0)
                    p=pos(1)+l;
                    if(p>l+1)
                        tmp1=shorterNames{p};
                        tmp2=tableLines{p};
                        
                        for pt=p:-1:l+2
                            shorterNames{pt}=shorterNames{pt-1};
                            tableLines{pt}=tableLines{pt-1};
                        end
                        
                        shorterNames{l+1}=tmp1;
                        tableLines{l+1}=tmp2;
                    end
                end
            end
            for l=1:length(tableLines)
                roiIndex=tableROIIndexes(l);
                fprintf(fid_table,'%s\n', tableLines{l});
            end
            fclose(fid_table);
        end
        function name=tableRoiName(obj,cellname)
            %This function output one line anatomy name formatted for a
            %table
            %one line roi name
            name=cellname{1};
            switch (name)
                case 'Left Cerebrum'
                    name='L';
                case 'Right Cerebrum'
                    name='R';
            end
            
            if(strcmpi(cellname{3},'undefined'))
                name=[name ' ' cellname{2}];
            else
                name=[name ' ' cellname{3}];
            end
            
            if(~strcmpi(cellname{5},'undefined'))
                name=[name ' (' cellname{5} ')'];
            end
            
            name=strrep(name,'brodmann area', 'BA');
        end
        function name=tableRoiShortName(obj,cellname)
            %one line roi name
            name=cellname{1};
            switch (name)
                case 'Left Cerebrum'
                    name='L';
                case 'Right Cerebrum'
                    name='R';
            end
            
            if(strcmpi(cellname{3},'undefined'))
                name=[name ' ' cellname{2}];
            else
                name=[name ' ' cellname{3}];
            end
        end
        
        %     if(~strcmpi(cellname{5},'undefined'))
        %         name=[name ' (' cellname{5} ')'];
        %     end
    end
    
    methods (Static=true)
        function tbl=getSPMSummaryTable(fname)
            msg=['Getting SPM summary table ' fname ' at ' datestr(clock)];
            cprintf('blue', '%s\n', msg);
            if(exist('SPM', 'var'))
                clear(SPM);
            end
            load(fname);
            
            dil=sprintf('%s\t','');
            
            folder=SPM.swd;
            ind=strfind(folder, 'IRB13278');
            folder=fullfile(CommonEnvMethods.getCOHPath(),'data', folder(ind:length(folder)));
            folder=strrep(folder,'\',filesep);
            tbl={};
            fnameSPM=fullfile(folder,'SPM.mat');
            cons=SPM.xCon;
            len=length(cons);
            for c=1:len
                if(c<1)
                    continue;
                end
                CommonMethods.closeAll();
                con=cons(c);
                msg=['Getting stats summary table ' con.Vspm.fname ' at ' datestr(clock)];
                cprintf('blue', '%s\n', msg);
                obj=SPM_statFileHandler([],fullfile(folder,con.Vspm.fname),fnameSPM);
                tblt=obj.getOutputTableText();
                tblt=CommonMethods.applyingLineIndentationToTable(tblt,1,dil);
                tblt=CommonMethods.insertRowToTable(tblt,[con.name dil con.Vspm.fname],1,dil);
                tbl=CommonMethods.mergeTables(tbl,tblt,dil);
            end
            tbl=CommonMethods.applyingLineIndentationToTable(tbl,1,dil);
            tbl=CommonMethods.insertRowToTable(tbl,fnameSPM,1,dil);
        end
        function exportSPMSummaryTable(fname,fid)
            %fname: file name of the "SPM.mat"
            %fid: file id to which the summary table is exported.
                 
            msg=['Getting SPM summary table ' fname ' at ' datestr(clock)];
            cprintf('blue', '%s\n', msg);
            if(exist('SPM', 'var'))
                clear(SPM);
            end
            
            if(~exist('fid','var'))
                fid=-1;
            end
            
            closefile=false;
            if(fid<1)
                [folder,name,ext]=fileparts(fname);
                fid=fopen(fullfile(folder, 'spmSummary.tsv'),'wt');
                closefile=true;
            end
            
            load(fname);
            
            dil=sprintf('%s\t','');
            
            folder=SPM.swd;
            ind=strfind(folder, 'IRB13278');
            folder=fullfile(CommonEnvMethods.getCOHPath(),'data', folder(ind:length(folder)));
            folder=strrep(folder,'\',filesep);
            tbl={};
            fnameSPM=fullfile(folder,'SPM.mat');
            cons=SPM.xCon;
            len=length(cons);
            for c=1:len
                con=cons(c);
                mask=fullfile(folder,con.Vspm.fname);
                if(c==1)
%                     newInstance=true;
%                     xjViewh=xjViewHandler.get_xjView(newInstance);
                    objt=xjViewHandler(mask);
                    objt.selectPositiveOnly();
                    xjViewFigHandle=objt.getXjviewFigHandle();
                end
                msg=['Getting stats summary table ' con.Vspm.fname ' at ' datestr(clock)];
                cprintf('blue', '%s\n', msg);
                obj=SPM_statFileHandler(xjViewFigHandle,mask,fnameSPM,1,0.001,1);
                tblt=obj.getOutputTableText();
                tblt=CommonMethods.applyingLineIndentationToTable(tblt,2,dil);
                tblt=CommonMethods.insertRowToTable(tblt,[' ' dil con.name dil con.Vspm.fname],1,dil);
                if(c==1)
                    tblt=CommonMethods.insertRowToTable(tblt,fnameSPM,1,dil);
                end
                SPM_statFileHandler.exportTableText(tblt,fid);
            end
            if(closefile)
                fclose(fid);
            end
            CommonMethods.closeAll('xjView');
            delete(xjViewFigHandle);
        end
        
        function exportTableText(tbl,fid)
            for i=1:length(tbl)
                fprintf(fid,'%s\n',tbl{i});
            end
        end
        
        function exportSPMTable_fid(fnameSPM, fid,xjView)
            SPM_statFileHandler.exportSPMSummaryTable(fnameSPM, fid, xjView);
        end
        function exportSPMTable(fnameSPM, fname)
            fid=fopen(fname,'wt');
            SPM_statFileHandler.exportSPMTable_fid(fnameSPM,fid);
            fclose(fid);
        end
        
        function len=exportSPMTable_logFile_par(fname, checkProgress)
            %'_par' do it with parfor
            %fname: the lob file name
            %%%There will be problem if the stat (T/F) maps have no
            %%%superathreshold voxels. 
            %
            clc; 
            if(~exist('fname','var'))
                fname='';
%               %% the log file (the output of group analysis)
            end
            if(isempty(fname))
                fname=fullfile('/home1/data/IRB13278/data/Resting/group/analysisLogs','Rest_ALFF_fALFF_ReHo_FreqSpecific_TIV_GroupAnalysis_Regressions_2017-10-27_log.tsv');
%               %% the log file (the output of group analysis)
            end
            
            if(~exist('checkProgress','var'))
                checkProgress=false;
            end
            tab=sprintf('%s\t','');
            fid=fopen(fname);

            lines = textscan(fid, '%s', 'Delimiter','\n');
            lines=lines{1};            
            fclose(fid);
            
            len=length(lines);
            fnames=cell(1,len);
            funNames=cell(1,len);
            for i=1:len
                strs=strsplit(lines{i},tab);
                if(length(strs) <2)
                    continue;
                end
                funNames{i}=strs{1};
                fnames{i}=strs{2};
            end
            
            [folder,name,ext]=fileparts(fname);
            flog=fullfile(folder,['SPM_Summaries_' name '_error.log']);
            
            fnames=fnames(cellfun(@(x)~isempty(strfind(x,filesep)), fnames));
            idx=cellfun(@(x)~SPM_statFileHandler.summarized(x),fnames);
            fnames_o=fnames;
            fnames=fnames(idx);
            len=length(fnames);
            if(checkProgress)
                return;
            end
%            parfor i=1:len
            for i=1:1
                fnameSPM=fullfile(fnames{i}, 'SPM.mat');
                cd(fnames{i});
                cprintf('blue','%s\n',['Summarizing ' num2str(i) '/' num2str(len) ' ' fnameSPM]);
                
                try
                    SPM_statFileHandler.exportSPMSummaryTable(fnameSPM);
                catch err
                    fidLog=fopen(flog,'at');
                    while(fidLog<0)
                        %flog is locked by other threads. Wait for 2
                        %seconds and open again.
                        pause(2);
                        fidLog=fopen(flog,'at');
                    end
                    msg=['Error with summarizing ' fnameSPM];
                    cprintf('red', '%s\n', msg);
                    fprintf(fidLog,'%s\n',msg);
                    fclose(fidLog);
                end
            end
            [folder,name,ext]=fileparts(fname);
            fname1=fullfile(folder,['SPM_Summaries_' name ext]);
            fid1=fopen(fname1, 'wt');
            fnames=fnames_o;
            for f=1:length(fnames)
                sname=fullfile(fnames{f},'spmSummary.tsv');
                fidt=fopen(sname);
                lines = textscan(fidt, '%s', 'Delimiter','\n');
                lines=lines{1};
                for l=1:length(lines)
                    fprintf(fid1,'%s\n',lines{l});
                end
                fclose(fidt);
            end
            fclose(fid1);
        end     
        
        function summarizeSPM(folder,fname)
            if(~exist('folder','var'))
                folder=pwd;
%                 root='/home1/data/IRB13278/data/VBM/group_gmc/anova_repeatedMeasure_ageTICV';
%                 folder={
% %                    fullfile(root,'gmcneuropsychGMCChangechangeCorrelations') ...                    
%                     fullfile(root,'gmcneuropsychCorrelationstp1') ...
%                     fullfile(root,'gmcneuropsychCorrelations_change2change') ...
%  %                   fullfile(root,'gmcneuropsychchangeCorrelations') ...
%                     }
            end
            if(~exist('fname','var'))
 %               fname=fullfile('/home1/data/IRB13278/data/VBM/group_gmc/SummaryFiles',['SPM_Summaries_' datestr(clock,'yy-mm-dd') '.tsv']);
                root=pwd;
                fname=fullfile(root,['SPM_Summaries_' datestr(clock,'yy-mm-dd') '_ALFF_fALFF_ReHo_all.tsv']);
            end
            if ~iscell(folder)
                folder={folder};
            end
            fnames={};
            for i=1:length(folder)
                fnames=horzcat(fnames, CommonMethods.getFileList_recursive(folder{i},'.mat'));
            end
            [folders,names,exts]=cellfun(@(x)fileparts(x),fnames,'UniformOutput',false);
            selection=cellfun(@(x)strcmp(x,'SPM'),names);
            folders=folders(selection);
%             maskName='mask_FA_gt_point1';
%             selection=cellfun(@(x)~isempty(strfind(x,maskName)),folders);
%            folders=folders(selection);
            SPM_statFileHandler.exportSPMTable_fileList_par(fname,folders);
        end
        
        function len=exportSPMTable_fileList_par(fnameOutput,SPMFolders)
            %'_par' do it with parfor
            %
            %SPMFolders: cellarray list of folders containing SPM.mat
            %fnameOutput: the output file
            clc; 
            if(~exist('fnameOutput','var'))
                fnameOutput='';
            end
            if(~exist('SPMFolders','var'))
                SPMFolders={pwd};
            end
            idx=cellfun(@(x)~SPM_statFileHandler.summarized(x),SPMFolders);
            SPMFolders_o=SPMFolders;
            SPMFolders=SPMFolders(idx);
            len=length(SPMFolders);
%            parfor i=1:len
            for i=1:len
                fnameSPM=fullfile(SPMFolders{i}, 'SPM.mat');
                cd(SPMFolders{i});
                cprintf('blue','%s\n',['Summarizing ' num2str(i) '/' num2str(len) ' ' fnameSPM]);
                try
                    SPM_statFileHandler.exportSPMSummaryTable(fnameSPM);
                catch err
                    flog=fullfile(pwd,'SPM_summaryErr.log');
                    fidLog=fopen(flog,'at');
                    while(fidLog<0)
                        %flog is locked by other threads. Wait for 2
                        %seconds and open again.
                        pause(2);
                        fidLog=fopen(flog,'at');
                    end
                    msg=['Error with summarizing ' fnameSPM ':' err.message];
                    cprintf('red', '%s\n', msg);
                    fprintf(fidLog,'%s\n',msg);
                    fclose(fidLog);
                end
            end
            if(isempty(fnameOutput))
                return;
            end
            fid1=fopen(fnameOutput, 'wt');
            for f=1:length(SPMFolders_o)
                sname=fullfile(SPMFolders_o{f},'spmSummary.tsv');
                fidt=fopen(sname);
                lines = textscan(fidt, '%s', 'Delimiter','\n');
                lines=lines{1};
                for l=1:length(lines)
                    fprintf(fid1,'%s\n',lines{l});
                end
                fclose(fidt);
            end
            fclose(fid1);
        end     
        
        function IS=summarized(folder)
            IS=false;
            date='02/22/2018';
            dtn=datenum(date,'mm/dd/yyyy');
%            return;
            fname=fullfile(folder,'spmSummary.tsv');
            fnameSPM=fullfile(folder,'SPM.mat');
            if(~exist(fname,'file')||~exist(fnameSPM,'file'))
                return;
            end
            lst=dir(fname);
            if(lst.bytes==0)
                return;
            end
            if(lst.datenum < dtn)
                return;
            end
            lstSPM=dir(fnameSPM);
            if(lstSPM.datenum>lst.datenum)
                return;
            end
            tmp=load(fnameSPM);
            noCon=length(tmp.SPM.xCon);
            lines=CommonMethods.getTextLines(fname);
            idx=cellfun(@(x)~isempty(strfind(x,'spmT_'))||~isempty(strfind(x,'spmF_')),lines);
            if(nnz(idx)<noCon)
                return;
            end
            %checking whether the last contrast was completed
            inds=find(idx);
            iI=inds(noCon);
            len=length(lines);
            if(len<iI+2)
                return;
            end
            
            %checking whether the last line is incomplete
            line=lines{iI};
            
            strsI=strsplit(line,CommonMethods.tab);
            strsF=strsplit(lines{len},CommonMethods.tab);
            if(length(strsI)>length(strsF))
                return;
            end
            IS=true;
        end
    end
    
    
end