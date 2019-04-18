classdef spmStatSummaryHandler < handle
    %SPMSTATSUMMARYHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        statSummaryNodes;%summary node structure array
        dataRoot;%root directory of SPM stat files: the common path among all directory
        subfolders_u;%unique subfolders
        conImgs_sf;%con images in each subfolder
        p_cluster;%cluster level corrected p value threshold
        p_v_FWE; %voxel level FWE correction p value threshold
        p_v_FDR; %voxel level FDR correction p value threshold
        p_cluster_m;%cluster level corrected p value threshold after multiple comparision adjustment
        p_v_FWE_m; %voxel level FWE correction p value threshold after multiple comparision adjustment
        p_v_FDR_m; %voxel level FDR correction p value threshold after multiple comparision adjustment
        p_c;%cluster level p values for multiple comparison correction
        p_vFWE;%voxel level FWE corrected p values for multiple comparison correction
        p_vFDR;%voxel level FDR corrected p values for multiple comparison correction
        summaryInds;%indexes statSummaryNodes
        imgNames;%image names
        columnIDs;%image names for this version
        rowIDs;%combination of subfolder and conImg
        selectionNode;%
        sigRowsOnly%outputhe summarys only if there is at least one significant results in a rwo..
        resultsOutput;%a structure conaing parameters related with the output table (cell array of string).
        summaryFile;%the SPM summary file
        multipleComparisonCorrOption; %a structure containing information need for multiple comparison adjustment
        version;
    end
    
    methods
        function obj=spmStatSummaryHandler(summaryFile,p_cluster, p_v_FWE, p_v_FDR,multipleComparisonCorrOption)
            %mcro: multiple comparison corection option
            version=1;
            obj.statSummaryNodes={};
            if(~exist('p_cluster','var'))
                p_cluster=[];
            end
            if(~exist('p_v_FWE','var'))
                p_v_FWE=[];
            end
            if(~exist('p_v_FDR','var'))
                p_v_FDR=[];
            end
            if(isempty(p_cluster))
                p_cluster=1.1;
            end
            if(isempty(p_v_FWE))
                p_v_FWE=1.1;
            end
            if(isempty(p_v_FDR))
                p_v_FDR=1.1;
            end
            if(~exist('summaryFile','var'))
                summaryFile='';
            end
            if(~exist('multipleComparisonCorrOption','var'))
                multipleComparisonCorrOption={};
            end
            
            if(isempty(summaryFile))
                summaryFile=fullfile(COH1_CommonMethods.getDatadir(), 'Resting', 'group', 'analysisLogs','SPM_Summaries_Rest_ALFF_fALFF_ReHo_GroupAnalysis_all_21-Sep-2017-12-50-02_log.tsv');
            end
            
            [folder, name, ~]=fileparts(summaryFile);
            try
                savedFile=fullfile(folder,['savedStatSummaryNodes_' name '.mat']);
%                delete(savedFile);
                if(CommonMethods.IsANewerFile(savedFile,summaryFile))
                    delete(savedFile);
                end
                load(fullfile(folder,['savedStatSummaryNodes_' name '.mat']));
                if(~isfield(obj,'version'))
                    obj.version=1;
                end
                if(version > obj.version)
                    error('The version of the saved object is not supported');
                    %Throwing an error to force it to rebuild the obj.
                end
                newCutoff=false;
                if(p_cluster~=obj.p_cluster)
                    obj.p_cluster=p_cluster;
                    newCutoff=true;
                end
                if(p_v_FWE~=obj.p_v_FWE)
                    obj.p_v_FWE=p_v_FWE;
                    newCutoff=true;
                end
                if(p_v_FDR~=obj.p_v_FDR)
                    obj.p_v_FDR=p_v_FDR;
                    newCutoff=true;
                end
                if(newCutoff)
                    % to be implemented
                end
                obj.summaryFile=summaryFile;
%                 obj.p_c=cellfun(@(x)min(x.p_c),obj.statSummaryNodes);
%                 obj.p_vFWE=cellfun(@(x)min(x.p_vFWE),obj.statSummaryNodes);
%                 obj.p_vFDR=cellfun(@(x)min(x.p_vFDR),obj.statSummaryNodes);
%                 save(fullfile(folder,['savedStatSummaryNodes_' name '.mat']),'obj');
%                obj.tabulateSummaryNodes();
            catch err
%                 obj.p_cluster=p_cluster;
%                 obj.p_v_FWE=p_v_FWE;
%                 obj.p_v_FDR=p_v_FDR;
%                 obj.p_cluster_m=p_cluster;
%                 obj.p_v_FWE_m=p_v_FWE;
%                 obj.p_v_FDR_m=p_v_FDR;
                %Keep everything when build it for the first time.
                obj.p_cluster=1.1;
                obj.p_v_FWE=1.1;
                obj.p_v_FDR=1.1;
                obj.p_cluster_m=1.1;
                obj.p_v_FWE_m=1.1;
                obj.p_v_FDR_m=1.1;
                msg='building summary nodes';
                cprintf('blue','%s\n',msg);
                obj.buildSummaryNodes(summaryFile);
                msg='tabulating summary nodes';
                cprintf('blue','%s\n',msg);
                obj.tabulateSummaryNodes();
                msg='building new summary node fields';
                cprintf('blue','%s\n',msg);
                obj.buildNewSummaryNodeFields();
                msg='checing SNs';
                cprintf('blue','%s\n',msg);
                obj.checkingSNs();
                msg='building Selection';
                cprintf('blue','%s\n',msg);
                obj.buildSelection();
                msg='saveing results';
                cprintf('blue','%s\n',msg);
                obj.p_c=cellfun(@(x)min(x.p_c),obj.statSummaryNodes);
                obj.p_vFWE=cellfun(@(x)min(x.p_vFWE),obj.statSummaryNodes);
                obj.p_vFDR=cellfun(@(x)min(x.p_vFDR),obj.statSummaryNodes);
                save(fullfile(folder,['savedStatSummaryNodes_' name '.mat']),'obj');
            end
            
            obj.p_cluster=p_cluster;
            obj.p_v_FWE=p_v_FWE;
            obj.p_v_FDR=p_v_FDR;
            obj.p_cluster_m=p_cluster;
            obj.p_v_FWE_m=p_v_FWE;
            obj.p_v_FDR_m=p_v_FDR;      
            obj.summaryFile=summaryFile;
            obj.multipleComparisonCorrOption=multipleComparisonCorrOption;
            fprintf(1,'Checking SNs\n');
            if(~isempty(obj.multipleComparisonCorrOption))
                obj.adjustMultipleComarison();
            else
                obj.p_cluster_m=p_cluster;
                obj.p_v_FWE_m=p_v_FWE;
                obj.p_v_FDR_m=p_v_FDR;
            end
            obj.checkingSNs();
            fprintf(1,'Buiding Table Data\n');
            obj.buildResultOutput();
            obj.resultsOutput.summaryFile=summaryFile;
        end
        
        function adjustMultipleComarison(obj)
            options=obj.multipleComparisonCorrOption;
            rows=length(obj.rowIDs);
            cols=length(obj.columnIDs);
            if(~isfield(obj.multipleComparisonCorrOption,'rowSelection'))
                rowSelection=true(1,rows);
            else
                if(length(obj.multipleComparisonCorrOption.rowSelection)==rows)
                    rowSelection=obj.multipleComparisonCorrOption.rowSelection;
                else
                    rowSelection=ture(1,rows);
                end
            end
            if(~isfield(obj.multipleComparisonCorrOption,'colSelection'))
                colSelection=true(1,cols);
            else                
                if(length(obj.multipleComparisonCorrOption.colSelection)==cols)
                    colSelection=obj.multipleComparisonCorrOption.colSelection;
                else
                    colSelection=ture(1,cols);
                end
            end
            inds=find(obj.summaryInds(rowSelection,colSelection));
            option=options.correctionOption;
            if(strcmp(option, 'Bonferroni-Holm'))
                if(obj.p_cluster <1)
                    obj.p_cluster_m=bonf_holm(obj.p_c(inds), obj.p_cluster);
                end
                if(obj.p_v_FWE <1)
                    obj.p_v_FWE_m=bonf_holm(obj.p_vFWE(inds), obj.p_v_FWE);
                end
                if(obj.p_v_FDR < 1)
                    obj.p_v_FDR_m=bonf_holm(obj.p_vFDR(inds),obj.p_v_FDR);
                end
            elseif(strcmp(option,'FDR'))
                report='no';
                method='pdep';
                if(obj.p_cluster <1)
                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(obj.p_c(inds),obj.p_cluster,method,report);
                    obj.p_cluster_m=crit_p;
                end
                if(obj.p_v_FWE <1)
                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(obj.p_vFWE(inds),obj.p_v_FWE,method,report);
                    obj.p_v_FWE_m=crit_p;
                    
                end
                if(obj.p_v_FDR <1)
                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(obj.p_vFDR(inds),obj.p_v_FDR,method,report);
                    obj.p_v_FDR_m=crit_p;
                end
            end
        end
        
        function tbl=getComparisonTable(obj)
            tbl=obj.resultsOutput.table;
            
        end
        
        function strs=strsplit(str, del)
            len1=length(del);
            len2=length(str);
            idx=strfind(str,del);
            iIs=[1 [idx+len1]];
            iFs=[[idx-1] len2];
            len=length(idx)+1;
            strs=arrayfun(@(x)str(iIs(x):iFs(x)),1:len,'UniformOutpt',false);
        end
        function buildSummaryNodes(obj,fname)
            %fname is a output file of spm stat summary
            tab=sprintf('%s\t','');
            fid=fopen(fname);
            lines = textscan(fid, '%s', 'Delimiter','\n');
            %%%textscan ignores the giginning white-space from each line
            lines=lines{1};
            lines=cellfun(@(x)strrep(x,'C:/Users/tjin/City of Hope', '/home1'),lines,'UniformOutput', false);
            fclose(fid);
            len=length(lines);
            fnames=cell(1,len);
            funNames=cell(1,len);
            inds=[];%the index of the lines containing 'SPM.mat' file name
            for i=1:len
                strs=strsplit(lines{i},tab);
                str=strs{1};
                [~,name,ext]=fileparts(str);
                if(strcmp([name ext],'SPM.mat'))
                    inds(end+1)=i;
                end
            end
            len=length(inds);
            lent=length(lines);
            
            for i0=1:len
                if(mod(i0,100)==1)
                    msg=['building the ' num2str(i0) 'th summary node (total: ' num2str(len) ')'];
                    cprintf('blue', '%s\n',msg);
                end
                if(i0==868)
                    a=1;
                end
                lI=inds(i0);
                if(i0<len)
                    lF=inds(i0+1)-1;
                else
                    lF=lent;
                end
                obj.addSummaryNode(lines,lI,lF);%the lines about a single SPM.mat
            end
        end
        function addSummaryNode(obj,lines,lI0,lF0)
            %building summary nodes and add to obj.statSummaryNodes
            tab=sprintf('%s\t','');
            fnameSPM=lines{lI0};
            idx=[];
            for l=lI0+1:lF0
                strs=strsplit(lines{l},tab);
                if(length(strs) <2)
                    continue;
                end
                if(~isempty(strfind(strs{2},'.nii')))
                    idx(end+1)=l;
                end
            end
            len=length(idx);
            for i=1:len
                lI=idx(i);
                if(i<len)
                    lF=idx(i+1)-1;
                else
                    lF=lF0;
                end
                obj.statSummaryNodes{end+1}=obj.buildSummaryNode(lines,lI,lF,fnameSPM);
            end
        end
        
        function sn=buildSummaryNode(obj,lines,lI,lF,fnameSPM)
            %this function build summary node
            tab=CommonMethods.tab;
            sn.fnameSPM=fnameSPM;
            strs=strsplit(fnameSPM,filesep);
            [~,pos]=ismember('group',strs);
            sn.imgName=strs{pos+1};
            sn.subfolder=CommonMethods.cell2str(strs(pos+2:length(strs)-1),filesep);
            
            strs=strsplit(lines{lI},tab);
            sn.conDisp=strs{1};
            sn.conImg=strs{2};
            sn.header={CommonMethods.replaceEmptyCells(strsplit(lines{lI+1},tab), ' ');CommonMethods.replaceEmptyCells(strsplit(lines{lI+2}), ' ')};
            %selecting significant clusters
            idx=lI+3:lF;
            p_c=arrayfun(@(x)obj.getp_cluster(lines{x}), idx);
            p_vFWE=arrayfun(@(x)obj.getp_vFWE(lines{x}), idx);
            p_vFDR=arrayfun(@(x)obj.getp_vFDR(lines{x}), idx);
%             selection=p_c <= obj.p_cluster;
%             idx=idx(selection);
            sn.data={};
            if(isempty(idx))
                sn.p_c=1.1;
                sn.p_vFWE=1.1;
                sn.p_vFDR=1.1;
            else
                sn.p_c=p_c;
                sn.p_vFWE=p_vFWE;
                sn.p_vFDR=p_vFDR;
            end
            len=length(idx);
            for i=1:len
                sn.data{end+1}=strsplit(lines{idx(i)},tab);
            end
        end
        
        function p=getp_cluster(obj,dataline)
            tab=sprintf('%s\t','');
            strs=strsplit(dataline,tab);
            try
                p=str2num(strs{3});
            catch err
                p=1.1;%caused by corrupted summary file. will be able to tell by p values
            end
        end
        
        function p=getp_vFWE(obj,dataline)
            tab=sprintf('%s\t','');
            strs=strsplit(dataline,tab);
            try
                p=str2num(strs{4});
            catch err
                p=1.1;%caused by corrupted summary file. will be able to tell by p values
            end
        end
        
        function p=getp_vFDR(obj,dataline)
            tab=sprintf('%s\t','');
            strs=strsplit(dataline,tab);
            try
                p=str2num(strs{5});
            catch err
                p=1.1;%caused by corrupted summary file. will be able to tell by p values
            end
        end
        
        function sortSubfolder(obj)
            %this function sorts subfolders according to the regression
            %types, rather than the alpbetical
            corrTypes={'corr_tp1' 'corr_Img2ScoreChange' 'ImgChange2Score_corr' 'corr_change2change'};
            corrTypes_order={'1corr_tp1' '2corr_Img2ScoreChange' '3ImgChange2Score_corr' '4corr_change2change'};;
            %            coorTypes_g={'Img2Score_correlation_tp1' 'Image2ScoreChange_correation' 'ImageChange2Score_correlation' 'ImageChange2ScoreChange_correlation'};
            subfolders_u=obj.subfolders_u;
            len=length(subfolders_u);
            sfu=cell(1,len);
            %            sfg=cell(1,len);
            for i=1:len
                subfolder=subfolders_u{i};
                for j=1:length(corrTypes)
                    if(isempty(strfind(subfolder,corrTypes{j})))
                        continue;
                    end
                    subfolder=strrep(subfolder,corrTypes{j},corrTypes_order{j});
                end
                sfu{i}=subfolder;
            end
            
            [~, io]=sort(sfu);
            obj.subfolders_u=subfolders_u(io);
            %             obj.conImgs_sf=conImgs_sf(io);
            %             obj.conDisps_sf=conDisps_sf(io);
        end
        function tabulateSummaryNodes(obj)
            subfolders=cellfun(@(x)x.subfolder,obj.statSummaryNodes, 'UniformOutput', false);
            conImgs=cellfun(@(x)x.conImg,obj.statSummaryNodes, 'UniformOutput', false);
            conDisps=cellfun(@(x)x.conDisp,obj.statSummaryNodes, 'UniformOutput', false);
            imgNames=cellfun(@(x)x.imgName,obj.statSummaryNodes, 'UniformOutput', false);
            obj.subfolders_u=unique(subfolders);
            len=length(obj.subfolders_u);
            obj.conImgs_sf=cell(1,len);
            
            obj.imgNames=unique(imgNames);
            obj.sortSubfolder();
            obj.columnIDs=obj.imgNames;
            obj.rowIDs={};
            for s=1:length(obj.subfolders_u)
                subfolder=obj.subfolders_u{s};
                idx=ismember(subfolders,subfolder);
                conImgst=unique(conImgs(idx));
                obj.conImgs_sf{s}=conImgst;
                for c=1:length(conImgst)
                    conImg=conImgst{c};
                    obj.rowIDs{end+1}={subfolder conImg};
                end
            end
            
            tmpSN=obj.statSummaryNodes{1};%a template summary node for making dummy nodes
            tmpSN.data={};
            
            rows=length(obj.rowIDs);
            cols=length(obj.columnIDs);
            summaryInds=zeros(rows,cols);
            for r=1:rows
                rc=obj.rowIDs{r};
                subfolder=rc{1};
                conImg=rc{2};
                for c=1:cols
                    imgName=obj.columnIDs{c};
                    idx=find(ismember(subfolders, subfolder).*ismember(conImgs, conImg).*ismember(imgNames, imgName));
%                    if(length(idx)==1)
                    if(~isempty(idx))
                        summaryInds(r,c)=idx(1);
                    else %adding dummy node
                        summaryInds(r,c)=length(obj.statSummaryNodes)+1;
                        sn=tmpSN;
                        sn.fnameSPM=strrep(sn.fnameSPM,sn.subfolder,subfolder);
                        sn.conImg=conImg;
                        sn.subfolder=subfolder;
                        sn.conDisp='dummy node';
                        obj.statSummaryNodes{end+1}=sn;
                    end
                end
            end
            obj.summaryInds=summaryInds;
        end
        function buildSelection(obj)
%             folder=fullfile(COH1_CommonMethods.getDatadir(), 'Resting', 'group', 'analysisLogs');
%             fname_inclusion1=fullfile(folder,'AnalysisTypes_Inclussion.xlsx');
%             fname_inclusion2=fullfile(folder,'AnalysisTypes_Inclussion_1sampleT.xlsx');	
%             [~,~,raw]=xlsread(fname_inclusion1);
%             exclusionColumn=find(cellfun(@(x)strcmp(x,'exclude'),raw(1,:)));
%             exclusion1=cellfun(@(x)strcmp(x,'a'), raw(2:size(raw,1),exclusionColumn));
%             [~,~,raw]=xlsread(fname_inclusion2);
%             exclusionColumn=find(cellfun(@(x)strcmp(x,'exclude'),raw(1,:)));
%             exclusion2=cellfun(@(x)strcmp(x,'a'), raw(2:size(raw,1),exclusionColumn));
%             exclusion=(exclusion1.*exclusion2)~=0;
            exclusion=false(1,length(obj.rowIDs));
            obj.selectionNode.exclusion=exclusion;
        end
        
        function checkingSNs(obj)
            %recomputing cluster selection 
            snst=obj.statSummaryNodes;
            
            cs=cellfun(@(x)length(x.data),snst);
            len=length(snst);
            sns=cell(1,len);
            for i=1:len
                sn=snst{i};   
                data=sn.data;
                if(~isempty(data))
                    cs=cellfun(@(x)str2num(x{2}),data);
                    cSize=cs(1);
                    ts=cellfun(@(x)str2num(x{6}),data);                    
                    t=ts(1);
                    
                    if(cSize>10000&&t<4)
                        sn.data={};
                    else
                        cSizes=cellfun(@(x)str2num(x{2}),data);
                        dd=diff([cSizes 0]);
                        inds=find(dd);%last index of each cluster
                        pcs=cellfun(@(x)str2num(x{3}),data);
                        pFWEs=cellfun(@(x)str2num(x{4}),data);
                        pFDRs=cellfun(@(x)str2num(x{5}),data);
                        selection=(pcs<=obj.p_cluster_m).*(pFWEs<=obj.p_v_FWE_m).*(pFDRs<=obj.p_v_FDR_m);
                        selection=selection>0;
                        %selecting the entire cluster if anay local peaks
                        %are selected
                        iI=1;
                        for t=1:length(inds)
                            idx=iI:inds(t);
                            iI=inds(t)+1;
                            if(any(selection(idx)))
                                selection(idx)=true;
                            end
                        end
                        sn.data=data(selection);
                        selection_c=false(size(selection));
                        selection_c(inds)=true;
                        sn.clusterSizes=cSizes((selection.*selection_c)>0);
                    end
                end
                sns{i}=sn;
            end
            obj.statSummaryNodes=sns;
        end
        function exportResultsComparison_list(obj)
            folder='/home1/data/IRB13278/data/Resting/group/analysisLogs';
            imgNames=obj.imgNames;
            
            imgLine=imgNames{1};
            imgLine_blank=' ';
            tab=CommonMethods.tab;
            for i=2:length(imgNames)
                imgLine=[imgLine tab imgNames{i}];
                imgLine_blank=[imgLine_blank tab ' '];
            end
            
            fo=fopen(fullfile(folder,'AnalysisTypes.tsv'),'wt');
            tab=CommonMethods.tab;
            line=['subfolder' tab 'con img' tab 'con description' tab imgLine];
            fprintf(fo,'%s\n',line);
            
            colIDs=obj.columnIDs;
            rowIDs=obj.rowIDs;
            len=length(rowIDs);
            exclusion=obj.selectionNode.exclusion;
            for i=1:len
                subfolder=rowIDs{i}{1};
                conImg=rowIDs{i}{2};
                sn=obj.statSummaryNodes{obj.summaryInds(i,1)};
                if(exclusion(i))
                    continue;
                end
                line=[subfolder tab conImg tab sn.conDisp tab imgLine_blank];
                fprintf(fo,'%s\n',line);
            end
            fclose(fo);
        end
        function exportResultsComparison(obj)
            folder='/home1/data/IRB13278/data/Resting/group/analysisLogs';
            imgNames=obj.imgNames;
            
            imgLine=imgNames{1};
            imgLine_blank=' ';
            tab=CommonMethods.tab;
            
            fo=fopen(fullfile(folder,'AnalysisTypes.tsv'),'wt');
            tab=CommonMethods.tab;
            
            colIDs=obj.columnIDs;
            rowIDs=obj.rowIDs;
            rows=length(rowIDs);
            cols=length(colIDs);
            tw=10;%table with for each summary node
            
            lineSegs=cellfun(@(x)CommonMethods.cell2tableline(x,1,tw,tab,' '), imgNames, 'UniformOutput', false);
            imgLine=CommonMethods.cell2str(lineSegs, tab);
            
            line=['subfolder' tab 'con img' tab 'con description' tab imgLine];
            fprintf(fo,'%s\n',line);
            
            exclusion=obj.selectionNode.exclusion;
            for r=1:rows
                subfolder=rowIDs{i}{1};
                conImg=rowIDs{i}{2};
                sn=obj.statSummaryNodes{obj.summaryInds(i,1)};
                if(exclusion(i))
                    continue;
                end
                line=[subfolder tab conImg tab sn.conDisp tab imgLine_blank];
                fprintf(fo,'%s\n',line);
            end
            fclose(fo);
        end
        
        function [tbl,obj]=getResultsComparison(obj)
            tbl=obj;
        end
        
        
        function tbl=buildResultOutput(obj)
            %To do: handle multiple comparison, the options will be passed
            %in as input argument to the constuctor. 
            
            colIDs=obj.columnIDs;
            imgNames=obj.imgNames;
            rowIDs=obj.rowIDs;
            
            if(~isfield(obj,'selectionNode'))
                obj.selectionNode=[];
            end
            
            if(~isfield(obj.selectionNode,'exclusion'))
                obj.selectionNode.exclusion=false(1,size(obj.summaryInds,1));
            end
            exclusion=obj.selectionNode.exclusion;
            inds_include=find(~exclusion);
            rows=length(inds_include);
            cols=length(colIDs);
            tw=10;%table with for each summary node
            obj.resultsOutput.snTableWidth=tw;
            tblResults=cell(rows,cols);
            tblSizes=zeros(rows,1);
            numImgs=length(imgNames);
            tbl_o=cell(rows,numImgs);
            rowNames=cell(rows,1);
            %the number of lines in result tables
            for r=1:rows
                ri=inds_include(r);
                sns=arrayfun(@(x)obj.statSummaryNodes{obj.summaryInds(ri,x)}, 1:cols, 'UniformOutput',false);
                tblt=obj.sns2str(sns,tw);
                tblResults(r,:)=tblt;
                tblSizes(r)=size(tblt{1},1);                
                rowNames{r}=strrep(sns{1}.conDisp, 'R_ALFF', 'rsfMRI');
            end
            summaryNodeLookupTable=zeros(sum(tblSizes),cols*tw);
            significant=false(sum(tblSizes),cols*tw);
            tbl=cell(sum(tblSizes),cols*tw);
            summaryNodeLookupTable_o=zeros(rows,numImgs);%for the overview table
            significant_o=false(rows,numImgs);
            preRow=0;
            
            for r=1:rows
                ri=inds_include(r);
                n=tblSizes(r);
                preCol=0;
%                tbl_o(r,1:2)={rowIDs{r}{1},rowIDs{r}{2}};
                for c=1:cols
                    tblt=tblResults{r,c};
                    ind=obj.summaryInds(ri,c);
                    summaryNodeLookupTable_o(r,c)=ind;
                    sn=obj.statSummaryNodes{ind};
                    
                    if(isempty(sn.data))
                        if(~isempty(strfind(sn.conDisp,'dummy')))
                            
                            tbl_o{r,c}='dummy';
                            sig=false;
                        else
                            tbl_o{r,c}=' ';
                            sig=false;
                        end
                    else
                        cs=[sn.clusterSizes 0];
                        idx=find(diff(cs));
                        cs=cs(idx);
                        tbl_o{r,c}=[num2str(sum(cs)) ' (' num2str(length(cs)) ')'];
                        significant_o(r,c)=true;
                        sig=true;
                    end
                    for r1=1:n
                        curRow=preRow+r1;
                        for c1=1:tw
                            curCol=preCol+c1;
                            summaryNodeLookupTable(curRow,curCol)=ind;
                            tbl{curRow,curCol}=tblt{r1,c1};
                            significant(curRow,curCol)=sig;
                        end
                    end
                    preCol=preCol+tw;
                end
                preRow=preRow+n;
            end
            
            obj.resultsOutput.table_detailed.data=tbl;
            obj.resultsOutput.table_detailed.type='detailed';
            obj.resultsOutput.table_detailed.colNames=arrayfun(@(x) CommonMethods.getDefaultColumnNames(x), 1:cols*tw, 'UniformOutput', false);
            obj.resultsOutput.table_detailed.rowNames=rowNames;
            obj.resultsOutput.table_detailed.summaryNodeLookupTable=summaryNodeLookupTable;
            obj.resultsOutput.table_detailed.tableSizes=tblSizes;
            obj.resultsOutput.table_detailed.significant=significant;
            obj.resultsOutput.table_detailed.p_cluster=obj.p_cluster;
            obj.resultsOutput.table_detailed.p_v_FWE=obj.p_v_FWE;
            obj.resultsOutput.table_detailed.p_v_FDR=obj.p_v_FDR;
            obj.resultsOutput.table_detailed.p_cluster_m=obj.p_cluster_m;
            obj.resultsOutput.table_detailed.p_v_FWE_m=obj.p_v_FWE_m;
            obj.resultsOutput.table_detailed.p_v_FDR_m=obj.p_v_FDR_m;
            obj.resultsOutput.table_detailed.summaryFile=obj.summaryFile;
                      
            obj.resultsOutput.table_overview.significant=significant_o;
            obj.resultsOutput.table_overview.type='overview';
            obj.resultsOutput.table_overview.data=tbl_o;
 %           obj.resultsOutput.table_overview.header=horzcat({'subflder' 'con discription' 'con img'}, imgNames);
            obj.resultsOutput.table_overview.colNames=imgNames;
            obj.resultsOutput.table_overview.rowNames=rowNames;
            obj.resultsOutput.table_overview.summaryNodeLookupTable=summaryNodeLookupTable_o;
            obj.resultsOutput.handleTableSelection=@obj.handleTableSelection;
            obj.resultsOutput.table_overview.p_cluster=obj.p_cluster;
            obj.resultsOutput.table_overview.p_v_FWE=obj.p_v_FWE;
            obj.resultsOutput.table_overview.p_v_FDR=obj.p_v_FDR;
            obj.resultsOutput.table_overview.p_cluster_m=obj.p_cluster_m;
            obj.resultsOutput.table_overview.p_v_FWE_m=obj.p_v_FWE_m;
            obj.resultsOutput.table_overview.p_v_FDR_m=obj.p_v_FDR_m;
            obj.resultsOutput.table_detailed.summaryFile=obj.summaryFile;
            obj.resultsOutput.table_overview.summaryFile=obj.summaryFile;
        end
        function buildNewSummaryNodeFields(obj,inds)
            if(~exist('ind','var'))
                inds=1:length(obj.statSummaryNodes);
            end
            for i=1:length(inds)
                ind=inds(i);
                sn=obj.statSummaryNodes{ind};
                sn.clusterSizes=cellfun(@(x)str2num(x{2}),sn.data);
                strs=strsplit(sn.subfolder, filesep);
                sn.analysisType=strs{1};
                obj.statSummaryNodes{ind}=sn;
            end
        end
        function tbl=sns2str(obj,sns,tw)
            %returns a cell array (ml,tw) cell array. ml is the maximum of
            %the number of rows needed for sns.
            if(~exist('tw','var'))
                tw=obj.resultsOutput.snTableWidth;
            end
            nls=cellfun(@(x)obj.getSummaryNodeSize(x),sns);
            ml=max(nls);
            tbl=cellfun(@(x)obj.sn2str(x,ml,tw),sns,'UniformOutput',false);
        end
        function tbl=sn2str(~,sn,ml,tw)
            %returns a (ml,tw) cellarray
            tbl=cell(ml,tw);
            tbl(1,:)=(CommonMethods.cell2tablecells({sn.conDisp sn.conImg ' ' sn.imgName},1,tw,' '));
            try
                tbl(2,:)=(CommonMethods.cell2tablecells(sn.header{1},3,tw,' '));
            catch err
                a=1;
            end
            strst=sn.header{2};
            if(strfind(strst{2}, '{mm}"'))
                strs{1}=[strst{1} strst{2}];
                strs(2:length(strst)-1)=strst(3:length(strst));
            else
                strs=strst;
            end
            tbl(3,:)=(CommonMethods.cell2tablecells(strs,2,tw,' '));
            len=length(sn.data);
            for i=4:ml
                ii=i-3;
                if(ii<=len)
                    tbl(i,:)=(CommonMethods.cell2tablecells(sn.data{ii},2,tw,' '));
                else
                    tbl(i,:)=(CommonMethods.cell2tablecells(' ',2,tw,' '));                    
                end
            end
        end
        function n=getSummaryNodeSize(~,sn)
            n=length(sn.data)+3;
        end
        function msg=getSNInfo(~,sn)
            folderStr=['folder: ' sn.subfolder];
            conDisp=['con: ' sn.conDisp];
            Img=['Img: ' sn.conImg];
            msg={folderStr conDisp Img};
        end
        function argout=handleTableSelection(obj,table,cellIndices,showStatImg,lockImg)
            if(~exist('lockImg','var'))
                lockImg=false;
            end
            if(~exist('showStatImg','var'))
                showStatImg=false;
            end
            rows=size(cellIndices,1);
            inds=arrayfun(@(x)table.summaryNodeLookupTable(cellIndices(x,1),cellIndices(x,2)),1:rows);
            inds=inds(inds>0);
            if(isempty(inds))
                argout=[];
                return;
            end
            
            sns=arrayfun(@(x)obj.statSummaryNodes{x},inds,'UniformOutput',false);
            argout.msg=obj.getSNInfo(sns{end});
            if(strcmp(table.type,'local'))
                argout.table=table;
            else
                tbl=obj.sns2str(sns);% a cellarray
                [argout.table.data, argout.table.summaryNodeLookupTable]=obj.makeLocalTable(tbl,inds,obj.resultsOutput.snTableWidth);
                argout.table.header=CommonMethods.getDefaultTableHeader(size(argout.table.data,2));
                argout.table.type='local';
            end
            
            if(showStatImg||lockImg)
                coorStr=table.data{cellIndices(end,1),cellIndices(end,2)};
                if(lockImg)
                    CommonMethods.updateSPMExt(coorStr);
                else
                    tmp=obj.displayStatImg(sns{end},coorStr);
                    argout.xjView=tmp.xjView;
                    argout.spmExt=tmp.spmExt;
                end
            end
        end
        function [data,lookupTable]=makeLocalTable(~,tbl,inds,tw)
            len=length(tbl);
            rows=size(tbl{1},1);
            cols=len*tw;
            lookupTable=zeros(rows,cols);
            data=cell(rows,cols);
            for i=1:len
                cI=(i-1)*tw+1;
                cF=i*tw;
                data(:,cI:cF)=tbl{i};
                lookupTable(:,cI:cF)=inds(i);
            end
        end
        function argout=displayStatImg(obj,sn,coorStr)
            [folder,~,~]=fileparts(sn.fnameSPM);
            statMap=fullfile(folder,sn.conImg);
            xjViewFigHandle=xjViewHandler.get_xjView();
            argout.handler=xjViewHandler(statMap,xjViewFigHandle,'+',min(sn.clusterSizes));
            argout.xjView=argout.handler.getXjviewFigHandle();
            
            argout.spmExt=CommonMethods.updateSPMExt(coorStr,sn.fnameSPM);
        end
    end
    
    
    
    %%%%%%%% Static functions
    methods (Static=true)
        function ind=getSummaryNodeIndex(so,imgName,subfolder,conImg)
            ind=find(ismember(so.imgNames,imgName).*ismember(so.subfolders,subfolder).*ismember(so.conImgs,conImg));
        end
        function tbl=getSummaryComparisonTable(so,inds,imgNames)
            tbls=cellfun(@(x)spmStatSummaryHandler.getSummaryTable(so.statSummaryNodes{x}),inds,'UniformOutput',false);
            if(exist('imgNames','var'))
                tbls=arrayfun(@(x)CommonMethods.insertRowToTable(tbls{x}, imgNames{x},1,CommonMethods.tab),1:length(inds),'UniformOutput',false);
            end
            tbl=tbls{1};
            sizes=cellfun(@(x)length(x),tbls);
            for i=2:length(tbls)
                tbl=CommonMethods.horzcatTables(tbl,tbls{i});
                if(length(tbl)>4)
                    strs=cellfun(@(x)strsplit(x,CommonMethods.tab),tbl, 'UniformOutput', false);
                    strsi=cellfun(@(x)strsplit(x,CommonMethods.tab),tbls{i}, 'UniformOutput', false);
                    t=1;
                end
            end
        end
        function tbl=getSummaryTable(sn)
            if(~exist('sn','var'))
                sn=[];
            end
            if(isempty(sn))
                tbl={};
            else
                tbl=vertcat(sn.conLine,sn.headerLines);
                tbl=vertcat(tbl,sn.dataLines);
                if(length(tbl)>3)
                    strs=cellfun(@(x)strsplit(x,CommonMethods.tab),tbl, 'UniformOutput', false);
                    t=1;
                end
            end
        end
        
       function [imgName,pos]=getImgName(sn,pos)
            fname=sn.fnameSPM;
            strs=strsplit(fname,filesep);
            if(~exist('pos','var'))
                [~,poss]=ismember('group',strs);
                pos=poss(1)+1;
            end
            imgName=spmStatSummaryHandler.getImgName_fnameSPM(fname,pos);
        end
        function [imgName, pos]=getImgName_fnameSPM(fname,pos)
            strs=strsplit(fname,filesep);
            if(~exist('pos','var'))
                [~,poss]=ismember('group',strs);
                pos=poss(1)+1;
            end
            imgName=strs{pos};
        end
        function conImg=getConImg(sn)
            strs=strsplit(sn.conLine,CommonMethods.tab);
            conImg=strs{3};
        end
        function conDisp=getConDisp(sn)
            strs=strsplit(sn.conLine,CommonMethods.tab);
            conDisp=strs{2};
        end
        function subfolder=getSubfolder(sn,pos)
            fname=sn.fnameSPM;
            idx=strfind(fname,filesep);
            iI=idx(pos);
            iF=idx(end);
            subfolder=fname(iI+1:iF-1);
        end
    end
end




