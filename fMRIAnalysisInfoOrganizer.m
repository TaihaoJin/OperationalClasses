classdef fMRIAnalysisInfoOrganizer <handle
    %UNTITLED9 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = public)
        AnalysisInfoFileName;
        AnalysisInfo;
        ExcludedScans;
        nameSuffix;
    end
    
    methods
        function obj = fMRIAnalysisInfoOrganizer(study,paradigm, dataRoot, anatRoot,ExcludedScans,repopulate, nameSuffix)
            repopulated=false;
            if(exist('nameModifier','var'))
                obj.nameSuffix=nameSuffix;
            else
                obj.nameSuffix='';
            end
            
            if(exist('study','var'))
                
                if(~isempty(obj.nameSuffix));
                    part=strcat('AnalysisInfo_',study,'_',paradigm,nameSuffix,'.mat');
                else
                    part=strcat('AnalysisInfo_',study,'_',paradigm,'.mat');
                end
                %            AnalysisInfoPath=[dataRoot filesep study filesep part];
                AnalysisInfoPath=fullfile(dataRoot,study,part);
                %  I had problem with fullfile( ... ), which should be equivalent
                obj.AnalysisInfoFileName=AnalysisInfoPath;
                
                if(exist('ExcludedScans', 'var'))
                    obj.ExcludedScans=ExcludedScans;
                else
                    obj.ExcludedScans={};
                end
                
                if(~exist('repopulate', 'var'))
                    repopulate=false;
                end
                if(repopulate)
                    obj.populateAnalysisInfoNode(study,paradigm, dataRoot, anatRoot);
                    repopulated=true;
                elseif(~obj.loadAnalysisIfoNode(AnalysisInfoPath))
                    obj.populateAnalysisInfoNode(study,paradigm, dataRoot, anatRoot);
                    repopulated=true;
                end
            else
                obj.AnalysisInfo=[];
            end
            if(repopulated)
                obj.saveAnalysisInfo();
            end
        end
        
        function populateAnalysisInfoNode(obj, study,paradigm, dataRoot, anatRoot)
            obj.AnalysisInfo=AnalysisInfoNode(study,paradigm,fullfile(dataRoot,study),fullfile(anatRoot,study));
            obj.retrieveScanInfo();
            obj.retrieveScanIDs();
            obj.populateSubjectAndScanNodes();
        end
        
        function status = loadAnalysisIfoNode(obj,path)
           if(exist(path,'file'))
               load (path, 'fi');
               obj.AnalysisInfo=fi;
               status=true;
           else
                status=false;
           end
        end
        
        function path = fMRIImgPath(obj,subjectID,scanID)
            %finding the fMRI image based on data organization pattern.
            if(strcmp(obj.AnalysisInfo.StudyName, 'IRB13278'))
                path=COH1_CommonMethods.getAnalysisdir(scanID,obj.AnalysisInfo.ParadigmName);
            else
                path=fullfile(obj.AnalysisInfo.dataDir, subjectID, scanID, obj.AnalysisInfo.ParadigmName, '1', 'analysis', 'f.nii');
            end
        end
        
        function path = anatImgPath(obj,subjectID,scanID)
            %finding the T1 image based on data organization pattern.
            path=fullfile(obj.anatDir,subjectID, scanID, 'T1', '1', 'analysis');
        end
        
        function retrieveScanInfo(obj)
            if(strcmp(obj.AnalysisInfo.StudyName,'IRB13278'))
                scaninfo=COH1_CommonMethods.getScanInfo(obj.AnalysisInfo.ParadigmName);
            elseif(strcmp(obj.AnalysisInfo.StudyName,'GSK')&&strcmp(obj.AnalysisInfo.ParadigmName,'VOD'))
                scaninfo.NoVolumes=210;
                scaninfo.NoSlices=33;
                scaninfo.RepetitionTime=2;
                scaninfo.NoRuns=6;
                scaninfo.sliceOrder=[1:2:scaninfo.NoSlices 2:2:(scaninfo.NoSlices-1)];
            else
                mysql('open','localhost','root','denali','MRI');
                if(strcmp(obj.AnalysisInfo.StudyName,'Prodrome_sf'))
                    sqlCommand=strcat('select * from scanParameters where StudyName=''',obj.AnalysisInfo.StudyName, '''and ParadigmName=''',obj.AnalysisInfo.ParadigmName,'''');
                else
                    sqlCommand=strcat('select * from scanParameters where StudyName=''',obj.AnalysisInfo.StudyName, '''and ParadigmName=''',obj.AnalysisInfo.ParadigmName,'''');
                end
                
                scaninfo = mysql(sqlCommand);
            end
            
            if(~isempty(scaninfo))
                obj.AnalysisInfo.ScanInfo = scaninfo;
                obj.AnalysisInfo.noVols = scaninfo.NoVolumes;
                obj.AnalysisInfo.noSlices = scaninfo.NoSlices;
                obj.AnalysisInfo.repetitionTime = scaninfo.RepetitionTime;
                obj.AnalysisInfo.noRuns=scaninfo.NoRuns;
            end
        end
        
        function retrieveScanIDs(obj)
            if(strcmp(obj.AnalysisInfo.StudyName,'IRB13278'))
                [scanIDs, table]=COH1_CommonMethods.retrieveScanIDs(obj.AnalysisInfo.ParadigmName);
            else                
                mysql('open','localhost','root','denali','MRI');
                if(strcmp(obj.AnalysisInfo.StudyName,'Prodrome_sf'))
                    %                sqlCommand=strcat('select distinct scanID from analyzeData where studyName=''',obj.AnalysisInfo.StudyName, '''and ParadigmName=''',obj.AnalysisInfo.ParadigmName,''' and include=1 and scanID like ''S1%a%'' order by scanID');
                    %                sqlCommand=strcat('select distinct scanID from analyzeData where studyName=''',obj.AnalysisInfo.StudyName, '''and ParadigmName=''',obj.AnalysisInfo.ParadigmName,''' and include=1 and scanID like ''S1%'' order by scanID');
                    sqlCommand=strcat('select distinct scanID from analyzeData where studyName=''',obj.AnalysisInfo.StudyName, '''and ParadigmName=''',obj.AnalysisInfo.ParadigmName,''' and scanID like ''S1%'' and include=1 order by scanID');
                else
                    sqlCommand=strcat('select distinct scanID from analyzeData where studyName=''',obj.AnalysisInfo.StudyName, '''and ParadigmName=''',obj.AnalysisInfo.ParadigmName,''' and include=1 order by scanID');
                end
                
                %info = mysql('select distinct scanID from analyzeData where studyName=''PTSD'' and paradigmName=''Rest_FA77'' and include=1 order by scanID');
                info = mysql(sqlCommand);
                scanIDs0=sort(unique({info.scanID}));
                scanIDs={};
                
                if(~isempty(obj.ExcludedScans))
                    for i=1:length(scanIDs0)
                        scanID=scanIDs0{i};
                        if(~isempty(find(ismember(obj.ExcludedScans,scanID), 1)))
                            continue;
                        end
                        scanIDs{end+1}=scanID;
                    end
                else
                    scanIDs=scanIDs0;
                end
                numScans=length(scanIDs);
                table=cell(numScans,4);
                
                dates=cell(numScans,1);
                sjIDs=cell(numScans,1);
                for s=1:numScans
                    scanID=scanIDs{s};
                    command=['select subjectID, scandate from scans where scanID = ''' sprintf(scanID) ''''];
                    info=mysql(command);
                    table{s,1}=info.subjectID;
                    table{s,2}=scanID;
                    table{s,3}=info.scandate;
                    table{s,4}=-1;
                    dates{s}=info.scandate;
                    sjIDs{s}=info.subjectID;
                end
                
                table0=table;
                
                [~, idx]=sort(sjIDs);
                for i=1:numScans
                    id=idx(i);
                    for j=1:4
                        table{i,j}=table0{id,j};
                    end
                    scanIDs{i}=table{i,2};
                    dates{i}=table{i,3};
                    sjIDs{i}=table{i,1};
                end
                %/the rows in table is sorted according to subjectID
                table0=table;
                
                for s=1:numScans
                    if(table{s,4}<0)
                        sjID=table{s,1};
                        inds=find(cellfun(@(x)strcmp(x,sjID), sjIDs));
                        
                        strs=dates(inds);
                        [~,idx]=sort(strs);
                        for i=1:length(inds)
                            id=idx(i);
                            table0{inds(id),4}=i;
                        end
                        
                        for i=1:length(inds)
                            id=idx(i);
                            for j=1:4
                                table{inds(i),j}=table0{inds(id),j};
                            end
                        end
                    end
                end
            end
            obj.AnalysisInfo.ScanIDs = scanIDs;
            obj.AnalysisInfo.subjectIDscanIDTable=table;
            %the four columns are subjectID, scanID, scanDate and scan
            %order. scan order is the order of the scan for the same
            %subject. scanOrder=1,2 .. are the first, second ... scan of
            %the same subject.
        end
        function subjectNodes=getSubjectNodes(obj)
            subjectNodes=obj.AnalysisInfo.SubjectNodes;
        end
        function populateSubjectAndScanNodes(obj)
           ScanNodes={};
           SubjectNodes={};
           subjectIDs={};
           scanIDs=obj.AnalysisInfo.ScanIDs;
           for i=1:length(scanIDs)
               if(i==length(scanIDs))
                   i=i;
               end
                sid=scanIDs{i};         
                if(strcmp(obj.AnalysisInfo.StudyName,'PTSD'))
                    sjid=sid;
                else
                    sjid=obj.getSubjectID(sid);
                end
                if(~strcmp(obj.AnalysisInfo.StudyName, 'GSK'))
                    if(~obj.HasScanImage(sjid,sid))
                        continue;
                    end
                end
                                
                ScanNode=obj.buildScanNode(sid);
                len=length(ScanNodes);
                ScanNodes{len+1}=ScanNode;
                fi=find(ismember(subjectIDs,sjid), 1);
                if(isempty(fi))
                    len=length(SubjectNodes);     
                    sjNode = SubjectNode(obj,sjid);
                    SubjectNodes{len+1}=sjNode;
                    subjectIDs{len+1}=sjid;
                else
                    sjNode = SubjectNodes{fi(1)};
                end
                sjNode.addScanNode(ScanNode);
           end
            
           for s=1:length(SubjectNodes)
               sjNode=SubjectNodes{1};
               scNodes=sjNode.ScanNodes;
               num=length(scNodes);
               
               for sc=1:num
                   in=1;
                   st=scNodes{sc};
                   dn=st.ScanDate;
               end
           end
           
            obj.AnalysisInfo.SubjectNodes=SubjectNodes;
            obj.AnalysisInfo.ScanNodes=ScanNodes;
            obj.AnalysisInfo.SubjectIDs=subjectIDs;
        end
        
        function has = HasScanImage(obj,subjectID,scanID)
            %Test whether has the scan for the paradigm
            has = exist(obj.fMRIImgPath(subjectID,scanID),'file');
            if(~has)
                cprintf('red','%s\n', [obj.fMRIImgPath(subjectID,scanID) ' does not exist']);
            end
        end  
        
        function scanNode = buildScanNode(obj,sid)
            [sjid, date,order]=obj.getSubjectID(sid);
            noRuns=obj.AnalysisInfo.noRuns;
            runDirs={};
            for r=1:noRuns
                runDir=fullfile(obj.AnalysisInfo.dataDir,sjid,sid,obj.AnalysisInfo.ParadigmName,num2str(r),'analysis');
                if(exist(fullfile(runDir,'f.nii'),'file'))
                    runDirs{length(runDirs)+1}=runDir;
                end
            end
            scanNode = ScanNode(sid,runDirs);
            scanNode.ScanDate=date;
            scanNode.ScanOrder=order;
            scanNode.SubjectID=sjid;
        end
        
        function saveAnalysisInfo(obj)     
            fi=obj.AnalysisInfo;
            save (obj.AnalysisInfoFileName, 'fi');
        end
        
        function storePreprocessingNode (obj, preproNode)
            preproNodes=obj.StudyInfo.PreprocessingNodes;
            num=length(preproNodes);
            for i=1:num
                if strcmp(preponode.Name,preproNodes{i}.Name) 
                    return;
                end                   
            end
            preproNodes{num+1}=preproNode;
            obj.StudyInfo.PreprocessingNodes=preproNodes;
        end  
        
        function registerPreprocessedImages(obj,proName,proNotes,ScanIDs)
            %ScanIDs is a CellArrayList
            proNode=obj.getPreprocessingNode(proName);
            if(~isa(proNode,PreprocessNode)) 
                proNode=PreprocessNode(proName,proNotes,1,1,1,1);
                obj.AnalysisInfo.PreprocessingNodes.add(proNode);
            end
            for i=0:ScanIDs.length()
                proNode.add(ScanIDs.get(i));
            end
        end
        
        function registerPreprocessNode(obj,preproNode0)
            %ScanIDs is a CellArrayList
            preproNode=obj.getPreprocessingNode(preproNode0.name);
            if(~isa(proNode,PreprocessNode)) 
                obj.AnalysisInfo.PreprossNodes.add(preproNode0);
                return;
            end
            ScanIDs=preproNode0.ScanIDs;
            for i=0:ScanIDs.length()
                preproNode.add(ScanIDs.get(i));
            end
        end
        
        function conn_imgs=getConnImages(obj,roiName)
            %to return connectivity images Rest_FA77 of all scans, grouped
            %by subjectID
            scanIDs=obj.AnalysisInfo.ScanIDs;
            dataDir=obj.AnalysisInfo.dataDir;
            Paradigm=obj.AnalysisInfo.ParadigmName;
            subjectNodes=obj.AnalysisInfo.SubjectNodes;
            conn_imgs={};
            for i=1:length(subjectNodes);
                subjectNode=subjectNodes{i};
                scanNodes=subjectNode.ScanNodes;
                for j=1:length(scanNodes);
                    scanID=scanNodes{j}.ScanID;
                    conn_imgs{end+1}=fullfile(dataDir,subjectNode.SubjectID,scanID,Paradigm,'1','analysis',['beta_' roiName, '.nii']);
                end
            end               
        end
        
        function [IDs, scanIDs, groups]=getSubjectInfo(obj)
            scanIDs=obj.AnalysisInfo.ScanIDs;
            dataDir=obj.AnalysisInfo.dataDir;
            ParaDigm=obj.AnalysisInfo.ParadigmName;
            len=length(scanIDs);
            IDs={' ' len};
            groups={' ' len};
            data_imgs={' ' len};
            for i=1:length(scanIDs);
                scanID=scanIDs{i};
                ID=obj.getSubjectID(scanID);
                
                if(strcmp(scanID(2),'0')||strcmp(scanID(2),'1'))
                    Group='0';
                else
                    Group='1';
                end
                IDs{i}=ID;
                scanIDs{i}=scanID;
                groups{i}=Group;
            end               
        end
        
        function subjectNode = getSubjectNode(obj, scanID)
            subjectNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(subjectNodes)
                if(~isempty(find(ismember(cellfun(@(x) x.ScanID, subjectNodes{i}.ScanNodes,'UniformOutput',false),scanID))))
                    subjectNode=subjectNodes{i};
                    return;
                end
            end          
        end
        
        function subjectNode = getSubjectNode_subjectID(obj, subjectID)
            subjectNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(subjectNodes)
                sNode=subjectNodes{i};
                if(strcmp(sNode.SubjectID, subjectID))
                    subjectNode=sNode;
                    return;
                end
            end          
        end
        
        function nonCompleter = getFirstNonCompleter(obj,cur_scanIDs)
            for i=1:length(cur_scanIDs)
                if(~fMRIAnalysisInfoOrganizer.completedScans(obj.getSubjectNode(cur_scanIDs(i)),{'a' 'b'}))
                    nonCompleter=cur_scanIDs{i};
                    return;
                end
            end
        end
        
        function Completer = getFirstCompleter(obj, cur_scanIDs)
            for i=1:length(cur_scanIDs)
                if(fMRIAnalysisInfoOrganizer.completedScans(obj.getSubjectNode(cur_scanIDs(i)),{'a' 'b'}))
                    Completer=cur_scanIDs{i};
                    return;
                end
            end
        end
         
        function imgPath = getHRFConImagePath(obj, scanID, ParadigmName, conNo, taskSubfolder)
            if(~exist('taskSubfolder','var'))
                taskSubfolder='stats_arf_aCompCor';
            end
            imgPath=fullfile(obj.AnalysisInfo.dataDir, obj.getSubjectID(scanID), scanID, ParadigmName,taskSubfolder, sprintf('swcon_%04d.img', conNo)); 
        end
         
        function no = getHRFConNumber(obj, scanID, ParadigmName, cname, taskSubfolder)
            if(~exist('taskSubfolder','var'))
                taskSubfolder='stats_arf_aCompCor';
            end
            statsdir=fullfile(obj.AnalysisInfo.dataDir, obj.getSubjectID(scanID), scanID, ParadigmName,taskSubfolder); 
            load(fullfile(statsdir, 'SPM.mat'));
            cons=SPM.xCon;
            cnames=arrayfun(@(x) x.name, cons, 'UniformOutput', false);
            idx=find(ismember(cnames, cname));
            if(length(idx)~=1)
                error('wrong contrast name for the function: "getHRFConNumber"');
            end
            no=idx(1);
        end
        
        function imgPath = getHRFConImagePath_twoRuns(obj, scanID, ParadigmName, conNo)
            imgPath=fullfile(obj.AnalysisInfo.dataDir, fMRIAnalysisInfoOrganizer.getSubjectID(scanID), scanID, ParadigmName,'stats_arf_aCompCor2Runs', sprintf('swcon_%04d.img', conNo)); 
        end
        
        function imgPath = getHRFConImagePath_old(obj, scanID, ParadigmName, conNo)
            imgPath=fullfile(obj.AnalysisInfo.dataDir, fMRIAnalysisInfoOrganizer.getSubjectID(scanID), scanID, ParadigmName,'stats_sarf_aCompCor_ARTMotion7_ARTReg', sprintf('wcon_%04d.img', conNo)); 
        end
        
        function dir = getAnalysisDir(obj, scanID,runNo)
            dataDir=obj.AnalysisInfo.dataDir;
            dir=fullfile(dataDir,obj.getSubjectID(scanID), scanID, obj.AnalysisInfo.ParadigmName,num2str(runNo),'analysis');
        end
        
        function selectSubjectNodes(obj,subjectIDs)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            selected={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                if(~isempty(find(ismember(subjectIDs,sNode.SubjectID))))
                    selected{end+1}=sNode;
                end
            end
            obj.AnalysisInfo.updateSubjectNodes(selected);
        end
        
        function selectSubjectNodes_HereAfter(obj,subjectID)
            %repopulated SubjectNodes of obj. SubjectNodes of AnalysisInfo
            %of obj will contain the subject nodes starting from subjectID.
            fi=find(ismember(obj.AnalysisInfo.SubjectIDs,subjectID));
            if(length(fi==1))
                snodes=obj.AnalysisInfo.SubjectNodes(fi(1):length(obj.AnalysisInfo.SubjectNodes));
            end
            obj.AnalysisInfo.updateSubjectNodes(snodes);
        end
        
        function excludeSubjectNodes(obj,subjectIDs)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            selected={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                if(isempty(find(ismember(subjectIDs,sNode.SubjectID))))
                    selected{end+1}=sNode;
                end
            end
            obj.AnalysisInfo.updateSubjectNodes(selected);
        end
        
        function excludeTime2OnlySubjectNodes(obj)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            selected={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                scanNodes=sNode.ScanNodes;
                scanID=scanNodes{1}.ScanID;
                label=obj.getScanLabel(scanID);
                if(strcmp(label,'b'))
                    continue;
                end
                selected{end+1}=sNode;
            end
            obj.AnalysisInfo.updateSubjectNodes(selected);
        end
        
        function selectSubjectNodes_scanLabels(obj,labels)
            %select the subject nodes whose contain all scan labels
            %specified by labels.
            selected={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                scanNodes=sNode.ScanNodes;
                scanLabels=cellfun(@(x)obj.getScanLabel(x.ScanID), scanNodes,'UniformOutput', false);
                if(~all(cellfun(@(x)any(ismember(scanLabels,x)), labels)))
                    continue;
                end
                selected{end+1}=sNode;
            end
            obj.AnalysisInfo.updateSubjectNodes(selected);
        end
        
        function excludeTime1OnlySubjectNodes(obj)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            selected={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                scanNodes=sNode.ScanNodes;
                if(length(scanNodes)<2)
                    continue;
                end
                scanID=scanNodes{end}.ScanID;
                label=obj.getScanLabel(scanID);
                if(strcmp(label,'a'))
                    continue;
                end
                selected{end+1}=sNode;
            end
            obj.AnalysisInfo.updateSubjectNodes(selected);
        end
        
        function excluded=excludeWrongVolumesRuns(obj)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            excluded={};
            adjusted={};
            noVols=obj.AnalysisInfo.noVols;
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                scanNodes=sNode.ScanNodes;
                valid=ones(1,length(scanNodes));
                for j=1:length(scanNodes)
                    selected={};
                    rundir=scanNodes{j}.RunDirs;
                    for r=1:length(rundir)
                        V=spm_vol(fullfile(rundir{r},'f.nii'));
                        if(length(V)~=noVols)
                            excluded{end+1}=rundir{r};
                        else
                            selected{end+1}=rundir{r};
                        end
                    end
                    sNode.ScanNodes{j}.RunDirs=selected;
                    if(isempty(selected))
                        valid(j)=0;
                    end
                end
                idx=find(valid);
                if(isempty(idx))
                    excluded{end+1}=sNode.SubjectID;
                    cprintf('blue', '%s is excluded because of wrong volume number.\n', sNode.SubjectID);
                    continue;
                end
                sNode.ScanNodes=sNode.ScanNodes(idx);
                adjusted{end+1}=sNode;
                cprintf('blue', '%s\n', sNode.SubjectID);
            end
            obj.AnalysisInfo.updateSubjectNodes(adjusted);
        end
        function excluded=excludeScans(obj,scanIDs)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            excluded={};
            adjusted={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                scanNodes=sNode.ScanNodes;
                valid=ones(1,length(scanNodes));
                for j=1:length(scanNodes)
                    scnode=scanNodes{j};
                    valid(j)=~any(ismember(scanIDs, scnode.ScanID));
                end
                idx=find(valid);
                if(isempty(idx))
                    excluded{end+1}=sNode.SubjectID;
                    cprintf('blue', '%s is excluded.\n', sNode.SubjectID);
                    continue;
                end
                sNode.ScanNodes=sNode.ScanNodes(idx);
                adjusted{end+1}=sNode;
%                 cprintf('blue', '%s\n', sNode.SubjectID);
            end
            obj.AnalysisInfo.updateSubjectNodes(adjusted);
        end
        
        function updateSubjectNodes(obj,subjectNodes)
            obj.AnalysisInfo.updateSubjectNodes(subjectNodes);
        end
        
        function selected = getScanIDs_scanLabel(obj,label)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            selected={};
            scanIDs=obj.AnalysisInfo.ScanIDs;
            for i=1:length(scanIDs)
                scanID=scanIDs{i};
                if(strcmp(obj.getScanLabel(scanID),label))
                    selected{end+1}=scanID;
                end
            end
            obj.AnalysisInfo.SubjectNodes=selected;
        end
        
        function sids = getSubjectIDsHavingNScans(obj,n)
            snodes=obj.AnalysisInfo.SubjectNodes;
            sids={};
            for i=1:length(snodes)
                snode=snodes{i};
                if(length(snode.ScanNodes)==n)
                    sids{end+1}=snode.SubjectID;
                end
            end
        end
        
        function varargout = getFileSizeOutliers(obj,fileName)
            %repopulated SubjectNodes of obj.AnalysisInfo with those whose
            %subject IDs are specified in the input subjectIDs
            sizes=[];
            files={};
            sjNodes=obj.AnalysisInfo.SubjectNodes;
            for i=1:length(sjNodes)
                sNode=sjNodes{i};
                runDirs=obj.getRunDirs(sNode);
                for r=1:length(runDirs)
                    file=fullfile(runDirs{r},fileName);
                    p=dir(file);
                    size=p.bytes;
                    files{end+1}=file;
                    sizes(end+1)=size;
                end
            end
            m=mode(sizes);       
            indexes= sizes~=m;
            varargout={files(indexes)};
        end
        
        function [subjectID, scanDate, scanOrder] = getSubjectID(obj,scanID)
            table=obj.AnalysisInfo.subjectIDscanIDTable;
            thiscan=arrayfun(@(x)strcmp(scanID,table{x,2}),1:size(table,1));
            inds=find(thiscan);
            if(length(inds)~=1)
                error('invalid scanID for the function getSubjectID');
            end
            ind=inds(1);
            subjectID=table{ind,1};
            scanDate=table{ind,3};
            scanOrder=table{ind,4};
        end  
        
        function subjectID = getSubjectIDs(obj)
            subjectID=cellfun(@(x) x.SubjectID, obj.AnalysisInfo.SubjectNodes, 'UniformOutput', false);
        end  
        
        function [betas, imgs] = getHRFBetaNames(obj,scanID,ParadigmName, taskfolder)
            %betas: a cell array of unique beta names
            %imgs: a cell array of imgs for each beta, i.e.,
            %Nos{1} is a cell array of images corresponding to the first
            %beta
            if(~exist('taskfolder', 'var'))
                taskfolder='stats_arf_aCompCor';
            end
            path=fullfile(obj.AnalysisInfo.dataDir, obj.getSubjectID(scanID), scanID, ParadigmName, taskfolder);
            load(fullfile(path, 'SPM.mat')); 
            [betas, Nos]=CommonMethods.getBetaNames_SPM(SPM);
            
            imgs={};
            for b=1:length(betas)
                cur_Nos=Nos{b};
                imgs{end+1}=arrayfun(@(x) fullfile(obj.AnalysisInfo.dataDir, obj.getSubjectID(scanID), scanID, ParadigmName,taskfolder, sprintf('swcon_%04d.img', x)), cur_Nos,'UniformOutput', false); 
            end
        end        
    end
    
    methods (Static = true)        
        
        function completed = completedScans(SubjectNode,scanLabels0)
            scanIDs=fMRIAnalysisInfoOrganizer.getScanIDs(SubjectNode);
            scanLabels=cellfun(@(x)fMRIAnalysisInfoOrganizer.getScanLabel(x),scanIDs, 'UniformOutput',false);
            for i=1:length(scanLabels0)
                if (isempty(find(ismember(scanLabels,scanLabels0{i}))))
                    completed = false;
                    return
                end
            end
            completed = true;
        end
        
        function scanIDs = getScanIDs(SubjectNode)
            scanIDs=cellfun(@(x)x.ScanID,SubjectNode.ScanNodes, 'UniformOutput',false);
        end
        
        function label = getScanLabel(scanID)
            indexes=find(isstrprop(scanID,'alpha'));
            label=scanID(indexes(2));
        end
        
        function proNode = getPreprocessingNode(obj, proName)
            proNodes=obj.AnalysisInfo.PreprocessingNodes;
            for i=0:proNodes.length();
                if(strcmp(proName,proNodes.get([i]).name))
                    proNode=proNodes.get([i]);
                    return;
                end
            end
            proNode='none';
        end
        
        function scanpairs = getScanIDsBeforAndAfterTreatment(scanIDs)
            before=0;%scanID before the treatment has been chosen
            scanpairs={};
            for s=1:length(scanIDs)
                if(before==0)
                    if(strcmp(fMRIAnalysisInfoOrganizer.getScanLabel(scanIDs{s}),'a'))
                       scanpairs{1}=scanIDs{s};
                        before=true;
                    end
                else
                    if(strcmp(fMRIAnalysisInfoOrganizer.getScanLabel(scanIDs{s}),'b'))
                        scanpairs{2}=scanIDs{s};
                        break;
                    end
                end
            end
        end
        
        function [rundirs, scanIDs, runNos]= getRunDirs(subjectNode,maxRunsPerScan)
            dirs={};
            if(isa(subjectNode,'ScanNode'))
                scanNodes={subjectNode};
            else
                scanNodes=subjectNode.ScanNodes;
            end
            scanIDs={};
            runNos={};
            noRuns=1;
            for i=1:length(scanNodes)
                for j=1:length(scanNodes{i}.RunDirs)
                    dirs{noRuns}=scanNodes{i}.RunDirs{j};
                    scanIDs{noRuns}=scanNodes{i}.ScanID;
%                     if(exist('maxRunsPerScan','var'))
%                         if(noRuns>=maxRunsPerScan)
%                             break;
%                         end
%                     end
                    file=fullfile(dirs{noRuns},'f.nii');
                    runNos{noRuns}=j;
                    if(~exist(file,'file'))
                        sizes(noRuns)=0;
                    else
                        s = dir(file);
                        sizes(noRuns)=s.bytes;
                    end
                    noRuns=noRuns+1;
                end
            end
            
            md=median(sizes);
            ds=abs(sizes-md);
            tol=md/10000000;
            rundirs=dirs(ds<tol);
            scanIDs=scanIDs(ds<tol);
            runNos=runNos(ds<tol);
        end
        
        %Extracting run number from a run dir
        function runNo = getRunNumber(rundir,ParadigmName)
            parts=strsplit(filesep, rundir);
            len=length(parts);
            if(strcmp(parts(len),'analysis')&&strcmp(parts(len-2),ParadigmName))
                runNo=parts(len-1);
            else
                runNo=NaN;
            end
            if(iscell(runNo))
                runNo=runNo{1};
            end
        end
        function scanID = getScanID(rundir,ParadigmName)
            parts=strsplit(filesep, rundir);
            len=length(parts);
            if(strcmp(parts(len-2),ParadigmName))
                scanID=parts(len-3);
            else
                scanID='';
            end
            if(iscell(scanID))
                scanID=scanID{1};
            end
        end
        function dummies = getRunDummyVariables(subjectNode,maxRunsPerScan)
            %returns dummy variables for each run according to the scan
            %number.
            scanNodes=subjectNode.ScanNodes;
            dummies=[];
            for i=1:length(scanNodes)
                noRuns=1;
                runDirs=fMRIAnalysisInfoOrganizer.getRunDirs(scanNodes{i});
                for j=1:length(runDirs)
                    dummies(end+1)=i;
                    if(nargin==2)
                        if(noRuns>=maxRunsPerScan)
                            break;
                        end
                    end
                    noRuns=noRuns+1;
                end
            end
        end
                
        function cons = getHRFConNames(obj,scanID,ParadigmName,taskSubfolder)
            load(fullfile(obj.AnalysisInfo.dataDir, obj.getSubjectID(scanID), scanID, ParadigmName,taskSubfolder, 'SPM.mat')); 
            cons=cellfun(@(x) x,{SPM.xCon.name}, 'UniformOutput',false);
        end
        
        function unprep = UnpreprocessedSubject(subjectNode)
            maxRuns=inf;
            rundirs=fMRIAnalysisInfoOrganizer.getRunDirs(subjectNode, maxRuns);
            unprep=false;
            for r=1:length(rundirs)
                if(~exist(rundirs{r},'file'))
                    unprep=true;
                    return;
                end
            end
        end
        
        function scanID = getScanID_dataFile(datafile,paradigmName)
            parts=strsplit(filesep,datafile);
            fi=find(ismember(parts,paradigmName));
            if(length(fi)~=1)
                cprintf('red','%s\n', 'wrong inputfile to extract scanID');
                scanID='';
                return;
            end
            scanID=parts{fi(1)-1};           
        end
        
        function runNo = getRunNo_dataFile(datafile,paradigmName)
            parts=strsplit(filesep,datafile);
            fi=find(ismember(parts,paradigmName));
            if(length(fi)~=1)
                cprintf('red','%s\n', 'wrong inputfile to extract scanID');
                runNo='';
                return;
            end
            runNo=parts{fi(1)+1};           
        end
        
        function subjectID = getSubjectID_dataFile(datafile,paradigmName)
            parts=strsplit(filesep,datafile);
            fi=find(ismember(parts,paradigmName));
            if(length(fi)~=1)
                cprintf('red','%s\n', 'wrong inputfile to extract scanID');
                subjectID='';
                return;
            end
            subjectID=parts{fi(1)-2};           
        end
    end
end

