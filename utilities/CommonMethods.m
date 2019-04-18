classdef CommonMethods
    %/home/JinT/matlab_projects/fMRI/OperationalClasses/utilities/CommonMethods
    %Author: Taihao Jin
    %Purpose: This class is to hold static utility methods.
    properties
    end
    
    properties (Constant = true)
        tab=sprintf('%s\t','');
    end
    
    methods (Static = true)
        
        function y = LinearInterpolation(x0, y0, x1, y1, x)
            %this method returns the value of y at x on a line define by
            %the two points (x0, y0) and (x1, y1).
            y=y0+(x-x0)*(y1-y0)/(x1-x0);
        end
        function cstrs = strsplit_CellArray(str, spliters)
            %this function split the string "str" using a cell array of
            %spliters.
            if(~iscell(spliters))
                spliters={spliters};
            end
            if(~iscell(str))
                cstrs={str};
            else
                cstrs=str;
            end
            
            for i=1:length(spliters)
                cstrs0=cstrs;
                cstrs={};
                spliter=char(spliters{i});
                for j=1:length(cstrs0)
                    cstrs=horzcat(cstrs, strsplit(cstrs0{j}, spliter));
                end
            end
            cstrs0=cstrs;
            cstrs={};
            for i=1:length(cstrs0)
                cstrs{i}=char(cstrs0{i});
            end
        end
        
        function pos = findstr_StrInCellArray(strs, subStr)
            %strs is a cell array of string, and subStr is a string
            pos=cellfun(@(x) findstr(x,subStr), strs, 'UniformOutput',false);
        end
        
        function pos = findstr_CellArrayInStr(str, subStrs)
            %str is a string and subStr is a cellarray of strings
            pos=cellfun(@(x) findstr(str,x ), subStrs, 'UniformOutput',false);
        end
        
        function ctable = getStringArrayFromCSVFile(infile)
            %this function returns the contents of a CSV file "infile" as a
            %cell array. Each element of the cell array is a cell array
            %containing the content of a single line in "infile".
            fid = fopen(infile);
            lines = textscan(fid, '%s', 'Delimiter','\n');
            lines=lines{1};
            w=length(strsplit(lines{1},','));
            ctable=cell(length(lines),w);
            
            for i=1:length(lines)
                cline=strsplit(char(lines{i}),',','CollapseDelimiters',false);
                cline=cellfun(@(x) strrep(strtrim(x),'"',''),cline,'UniformOutput',false);
                ctable(i,1:length(cline))=cline;
            end
            fclose(fid);
        end
        function ctable = readDelimitedTextFile(infile,delimiter,charsToStrip)
            %this function returns a 2D cell array. It assumes each line in
            %the file containes equal number of delimiters.
            if(~exist('delimiter','var'))
                delimiter=',';
            end
            if(~exist('charsToStrip','var'))
                charsToStrip='';
            end
            
            if(~isempty(charsToStrip))
                pat=charsToStrip(1);
                for i=2:length(charsToStrip)
                    pat=[pat '|' charsToStrip(i)];
                end
            else
                pat='';
            end
            fid = fopen(infile);
            lines = textscan(fid, '%s', 'Delimiter','\n');
            lines=lines{1};
            w=length(strsplit(lines{1},delimiter));
            
            ctable=cell(length(lines),w);
            
            for i=1:length(lines)
                cline=strsplit(char(lines{i}),delimiter,'CollapseDelimiters',false);
                if(~isempty(pat))
                    cline=cellfun(@(x) regexprep(x,pat,''),cline,'UniformOutput',false);
                end
                ctable(i,1:length(cline))=cline;
            end
            fclose(fid);
        end
        
        function lines=getTextLines(fname)
            try
                fid = fopen(fname);
                lines = textscan(fid, '%s', 'Delimiter','\n');
                lines=lines{1};
            catch err
                lines=[];
                return;
            end
        end
        
        function lines=writeTextLines(fname,lines)
            folder=fileparts(fname);
            if(~exist(folder,'file'))
                mkdir(folder,'file');
            end
            fid=fopen(fname,'wt');
            for l=1:length(lines)
                fprintf(fid,'%s\n',lines{l});
            end
            fclose(fid);
        end
        
        function output = extractPATHDQA(paradigmName,parName)
            %this function generates a csv file containing "parName"
            %of all scans in PATHD.
            studyName='PATHD';
            %            parName = 'mean_sfnr_middle_slice';
            if(~exist('paradigmName','var'))
                paradigmName='Rest_FA77';
            end
            
            if(~exist('parName','var'))
                parName = 'mean_sfnr_middle_slice';
            end
            mysql('open','localhost','root','denali','MRI');
            sqlCommand=strcat('select * from dataQA where StudyName=''',studyName, '''and ParadigmName=''',paradigmName,'''');
            
            selection=mysql(sqlCommand);
            spreadsheet=['/home2/data/process/PATHD/GroupAnalysis/QA/EMOID_' parName '.csv'];
            fid = fopen(spreadsheet,'wt');
            
            header='studyName, scanID, paradigmName, seriesNumber, mean_sfnr_middle_slice';
            line='PATHD, ';
            parNames=['studyName, ' 'scanID, ' 'paradigmName, ' 'seriesNumber, ' 'mean_sfnr_middle_slice \n'];
            
            fprintf(fid,parNames);
            output=cell(length(selection),length(strsplit(parNames,',')));
            for i=1:length(selection)
                s=selection(i);
                line=[s.studyName ', ' s.scanID ', ' s.paradigmName ', ' num2str(s.seriesNumber) ', ' num2str(s.mean_sfnr_middle_slice) '\n'];
                fprintf(fid,line);
                output(i,:)={s.studyName  s.scanID  s.paradigmName  num2str(s.seriesNumber) num2str(s.mean_sfnr_middle_slice)};
            end
            fclose(fid);
        end
        
        function populateQAPars_PATHD_Rest()
            infile='/home2/data/process/PATHD/PATHD_restsubj_SUBJLIST_1214.csv';
            data=CommonMethods.getStringArrayFromCSVFile(infile);
            table=CommonMethods.extractPATHDQA('Rest_FA77');
            subjectIDs1=data(:,1);
            for i=2:length(subjectIDs1)
                n=str2num(subjectIDs1{i});
                if(isnumeric(n))
                    s=sprintf('P%04d',n);
                    subjectIDs1{i}=s;
                end
            end
            scanIDs2=table(:,2);
            pars=table(:,5);
            InfoOrganizer=PATHD_getInfoOrganizer('Rest_FA77');
            subjectNodes=InfoOrganizer.AnalysisInfo.SubjectNodes;
            dim=size(data);
            w0=dim(2);
            newData=cell(dim(1),dim(2)+2);
            newData(:,1:dim(2))=data;
            missing1={};
            missing2={};
            newData{1,4}='sfnr_scana';
            newData{1,5}='sfnr_scanb';
            newData{2,4}='100.1';
            for s=1:length(subjectNodes)
                sjNode=subjectNodes{s};
                sjID=sjNode.SubjectID;
                fi=find(ismember(subjectIDs1,sjID));
                if(length(fi)~=1)
                    missing1{end+1}=sjID;
                    continue;
                end
                index1=fi(1);
                scNodes=sjNode.ScanNodes;
                for sc=1:length(scNodes)
                    sNode=scNodes{sc};
                    scanID=sNode.ScanID;
                    fi=find(ismember(scanIDs2,scanID));
                    if(length(fi)~=1)
                        missing2{end+1}=scanID;
                        if(strcmp(sjID,'P4879'))
                            newData{index1,4}='49.8';
                            newData{index1,5}='59.1';
                        end
                        continue;
                    end
                    index2=fi(1);
                    par=pars{index2};
                    newData{index1,w0+sc}=par;
                end
            end
            rowNames={};
            columnNames={};
            outdir='/home2/data/process/PATHD/QA';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            outfile=fullfile(outdir,'PATHD_restsubj_SUBJLIST_1214_sfnr.csv');
            CommonMethods.exportTable_CSV(outfile,columnNames,rowNames,newData);
        end
        
        function populateQAPars_PATHD_Task()
            dir = '/home2/data/process/PATHD/QA';
            fName= 'pathd_task_subjects';
            paradigmName='EMOID';
            parName='sfnr';
            
            infile=fullfile(dir,[fName '.csv']);
            data=CommonMethods.getStringArrayFromCSVFile(infile);
            table=CommonMethods.extractPATHDQA(paradigmName);
            subjectIDs1=data(:,1);
            for i=1:length(subjectIDs1)
                n=str2num(subjectIDs1{i});
                if(isnumeric(n))
                    s=sprintf('P%04d',n);
                    subjectIDs1{i}=s;
                end
            end
            scanIDs2=table(:,2);
            pars=table(:,5);
            InfoOrganizer=PATHD_getInfoOrganizer(paradigmName);
            subjectNodes=InfoOrganizer.AnalysisInfo.SubjectNodes;
            dim=size(data);
            w0=dim(2);
            
            noRun=3;
            if(strcmp(paradigmName,'EMOWM'))
                noRun=8;
            end
            
            addCols=2*noRun;
            
            newData=cell(dim(1),dim(2)+addCols);
            newData(:,1:dim(2))=data;
            missing1={};
            missing2={};
            
            columnNames=cell(1,w0+addCols);
            columnNames(1:2)={'subject' 'group'};
            labels='ab';
            n=w0;
            for s=1:2
                for r=1:noRun
                    n=n+1;
                    columnNames{1,n}=[parName '_' paradigmName labels(s) num2str(r)];
                end
            end
            
            %            newData{2,4}='100.1';
            for s=1:length(subjectNodes)
                sjNode=subjectNodes{s};
                sjID=sjNode.SubjectID;
                fi=find(ismember(subjectIDs1,sjID));
                if(length(fi)~=1)
                    missing1{end+1}=sjID;
                    continue;
                end
                index1=fi(1);
                scNodes=sjNode.ScanNodes;
                col=w0;
                for sc=1:length(scNodes)
                    sNode=scNodes{sc};
                    scanID=sNode.ScanID;
                    fi=find(ismember(scanIDs2,scanID));
                    if(length(fi)==0)
                        missing2{end+1}=scanID;
                        col=col+noRun;
                        if(strcmp(scanID,'P3028b')&&strcmp(paradigmName,'EMOID'))
                            newData(index1,3:5)={'89.1' '89.0' '75.4'};
                            col=col+3;
                        end
                        continue;
                    end
                    par=pars(fi);
                    newData(index1,col+1:col+length(par))=par;
                    col=col+noRun;
                end
            end
            rowNames={};
            if(~exist(dir,'file'))
                mkdir(outdir);
            end
            
            outfile=fullfile(dir, [fName '_' paradigmName '_' parName '.csv']);
            CommonMethods.exportTable_CSV(outfile,columnNames,rowNames,newData);
        end
        
        function st=coor2str(coor)
            st='';
            for i=1:length(coor)
                c=coor(i);
                if(c<=0)
                    st=[st 'n' num2str(-c)];
                else
                    st=[st 'p' num2str(c)];
                end
            end
        end
        
        function extractGroupROIs()
            %this function perform ROI extractions on the scans in a group.
            %the group, contrast and rois are hard coded. It should be
            %rewritten to more general form if this type of tasks will be
            %performed frequently.
            
            % OUTPUT:
            % p.mean_act:     matrix of mean activations
            % p.median_act:   matrix of median activations
            % p.n_act:        number of voxels in mask (non NaN values)
            % p.min_act:      matrix of minimum (peak values) activations
            % p.max_act:      matrix of maximum (peak values) activations
            % p.v_roi:        SPM V structure for the ROIs
            % p.v_data:       SPM V structure for the data images
            % p.data:         All voxel values in ROI from which stats were extracted.
            %
            % p.int.mean(count,r)
            % p.int.median(count,r)
            % p.int.max(count,r)
            % p.int.min(count,r)
            % p.int.n(count,r)
            
            ParadigmName='EMOID'
            InfoOrganizer=PATHD_getInfoOrganizer(ParadigmName);
            groups=PATHD_Grouping.buildGroups_excluding_NoexistingData(ParadigmName);
            con='nobserveT1';
            %            cur_scanIDs=horzcat(groups{1}.ScanIDs, groups{2}.ScanIDs);
            cur_scanIDs=groups{1}.ScanIDs;
            
            scanID=InfoOrganizer.getFirstNonCompleter(groups{1}.ScanIDs);
            if(iscell(scanID))
                scanID=char(scanID);
            end
            cons = InfoOrganizer.getHRFConNames(InfoOrganizer,scanID,ParadigmName);
            fi=find(ismember(cons, con));
            c=fi(1);
            
            data_imgs = cellfun(@(x) InfoOrganizer.getHRFConImagePath(x, ParadigmName, c),cur_scanIDs,'UniformOutput',false);
            data_imgs_old = cellfun(@(x) InfoOrganizer.getHRFConImagePath_old(x, ParadigmName, c),cur_scanIDs,'UniformOutput',false);
            data_imgs_twoRuns = cellfun(@(x) InfoOrganizer.getHRFConImagePath_twoRuns(x, ParadigmName, c),cur_scanIDs,'UniformOutput',false);
            
            noneExist=[];
            img=data_imgs{1};
            for i=1:length(data_imgs_old);
                imgo=data_imgs_old{i};
                if(~exist(imgo,'file'))
                    noneExist(end+1)=i;
                    data_imgs_old{i}=img;
                end
            end
            
            %           roi_imgs={'/home2/data/process/PATHD/ROIs/Amygdala3d_d0_b_TD.nii'};
            stats = CommonMethods.statExtraction(char(roi_imgs), char(data_imgs));
            
            outdir='/home2/data/process/PATHD/GroupAnalysis/ROIStats/Trouble_Shooting';
            outfile=fullfile(outdir,[con '_' groups{1}.Name '_LateralVentricle_from_spmT_0001.csv']);
            
            stats = CommonMethods.statExtraction(char(roi_imgs), char(data_imgs));
            stats_old = CommonMethods.statExtraction(char(roi_imgs), char(data_imgs_old));
            stats_twoRuns = CommonMethods.statExtraction(char(roi_imgs), char(data_imgs_twoRuns));
            %            stats = roi_stat_extraction(roi_imgs, data_imgs, varargin);
            fid=fopen(outfile,'wt');
            %stats = roi_stat_extraction(roi_imgs, data_imgs);
            header='files, max, min,mean,median, files, max, min,mean,median\n';
            fprintf(fid,header);
            for i=1:length(data_imgs);
                line=[data_imgs{i} ', ' num2str(stats.max_act(i,1)) ', ' num2str(stats.min_act(i,1)) ', ' num2str(stats.mean_act(i,1)) ', ' num2str(stats.median_act(i,1)) ];
                fprintf(fid, line);
                line=[', ' num2str(stats_old.max_act(i,1)) ', ' num2str(stats_old.min_act(i,1)) ', ' num2str(stats_old.mean_act(i,1)) ', ' num2str(stats_old.median_act(i,1))];
                fprintf(fid, line);
                line=[', ' num2str(stats_twoRuns.max_act(i,1)) ', ' num2str(stats_twoRuns.min_act(i,1)) ', ' num2str(stats_twoRuns.mean_act(i,1)) ', ' num2str(stats_twoRuns.median_act(i,1)) '\n'];
                fprintf(fid, line);
            end
            fclose(fid);
        end
        
        
        
        function p = statExtraction(roi_imgs, data_imgs, varargin)
            
            %--------------------------------------------------------------------------
            % roi_stat_extraction
            %
            % Extracts mean, median, min, max, and intervals of ROIs
            %
            % USAGE:
            % p = roi_stat_extraction(roi_imgs, data_imgs, interval)
            %
            % INPUT:
            % roi_imgs:     char array of ROI image filenames. Means and medians are
            %               extracted where ROI~=0
            % data_imgs:    char array of data image filenames for which the means and
            %               medians are extracted
            % OPTIONAL INPUTS:
            % interval:     [min_value max_value] specifies a specific percentage interval to
            %               extract statistics from. E.g. [90 100] will extract
            %               mean/median from the 90-100 percentile.
            % direction:    'negative' or'positive' to extract negative or positive
            %               values. The default is both.
            %
            % OUTPUT:
            % p.mean_act:     matrix of mean activations
            % p.median_act:   matrix of median activations
            % p.n_act:        number of voxels in mask (non NaN values)
            % p.min_act:      matrix of minimum (peak values) activations
            % p.max_act:      matrix of maximum (peak values) activations
            % p.v_roi:        SPM V structure for the ROIs
            % p.v_data:       SPM V structure for the data images
            % p.data:         All voxel values in ROI from which stats were extracted.
            %
            % p.int.mean(count,r)
            % p.int.median(count,r)
            % p.int.max(count,r)
            % p.int.min(count,r)
            % p.int.n(count,r)
            %
            % see also: roi_stats_to_csv
            
            
            % load the default spm settings
            spm('defaults', 'fmri');
            
            for arg=1:size(varargin,2)
                if ischar(varargin{arg})
                    direction = varargin{arg};
                elseif isnumeric(varargin{arg}) && size(varargin{arg},2) == 2
                    interval = varargin{arg};
                end
            end
            
            
            hold = -2; % -2=second order sinc interpolation
            
            p.int = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.v_roi  = spm_vol(roi_imgs);
            p.v_data = spm_vol(data_imgs);
            if(iscell(p.v_roi))
                p.v_roi=cell2mat(p.v_roi);
            end
            if(iscell(p.v_data))
                p.v_data=cell2mat(p.v_data);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through the ROIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            maxroilen = size(char(roi_imgs), 2);
            
            for r = 1:size(roi_imgs,1)
                tic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xyz mm of non-zero points in roi img
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Y,XYZ] = spm_read_vols(p.v_roi(r));
                tmp     = find(Y);
                XYZ     = [XYZ(:, tmp); ones(1, length(tmp))];
                clear Y
                status_str  = sprintf(['\nROI: %s - Extracting %d voxels from %d images'], p.v_roi(r).fname(max(strfind(p.v_roi(r).fname, '/'))+1:end), size(XYZ,2), length(p.v_data));
                cprintf('blue','%s\n',status_str);
                pause(0.01);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % loop through the data images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                count = 1;
                for i = 1:length(p.v_data)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find voxel elements for which to extract stats
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    ixyz = p.v_data(i).mat \ XYZ;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % then sample the data from the volume at those locations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if(i==8)
                        i=i;
                    end
                    
                    data = spm_sample_vol(p.v_data(i),ixyz(1,:),ixyz(2,:),ixyz(3,:),hold);
                    
                    if exist('direction', 'var')
                        if strcmp(direction,'positive')
                            data = data(find(data > 0));
                        elseif strcmp(direction,'negative')
                            data = data(find(data < 0));
                        end
                        if isempty(data), data = NaN; end;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find means/medians
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    p.mean_act(count,r)   = nanmean(data);
                    p.median_act(count,r) = nanmedian(data);
                    p.max_act(count,r)    = nanmax(data);
                    p.min_act(count,r)    = nanmin(data);
                    p.n_act(count,r)      = length(find(~isnan(data)));
                    p.nZeros_act(count,r)      = length(find(abs(data)<0.000000001));
                    
                    p.data{count, r,:}    = data;
                    
                    if exist('interval', 'var')
                        
                        data = data(find(~isnan(data)));
                        
                        sorted_data = sort(data);
                        
                        %fprintf([num2str(start_indx) ' ' num2str(end_indx) '\n']);
                        
                        
                        start_indx = round(length(data)*interval(1)/100);
                        end_indx = round(length(data)*interval(2)/100);
                        
                        interval_data = sorted_data(start_indx:end_indx);
                        
                        p.int.mean(count,r)   = nanmean(interval_data);
                        p.int.median(count,r) = nanmedian(interval_data);
                        p.int.max(count,r)    = nanmax(interval_data);
                        p.int.min(count,r)    = nanmin(interval_data);
                        p.int.n(count,r)      = length(find(~isnan(interval_data)));
                    end
                    count = count+1;
                end
                toc
            end
        end
        
        function data = getVols(roi_img, data_img)
            %this function read the voxel values of data_img at every voxel
            %coordinate of roi_img
            % load the default spm settings
            
            %   Return voxel values from an image volume - a compiled routine
            %   FORMAT X = spm_sample_vol(V,x,y,z,hold)
            %   V        -  spm_vol structure
            %   x        -  matrix of x coordinates {voxels}
            %   y        -  matrix of y coordinates {voxels}
            %   z        -  matrix of z coordinates {voxels}
            %   hold     -  interpolation method for the resampling:
            %                0         : Zero-order hold (nearest neighbour)
            %                1         : First-order hold (trilinear interpolation)
            %                2->127    : Higher order Lagrange (polynomial) interpolation
            %                            using different holds (second-order upwards)
            %                -127 - -1 : Different orders of sinc interpolation
            %
            %   X        -  output image
            %
            %   FORMAT [X,dX,dY,dZ] = spm_sample_vol(V,x,y,z,hold)
            %   Similar to above, except that the derivatives in the three orthogonal
            %   directions are also returned.
            
            hold = 0; % -2=second order sinc interpolation
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(ischar(roi_img))
                roi_img=spm_vol(roi_img);
                roi_img=roi_img(1);
            end
            p.v_roi  = roi_img;
            p.v_data = spm_vol(data_img);
            if(length(p.v_data)>1)
                p.v_data=p.v_data(1);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % xyz mm of non-zero points in roi img
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [~,XYZ] = spm_read_vols(p.v_roi);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find voxel elements for which to extract stats
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XYZ=[XYZ;ones(1,size(XYZ,2))];
            ixyz = p.v_data.mat \ XYZ;
            ixyz=ixyz(1:3,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % then sample the data from the volume at those locations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ind1=ixyz(1,:);
            ind2=ixyz(2,:);
            ind3=ixyz(3,:);
            datat = spm_sample_vol(p.v_data,ind1,ind2,ind3,hold);
            data=nan(p.v_roi.dim);
            dim=p.v_roi.dim;
            valid=((ind1.*ind2.*ind3)>0).*ind1<=dim(1).*ind2<=dim(2).*(ind3<=dim(3))>0;
            data(valid)=datat(valid);
        end
        
        function voxels = getVoxels_1Img(data_img, coors)
            % coors should be a 3 by N (number of voxels) matrix;
            %this function read the voxel values of data_img at every voxel
            %coordinate specified by coors
            % load the default spm settings
            
            %   Return voxel values from an image volume - a compiled routine
            %   FORMAT X = spm_sample_vol(V,x,y,z,hold)
            %   V        -  spm_vol structure
            %   x        -  matrix of x coordinates {voxels}
            %   y        -  matrix of y coordinates {voxels}
            %   z        -  matrix of z coordinates {voxels}
            %   hold     -  interpolation method for the resampling:
            %                0         : Zero-order hold (nearest neighbour)
            %                1         : First-order hold (trilinear interpolation)
            %                2->127    : Higher order Lagrange (polynomial) interpolation
            %                            using different holds (second-order upwards)
            %                -127 - -1 : Different orders of sinc interpolation
            %
            %   X        -  output image
            %
            %   FORMAT [X,dX,dY,dZ] = spm_sample_vol(V,x,y,z,hold)
            %   Similar to above, except that the derivatives in the three orthogonal
            %   directions are also returned.
            
            hold = 0; % -2=second order sinc interpolation
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.v_data = spm_vol(data_img);
            if(length(p.v_data)>1)
                p.v_data=p.v_data(1);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % xyz mm of non-zero points in roi img
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XYZ = coors;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % find voxel elements for which to extract stats
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            XYZ=[XYZ;ones(1,size(XYZ,2))];
            ixyz = p.v_data(1).mat \ XYZ;
            ixyz=ixyz(1:3,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % then sample the data from the volume at those locations
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ind1=ixyz(1,:);
            ind2=ixyz(2,:);
            ind3=ixyz(3,:);
            voxels = spm_sample_vol(p.v_data,ind1,ind2,ind3,hold);
        end
        
        function [p, valid] = statExtraction_intersectMask(roi_imgs, data_imgs, mask_img, varargin)
            
            %--------------------------------------------------------------------------
            % roi_stat_extraction
            %
            % Extracts mean, median, min, max, and intervals of ROIs
            %
            % USAGE:
            % p = roi_stat_extraction(roi_imgs, data_imgs, interval)
            %
            % INPUT:
            % roi_imgs:     char array of ROI image filenames. Means and medians are
            %               extracted where ROI~=0
            % data_imgs:    char array of data image filenames for which the means and
            %               medians are extracted
            % OPTIONAL INPUTS:
            % interval:     [min_value max_value] specifies a specific percentage interval to
            %               extract statistics from. E.g. [90 100] will extract
            %               mean/median from the 90-100 percentile.
            % direction:    'negative' or'positive' to extract negative or positive
            %               values. The default is both.
            %
            % OUTPUT:
            % p.mean_act:     matrix of mean activations
            % p.median_act:   matrix of median activations
            % p.n_act:        number of voxels in mask (non NaN values)
            % p.min_act:      matrix of minimum (peak values) activations
            % p.max_act:      matrix of maximum (peak values) activations
            % p.v_roi:        SPM V structure for the ROIs
            % p.v_data:       SPM V structure for the data images
            % p.data:         All voxel values in ROI from which stats were extracted.
            %
            % p.int.mean(count,r)
            % p.int.median(count,r)
            % p.int.max(count,r)
            % p.int.min(count,r)
            % p.int.n(count,r)
            %
            % see also: roi_stats_to_csv
            
            
            % load the default spm settings
            spm('defaults', 'fmri');
            
            for arg=1:size(varargin,2)
                if ischar(varargin{arg})
                    direction = varargin{arg};
                elseif isnumeric(varargin{arg}) && size(varargin{arg},2) == 2
                    interval = varargin{arg};
                end
            end
            
            if(~exist('mask_img','var'))
                mask_img='';
            end
            
            if(exist(mask_img','file'))
                Vm=spm_vol(mask_img);
            else
                Vm=[];
            end
            
            hold = -2; % -2=second order sinc interpolation
            
            p.int = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.mask=mask_img;
            p.v_roi  = spm_vol(roi_imgs);
            p.v_data = spm_vol(data_imgs);
            if(iscell(p.v_roi))
                p.v_roi=cell2mat(p.v_roi);
            end
            if(iscell(p.v_data))
                p.v_data=cell2mat(p.v_data);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through the ROIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            maxroilen = size(char(roi_imgs), 2);
            
            for r = 1:size(roi_imgs,1)
                tic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xyz mm of non-zero points in roi img
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Y,XYZ] = spm_read_vols(p.v_roi(r));
                tmp     = find(Y);
                nRoio=length(tmp);
                XYZ     = [XYZ(:, tmp); ones(1, length(tmp))];
                clear Y
                
                if(~isempty(Vm))
                    ixyzm=Vm.mat \ XYZ;
                    datam=spm_sample_vol(Vm,ixyzm(1,:),ixyzm(2,:),ixyzm(3,:),0);
                    temp=find(datam);
                    XYZ=XYZ(:,temp);
                end
                
                count = 1;
                if(isempty(XYZ))
                    for i = 1:length(p.v_data)
                        p.mean_act(count,r)   = nan;
                        p.median_act(count,r) = nan;
                        p.max_act(count,r)    = nan;
                        p.min_act(count,r)    = nan;
                        p.n_act(count,r)      = 0;%the number of non-nan numbers
                        p.nZeros_act(count,r)      = 0;%the number of zeros
                        p.nRoinan_act(count,r)=0;%the number of nans
                        p.nRoi_act(count,r)=0; %the roi size after intersect with the mask
                        p.nRoio_act(count,r)=nRoio;%the roi size before intersect with the mask
                        p.data{count, r,:}    = [];
                        count=count+1;
                    end
                    continue;
                end
                status_str  = sprintf(['\nROI: %s - Extracting %d voxels from %d images'], p.v_roi(r).fname(max(strfind(p.v_roi(r).fname, '/'))+1:end), size(XYZ,2), length(p.v_data));
                cprintf('blue','%s\n',status_str);
                pause(0.01);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % loop through the data images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                count = 1;
                for i = 1:length(p.v_data)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find voxel elements for which to extract stats
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    ixyz = p.v_data(i).mat \ XYZ;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % then sample the data from the volume at those locations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if(i==8)
                        i=i;
                    end
                    
                    data = spm_sample_vol(p.v_data(i),ixyz(1,:),ixyz(2,:),ixyz(3,:),hold);
                    
                    if exist('direction', 'var')
                        if strcmp(direction,'positive')
                            data = data(find(data > 0));
                        elseif strcmp(direction,'negative')
                            data = data(find(data < 0));
                        end
                        if isempty(data), data = NaN; end;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find means/medians
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    p.mean_act(count,r)   = nanmean(data);
                    p.median_act(count,r) = nanmedian(data);
                    p.max_act(count,r)    = nanmax(data);
                    p.min_act(count,r)    = nanmin(data);
                    p.n_act(count,r)      = length(find(~isnan(data)));
                    p.nZeros_act(count,r)      = length(find(abs(data)<0.000000001));
                    p.nRoinan_act(count,r)=length(isnan(data));
                    p.nRoi_act(count,r)=length(data); %the roi size after intersect with the mask
                    p.nRoio_act(count,r)=nRoio;%the roi size before intersect with the mask
                    p.data{count, r,:}    = data;
                    
                    if exist('interval', 'var')
                        
                        data = data(find(~isnan(data)));
                        
                        sorted_data = sort(data);
                        
                        %fprintf([num2str(start_indx) ' ' num2str(end_indx) '\n']);
                        
                        
                        start_indx = round(length(data)*interval(1)/100);
                        end_indx = round(length(data)*interval(2)/100);
                        
                        interval_data = sorted_data(start_indx:end_indx);
                        
                        p.int.mean(count,r)   = nanmean(interval_data);
                        p.int.median(count,r) = nanmedian(interval_data);
                        p.int.max(count,r)    = nanmax(interval_data);
                        p.int.min(count,r)    = nanmin(interval_data);
                        p.int.n(count,r)      = length(find(~isnan(interval_data)));
                    end
                    
                    
                    count = count+1;
                end
                toc
            end
        end
        
        %         function intersectROIs = rois_masksIntersect(roi_imgs, mask_imgs)
        %             intersectROIs=[];
        %             %implement it later
        %         end
        
        function data = statExtraction_gettingDataMatrix(p)
            %this function rearranges the data in the roistats p into a
            %cellarray of data matrix. each cell in the cellarray data is a
            %data matrix from a single roi. The number of rows of the data matrix (a cell) equals the number of the data images,
            %and the nuber of columns is the number of voxels in the roi image. The i-th row and j-th column is
            %the voxel value of the j-th voxel (of the roi) in the i-th
            %image file.
            data0=p.data;
            %dim(1) is for different data files, and dim(2) are for
            %different rois.
            dim=size(data0);
            data=cell(1,dim(2));
            for r=1:dim(2)
                datar=zeros(dim(1),length(data0{1,r}));
                for i=1:dim(1)
                    datar(i,:)=cell2mat(data0(i,r));
                end
                data{r}=datar;
            end
        end
        
        function p = statExtraction_singleVoxel(voxels, data_imgs, varargin)
            
            %--------------------------------------------------------------------------
            % roi_stat_extraction from single voxels;
            %
            
            
            % load the default spm settings
            spm('defaults', 'fmri');
            
            for arg=1:size(varargin,2)
                if ischar(varargin{arg})
                    direction = varargin{arg};
                elseif isnumeric(varargin{arg}) && size(varargin{arg},2) == 2
                    interval = varargin{arg};
                end
            end
            
            
            hold = -2; % -2=second order sinc interpolation
            
            p.int = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.v_data = spm_vol(data_imgs);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through the ROIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            voxelLabels=cell(1,size(voxels,1));
            for it=1:size(voxels,2)
                voxelLabels{it}=['(' num2str(voxels(1,it)) ', ' num2str(voxels(2,it)) ', ' num2str(voxels(3,it)) ')'];
            end
            
            for r = 1:size(voxels,2)
                tic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xyz mm of non-zero points in roi img
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                XYZ     = [voxels; ones(1, size(voxels,2))];
                status_str  = sprintf(['\nROI: %s ' voxelLabels{r}]);
                cprintf('blue','%s\n',status_str);
                pause(0.01);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % loop through the data images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                count = 1;
                for i = 1:length(p.v_data)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find voxel elements for which to extract stats
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    ixyz = p.v_data(i).mat \ XYZ;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % then sample the data from the volume at those locations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if(i==8)
                        i=i;
                    end
                    
                    data = spm_sample_vol(p.v_data(i),ixyz(1,r),ixyz(2,r),ixyz(3,r),hold);
                    
                    if exist('direction', 'var')
                        if strcmp(direction,'positive')
                            data = data(find(data > 0));
                        elseif strcmp(direction,'negative')
                            data = data(find(data < 0));
                        end
                        if isempty(data), data = NaN; end;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find means/medians
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    p.mean_act(count,r)   = nanmean(data);
                    p.data{count, r,:}    = data;
                    
                    if exist('interval', 'var')
                        
                        data = data(find(~isnan(data)));
                        
                        sorted_data = sort(data);
                        
                        %fprintf([num2str(start_indx) ' ' num2str(end_indx) '\n']);
                        
                        
                        start_indx = round(length(data)*interval(1)/100);
                        end_indx = round(length(data)*interval(2)/100);
                        
                        interval_data = sorted_data(start_indx:end_indx);
                        
                        p.int.mean(count,r)   = nanmean(interval_data);
                    end
                    
                    
                    count = count+1;
                end
                toc
            end
        end
        
        function p = statExtraction_Timeseries_withFirstPC(roi_imgs, data_img, compare_spm)
            % Modified by Taihao Jin to include the first principal
            % component signal
            %--------------------------------------------------------------------------
            % roi_stat_extraction
            %
            % Extracts mean, median, min, max, and intervals of ROIs
            %
            % USAGE:
            % p = roi_stat_extraction(roi_imgs, data_imgs, interval)
            %
            % INPUT:
            % roi_imgs:     char array of ROI image filenames. Means and medians are
            %               extracted where ROI~=0
            % data_imgs:    char array of data image filenames for which the means and
            %               medians are extracted
            % OPTIONAL INPUTS:
            % interval:     [min_value max_value] specifies a specific percentage interval to
            %               extract statistics from. E.g. [90 100] will extract
            %               mean/median from the 90-100 percentile.
            % direction:    'negative' or'positive' to extract negative or positive
            %               values. The default is both.
            %
            % OUTPUT:
            % p.mean_act:     matrix of mean activations
            % p.median_act:   matrix of median activations
            % p.n_act:        number of voxels in mask (non NaN values)
            % p.min_act:      matrix of minimum (peak values) activations
            % p.max_act:      matrix of maximum (peak values) activations
            % p.v_roi:        SPM V structure for the ROIs
            % p.v_data:       SPM V structure for the data images
            % p.data:         All voxel values in ROI from which stats were extracted.
            %
            % p.int.mean(count,r)
            % p.int.median(count,r)
            % p.int.max(count,r)
            % p.int.min(count,r)
            % p.int.n(count,r)
            %
            % see also: roi_stats_to_csv
            
            
            % load the default spm settings
            spm('defaults', 'fmri');
            
            hold = -2; % -2=second order sinc interpolation
            
            p.int = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.v_roi  = spm_vol(roi_imgs);
            p.v_data = spm_vol(data_img);
            vol=spm_read_vols(p.v_data);
            dim=size(vol);
            volume=dim(1)*dim(2)*dim(3);
            dataSize=length(p.v_data);
            p.RoiSizes=-1*ones(1,size(roi_imgs,1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through the ROIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for r = 1:size(roi_imgs,1)
                tic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xyz mm of non-zero points in roi img
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Y,XYZ] = spm_read_vols(p.v_roi(r));
                tmp     = find(Y);
                XYZ     = [XYZ(:, tmp); ones(1, length(tmp))];
                clear Y
                cprintf('blue','%s\n',sprintf('\nROI: %s - Extracting %d voxels from %d images', p.v_roi(r).fname(max(strfind(p.v_roi(r).fname, '/'))+1:end), size(XYZ,2), length(p.v_data)));
                ixyz = round(p.v_data(1).mat \ XYZ);
                inds=sub2ind(dim,ixyz(1,:),ixyz(2,:),ixyz(3,:));
                inds0=unique(inds);
                p.RoiSizes(r)=length(inds0);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % extracting time series
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                count = 1;
                
                %number of data images
                %number of voxels in the ROI
                roiSize=length(inds0);
                
                Mdata=zeros(dataSize, roiSize);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find voxel elements for which to extract stats
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for i = 1:dataSize
                    offset=(i-1)*volume;
                    indst=inds0+offset;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % then sample the data from the volume at those locations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    data = vol(indst)';
                    test=testMat(i);
                    %                    data = spm_sample_vol(p.v_data(i),ixyz(1,:),ixyz(2,:),ixyz(3,:),hold);
                    
                    %                     if exist('direction', 'var')
                    %                         if strcmp(direction,'positive')
                    %                             data = data(find(data > 0));
                    %                         elseif strcmp(direction,'negative')
                    %                             data = data(find(data < 0));
                    %                         end
                    %                         if isempty(data), data = NaN; end;
                    %                     end
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find means/medians
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    p.mean_act(count,r)   = nanmean(data);
                    p.median_act(count,r) = nanmedian(data);
                    p.max_act(count,r)    = nanmax(data);
                    p.min_act(count,r)    = nanmin(data);
                    p.n_act(count,r)      = length(find(~isnan(data)));
                    p.nZeros_act(count,r)      = length(find(abs(data)<0.000000001));
                    
                    p.data{count, r,:}    = data;
                    
                    if exist('interval', 'var')
                        
                        data = data(find(~isnan(data)));
                        
                        sorted_data = sort(data);
                        
                        %fprintf([num2str(start_indx) ' ' num2str(end_indx) '\n']);
                        
                        
                        start_indx = round(length(data)*interval(1)/100);
                        end_indx = round(length(data)*interval(2)/100);
                        
                        interval_data = sorted_data(start_indx:end_indx);
                        
                        p.int.mean(count,r)   = nanmean(interval_data);
                        p.int.median(count,r) = nanmedian(interval_data);
                        p.int.max(count,r)    = nanmax(interval_data);
                        p.int.min(count,r)    = nanmin(interval_data);
                        p.int.n(count,r)      = length(find(~isnan(interval_data)));
                    end
                    
                    
                    count = count+1;
                    Mdata(i,:)=data;
                    
                    offset=offset+volume;
                end
                
                nanColumns=arrayfun(@(x) any(isnan(Mdata(:,x))), 1:roiSize);
                numericColumns=find(~nanColumns);
                len=length(numericColumns);
                if(len>1)
                    %The extracted data that do not contains Nan at any rows
                    Mdata1=zeros(dataSize,len);
                    
                    for col=1:len
                        Mdata1(:,col)=Mdata(:,numericColumns(col));
                    end
                    [Signal, pc, S] = CommonMethods.getPCs(Mdata1);
                    %signal: mXn matrix of projected data
                    %pc: each column is a PC, the eighen vector of the covarience
                    %matrix
                    %S: nX1 matrix of Variance
                    
                    d=sign(sum(pc(:,1)));
                    PC=Signal(:,1)*d;
                else
                    PC=nan*ones(dataSize);
                end
                
                if(compare_spm)
                    [coeff, score] = princomp(Mdata1);
                    d=sign(sum(coeff(:,1)));
                    PC_matlab=d*score(:,1);
                    [coeff_spm, score_spm] = CommonMethods.spm_princomp(Mdata1);
                    d=sign(sum(coeff_spm(:,1)));
                    PC_spm=d*score_spm(:,1);
                    
                    vLoading=zeros(dim(1:3));
                    vLoading_spm=zeros(dim(1:3));
                    for it=1:length(inds0)
                        vLoading(inds0(it))=coeff(it,1);
                        vLoading_spm(inds0(it))=coeff_spm(it,1);
                    end
                    p.vLoading_act{r}=vLoading;
                    p.vLoading_spm_act{r}=vLoading_spm;
                end;
                
                for row=1:dataSize
                    p.PC_act(row,r)=PC(row);
                    if(compare_spm)
                        p.PC_act_spm(row,r)=PC_spm(row);
                        p.PC_act_matlab(row,r)=PC_matlab(row);
                    end
                end
            end
            
            function mat=testMat(i)
                pos=[1 4 7 35 34 15];
                num=length(pos);
                mat=nan*ones(1,num);
                for pt=1:num
                    mat(pt)=vol(inds(pt)+offset)-vol(ixyz(1,pt),ixyz(2,pt),ixyz(3,pt),i);
                end
            end
        end
        
        function p = statExtraction_withFirstPC(roi_imgs, data_imgs, varargin)
            % Modified by Taihao Jin to include the first principal
            % component score
            % p.PC_act, a matrix holding first PC component scores. The
            % number of columns is the number of roi_imgs and the number of
            % rows are the number of volumns in data_imgs. It mainly used
            % for extracting a single 4D image or a collection of 3D
            % images.
            
            
            %--------------------------------------------------------------------------
            % roi_stat_extraction
            %
            % Extracts mean, median, min, max, and intervals of ROIs
            %
            % USAGE:
            % p = roi_stat_extraction(roi_imgs, data_imgs, interval)
            %
            % INPUT:
            % roi_imgs:     char array of ROI image filenames. Means and medians are
            %               extracted where ROI~=0
            % data_imgs:    char array of data image filenames for which the means and
            %               medians are extracted
            % OPTIONAL INPUTS:
            % interval:     [min_value max_value] specifies a specific percentage interval to
            %               extract statistics from. E.g. [90 100] will extract
            %               mean/median from the 90-100 percentile.
            % direction:    'negative' or'positive' to extract negative or positive
            %               values. The default is both.
            %
            % OUTPUT:
            % p.mean_act:     matrix of mean activations
            % p.median_act:   matrix of median activations
            % p.n_act:        number of voxels in mask (non NaN values)
            % p.min_act:      matrix of minimum (peak values) activations
            % p.max_act:      matrix of maximum (peak values) activations
            % p.v_roi:        SPM V structure for the ROIs
            % p.v_data:       SPM V structure for the data images
            % p.data:         All voxel values in ROI from which stats were extracted.
            %
            % p.int.mean(count,r)
            % p.int.median(count,r)
            % p.int.max(count,r)
            % p.int.min(count,r)
            % p.int.n(count,r)
            %
            % see also: roi_stats_to_csv
            
            
            % load the default spm settings
            spm('defaults', 'fmri');
            
            for arg=1:size(varargin,2)
                if ischar(varargin{arg})
                    direction = varargin{arg};
                elseif isnumeric(varargin{arg}) && size(varargin{arg},2) == 2
                    interval = varargin{arg};
                end
            end
            
            
            hold = -2; % -2=second order sinc interpolation
            
            p.int = [];
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            p.v_roi  = spm_vol(roi_imgs);
            p.v_data = spm_vol(data_imgs);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through the ROIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            maxroilen = size(char(roi_imgs), 2);
            
            for r = 1:size(roi_imgs,1)
                tic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xyz mm of non-zero points in roi img
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Y,XYZ] = spm_read_vols(p.v_roi(r));
                tmp     = find(Y);
                XYZ     = [XYZ(:, tmp); ones(1, length(tmp))];
                clear Y
                status_str  = sprintf(['\nROI: %s - Extracting %d voxels from %d images'], p.v_roi(r).fname(max(strfind(p.v_roi(r).fname, '/'))+1:end), size(XYZ,2), length(p.v_data));
                pause(0.01);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % loop through the data images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                count = 1;
                
                %number of data images
                dataSize=length(p.v_data);
                %number of voxels in the ROI
                roiSize=length(XYZ);
                
                Mdata=zeros(dataSize, roiSize);
                
                for i = 1:length(p.v_data)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find voxel elements for which to extract stats
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    ixyz = p.v_data(i).mat \ XYZ;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % then sample the data from the volume at those locations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    data = spm_sample_vol(p.v_data(i),ixyz(1,:),ixyz(2,:),ixyz(3,:),hold);
                    
                    if exist('direction', 'var')
                        if strcmp(direction,'positive')
                            data = data(find(data > 0));
                        elseif strcmp(direction,'negative')
                            data = data(find(data < 0));
                        end
                        if isempty(data), data = NaN; end;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find means/medians
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    p.mean_act(count,r)   = nanmean(data);
                    p.median_act(count,r) = nanmedian(data);
                    p.max_act(count,r)    = nanmax(data);
                    p.min_act(count,r)    = nanmin(data);
                    p.n_act(count,r)      = length(find(~isnan(data)));
                    p.nZeros_act(count,r)      = length(find(abs(data)<0.000000001));
                    
                    p.data{count, r,:}    = data;
                    
                    if exist('interval', 'var')
                        
                        data = data(find(~isnan(data)));
                        
                        sorted_data = sort(data);
                        
                        %fprintf([num2str(start_indx) ' ' num2str(end_indx) '\n']);
                        
                        
                        start_indx = round(length(data)*interval(1)/100);
                        end_indx = round(length(data)*interval(2)/100);
                        
                        interval_data = sorted_data(start_indx:end_indx);
                        
                        p.int.mean(count,r)   = nanmean(interval_data);
                        p.int.median(count,r) = nanmedian(interval_data);
                        p.int.max(count,r)    = nanmax(interval_data);
                        p.int.min(count,r)    = nanmin(interval_data);
                        p.int.n(count,r)      = length(find(~isnan(interval_data)));
                    end
                    
                    
                    count = count+1;
                    Mdata(i,:)=data;
                end
                nanColumns=arrayfun(@(x) any(isnan(Mdata(:,x))), 1:roiSize);
                numericColumns=find(~nanColumns);
                len=length(numericColumns);
                %[coeff, score, latent, tsquare] = princomp(x,econFlag)
                if(len>1)
                    %The extracted data that do not contains Nan at any rows
                    Mdata1=zeros(dataSize,len);
                    
                    for col=1:len
                        Mdata1(:,col)=Mdata(:,numericColumns(col));
                    end
                    [coeff, score] = princomp(Mdata1);
                    PC=score(:,1);
                    d=sign(sum(coeff(:,1)));
                    PC=d*PC;
                else
                    PC=nan*ones(dataSize);
                end
                
                for row=1:dataSize
                    p.PC_act(row,r)=PC(row);
                end
                
                toc
            end
        end
        
        function vols = getVoxelValues(roi_imgs, data_imgs,hold)
            %Taihao
            %This function returns the voxel values of data_imgs at the
            %voxel positions of roi_imgs
            %the output vols is a 2d cellarray. The cell at i row and j
            %colunm is the voxel value of the j-th data_img at the i-th
            %roi_img
            if(~iscell(roi_imgs))
                roi_imgs={roi_imgs};
            end
            
            if(~iscell(data_imgs))
                data_imgs={data_imgs};
            end
            
            vols=cell(length(roi_imgs), length(data_imgs));
            
            if(~exist('hold','var'))
                hold = 1; % trilinear
            end
            
            % hold   -  sets the interpolation method for the resampling.
            %           0          Zero-order hold (nearest neighbour).
            %           1          First-order hold (trilinear interpolation).
            %           2->127     Higher order Lagrange (polynomial) interpolation using
            %                      different holds (second-order upwards).
            %          -127 - -1   Different orders of sinc interpolation.
            % X      -  output image
            
            % load the default spm settings
            spm('defaults', 'fmri');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % use SPM functions to read image information
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if(iscell(roi_imgs))
                roi_imgs=char(roi_imgs);
            end
            
            if(iscell(data_imgs))
                data_imgs=char(data_imgs);
            end
            
            v_roi=spm_vol(roi_imgs);
            v_data=spm_vol(data_imgs);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % loop through the ROIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for r = 1:length(v_roi)
                tic
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % xyz mm of non-zero points in roi img
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [vol,XYZ] = spm_read_vols(spm_vol(v_roi(r)));
                idx_roi     = find(vol);
                XYZ     = [XYZ(:, idx_roi); ones(1, length(idx_roi))];
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % loop through the data images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                for i = 1:length(v_data)
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % find voxel elements for which to extract stats
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    ixyz = v_data(i).mat \ XYZ;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % then sample the data from the volume at those locations
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    data = spm_sample_vol(v_data(i),ixyz(1,:),ixyz(2,:),ixyz(3,:),hold);
                    voli=vol;
                    voli(idx_roi)=data;
                    vols{r,i}=voli;
                end
            end
        end
        
        function vols=getVoxelValues_coors(XYZ,data_imgs,hold)
            %XYZ: the coordinates, each column is a x, y, z coordinates
            %data_imgs, a cellarray of data images
            
            % hold   -  sets the interpolation method for the resampling.
            %           0          Zero-order hold (nearest neighbour).
            %           1          First-order hold (trilinear interpolation).
            %           2->127     Higher order Lagrange (polynomial) interpolation using
            %                      different holds (second-order upwards).
            %          -127 - -1   Different orders of sinc interpolation.
            %
            %vols: voxel values
            if(~iscell(data_imgs))
                data_imgs={data_imgs};
            end
            
            vols=cell(size(XYZ,2), length(data_imgs));
            
            if(~exist('hold','var'))
                hold = 1; % trilinear
            end
            v_data=spm_vol(data_imgs);
            
            for i = 1:length(v_data)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find voxel elements for which to extract stats
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                ixyz = v_data(i).mat \ XYZ;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % then sample the data from the volume at those locations
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                data = spm_sample_vol(v_data(i),ixyz(1,:),ixyz(2,:),ixyz(3,:),hold);
                voli=vol;
                voli(idx_roi)=data;
                vols{r,i}=voli;
            end
        end
        
        function fdir = getFileDir(fpath)
            %this function returns the directory part of "fpath".
            parts=strsplit(fpath, filesep);
            fdir=[filesep parts{1}];
            for i=1:length(parts)-1
                fdir=fullfile(fdir,parts{i});
            end
        end
        
        function fname = getFileName(fpath)
            %this function returns the file name of "fpath".
            fname=[];
            if(~isempty(fpath))
                [path, fname, ext] = fileparts(fpath);
            end
        end
        
        function fnames = getFileNames(fpaths)
            %this function returns the file name of "fpath".
            fnames=cellfun(@(x) CommonMethods.getFileName(x), fpaths, 'UniformOutput', false);
        end
        
        function ext = getFileExtension(fpath)
            %this function returns the file extension of "fpath".
            parts=strsplit(fpath,filesep);
            parts=strsplit(parts{end}, '.');
            len=length(parts);
            if(len==1)
                ext='';
                return;
            end
            ext=parts{len};
        end
        
        function path = setFileExtension(fpath,ext)
            %this function sets file extension as ext.
            indexes=strfind(fpath,'.');
            if(isempty(indexes))
                path=[fpath '.' ext];
            end
            ind=indexes(end);
            path=[fpath(1:ind) ext];
        end
        
        function convertImag(infile,outfile)
            vol=spm_vol(infile);
            V=spm_read_vols(vol);
            vol.fname=outfile;
            spm_write_vols(V,vol);
        end
        
        function sizes=buildGroupAcompCorNoiseMask()
            %this function build group acompCor Noise Mask
            %it reads in aCompCor noise mask files of all scan in a group,
            %to create group mask according to the specified cutoff value.
            
            ParadigmName='EMOID';
            InfoOrganizer=PATHD_getInfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            scanIDs=InfoOrganizer.getScanIDs_scanLabel('a');
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','swmaskU1.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            
            sizes=zeros(1,length(maskImgs));
            sizes(1)=nnz(A);
            for s=2:length(scanIDs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                sizes(s)=nnz(A0);
                A=A+A0;
                u=unique(A0);
                cprintf('blue','%s\n', [scanIDs{s} 'number of unique elements: ' CommonMethods.array2str_CSV(u)]);
                if(length(u)>2)
                    s=s;
                end
            end
            B=A>0.99999*length(scanIDs);
            C=A>0.000001*length(scanIDs);
            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'noise1Maks_Intersection.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'noise1Maks_Union.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'noise1MaksSum.nii');
            spm_write_vol(V,A);
            
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','swmaskU2.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            for s=2:length(scanIDs)
                A=A+spm_read_vols(spm_vol(maskImgs{s}));
                u=unique(A0);
                cprintf('blue','%s\n', [scanIDs{s} 'number of unique elements: ' CommonMethods.array2str_CSV(u)]);
                if(length(u)>2)
                    s=s;
                end
            end
            B=A>0.99999*length(scanIDs);
            C=A>0.000001*length(scanIDs);
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'noise2Maks_Intersection.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'noise2Maks_Union.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'noise2MaksSum.nii');
            spm_write_vol(V,A);
            
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU1.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            for s=2:length(scanIDs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                A=A+A0;
                u=unique(A0);
                cprintf('blue','%s\n', [scanIDs{s} 'number of unique elements: ' CommonMethods.array2str_CSV(u)]);
                if(length(u)>2)
                    s=s;
                end
            end
            B=A>0.99999*length(scanIDs);
            C=A>0.00001*length(scanIDs);
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'noise1Maks_Intersection_Unsmoothed.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'noise1Maks_Union_Unsmoothed.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'noise1MaksSum_Unsmoothed.nii');
            spm_write_vol(V,A);
            
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU2.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            for s=2:length(scanIDs)
                A=A+spm_read_vols(spm_vol(maskImgs{s}));
                u=unique(A0);
                cprintf('blue','%s\n', [scanIDs{s} 'number of unique elements: ' CommonMethods.array2str_CSV(u)]);
                if(length(u)>2)
                    s=s;
                end
            end
            B=A>0.99999*length(scanIDs);
            C=A>0.00001*length(scanIDs);
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'noise2Maks_Intersection_Unsmoothed.nii');
            spm_write_vol(V,A);
            V.fname = fullfile(outdir,'noise2Maks_Union_Unsmoothed.nii');
            spm_write_vol(V,A);
            V.fname = fullfile(outdir,'noise2MaksSum_Unsmoothed.nii');
            spm_write_vol(V,A);
        end
        
        function sizes=buildGroupeSPMMask()
            %this function build group acompCor Noise Mask
            %it reads in aCompCor noise mask files of all scan in a group,
            %to create group mask according to the specified cutoff value.
            
            ParadigmName='EMOID';
            InfoOrganizer=PATHD_getInfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            scanIDs=InfoOrganizer.getScanIDs_scanLabel('a');
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'stats_arf_aCompCor','swmaskSPM.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            sizes=zeros(1,length(maskImgs));
            sizes(1)=nnz(A);
            for s=2:length(scanIDs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                sizes(s)=nnz(A0);
                A=A+A0;
            end
            B=A>0.99999*length(scanIDs);
            C=A>0.00001*length(scanIDs);
            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'MaksSPM_Intersection.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'MaksSPM_Union.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'MaksSPM_Sum.nii');
            spm_write_vol(V,A);
            
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'stats_arf_aCompCor','wmaskSPM.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            for s=2:length(scanIDs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                A=A+A0;
                u=unique(A0);
                cprintf('blue','%s\n', [scanIDs{s} 'number of unique elements: ' CommonMethods.array2str_CSV(u)]);
                if(length(u)>2)
                    s=s;
                end
            end
            B=A>0.999*length(scanIDs);
            C=A>0.001*length(scanIDs);
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'MaksSPM_Intersection_Unsmoothed.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'MaksSPM_Union_Unsmoothed.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'MaksSPMSum_Unsmoothed.nii');
            spm_write_vol(V,A);
        end
        
        function sizes=buildGroupAcompCorMask()
            %this function build a group mask file by excluding the
            %aCompCor noise mask from the SPM mask.
            %the hard coded parameters should be generalize if this type of
            %tasks need to be performed frequently.
            ParadigmName='EMOID';
            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            InfoOrganizer=PATHD_getInfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            %            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU1.nii'),scanIDs,'UniformOutput', false);
            scanIDs=InfoOrganizer.getScanIDs_scanLabel('a');
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'stats_arf_aCompCor','swmaskSPM-aCompCor1.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            sizes=zeros(1,length(maskImgs));
            sizes(1)=nnz(A);
            for s=2:length(scanIDs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                sizes(s)=nnz(A0);
                A=A+A0;
            end
            %             A=A/length(scanIDs);
            %             name='maskSPM-aCompCor1_';
            %             for i=80:2:100
            %                 B=A>i*0.01;
            %                 V.fname = fullfile(outdir,[name num2str(i) '.nii']);
            %                 spm_write_vol(V,B);
            %             end
            B=A>0.9999*length(scanIDs);
            C=A>0.0001*length(scanIDs);
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'maskSPM-aCompCor_Intersection.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'maskSPM-aCompCor_Union.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'maskSPM-aCompCor_Sum.nii');
            spm_write_vol(V,A);
        end
        
        function sizes=buildGroupFreeSurferMask()
            %this function build a group mask file by excluding the
            %aCompCor noise mask from the SPM mask.
            %the hard coded parameters should be generalize if this type of
            %tasks need to be performed frequently.
            ParadigmName='EMOID';
            %            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            outdir='/home5/data/process/PTST/GroupAnalysis/GroupMasks/EMOID';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            InfoOrganizer=getPTSD_InfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            %            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU1.nii'),scanIDs,'UniformOutput', false);
            scanIDs=InfoOrganizer.AnalysisInfo.ScanIDs;
            maskImgs=cellfun(@(x) fullfile(dataDir,x,x,ParadigmName,'1','analysis',['wr' x '_Amygdala_whole.nii']),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            A=A~=0;
            sizes=zeros(1,length(maskImgs));
            sizes(1)=nnz(A);
            As=cell(1,length(scanIDs));
            len=length(scanIDs);
            As{1}=A;
            for s=2:length(scanIDs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                A0=A0~=0;
                sizes(s)=nnz(A0);
                A=A+A0;
                As{s}=A0;
            end
            %             A=A/length(scanIDs);
            %             name='maskSPM-aCompCor1_';
            %             for i=80:2:100
            %                 B=A>i*0.01;
            %                 V.fname = fullfile(outdir,[name num2str(i) '.nii']);
            %                 spm_write_vol(V,B);
            %             end
            overlapping=ones(1,len*(len-1)/2);
            JaccardIndex=ones(1,len*(len-1)/2);
            noPairs=len*(len-1)/2;
            pairIndexes=ones(noPairs,2);
            n=0;
            for i=1:len-1
                Ai=As{i};
                nnzi=sizes(i);
                for j=i+1:len
                    n=n+1;
                    Aj=As{j};
                    nnzj=sizes(j);
                    Aij=Ai.*Aj;
                    nnzij=nnz(Aij);
                    overlapping(n)=nnzij/min(nnzi,nnzj);
                    JaccardIndex(n)=nnzij/nnz(Ai+Aj);
                    pairIndexes(n,:)=[i j];
                end
            end
            [OFS, PIS]=sort(overlapping);
            
            figure();
            hist(JaccardIndex);
            title(gca,'Pairwise Jaccard Index','FontSize', 24);
            set(gca,'FontSize',18);
            
            figure();
            hist(overlapping);
            title(gca,'Pairwise Overlapping Fraction','FontSize', 24);
            set(gca,'FontSize',18);
            %%
            
            B=A==length(scanIDs);
            C=A~=0;
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'maskAmygdala_whole_Intersection.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'maskAmygdala_whole_Union.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'Amygdala_whole_Sum.nii');
            spm_write_vol(V,A);
            
            fname='pairwiseOverlapping_Amygdala_whole.cvs';
            columnNames={'scanID1' 'scanID2' 'size1' 'size2' 'OverlappingFraction' 'JaccardIndex'};
            rowNames={};
            table=cell(noPairs,length(columnNames));
            for n0=1:noPairs
                n=PIS(n0);
                i=pairIndexes(n,1);
                j=pairIndexes(n,2);
                if(n0==2||n0==(noPairs-3))
                    O=CommonMethods.extractMaskOverlapping(As{i},As{j});
                    CommonMethods.displayVolAsImage(O,3,1);
                    v1=V;
                    fnamev=['overlapping_' scanIDs{i} '_' scanIDs{j} '.nii'];
                    v1.fname=fullfile(outdir,fnamev);
                    v1.dim=size(O);
                    spm_write_vol(v1,O);
                    figFile=fullfile(outdir, CommonMethods.setFileExtension(fnamev,'tif'));
                    print ('-dtiff', figFile);
                end
                table(n0,:)={scanIDs(i) scanIDs(j) num2str(sizes(i)) num2str(sizes(j)) num2str(overlapping(n)) num2str(JaccardIndex(n))};
            end
            CommonMethods.exportTable_CSV(fullfile(outdir,fname),columnNames,rowNames,table);
        end
        
        
        function sizes=displaySPMMaskOverlapping()
            %this function build a group mask file by excluding the
            %aCompCor noise mask from the SPM mask.
            %the hard coded parameters should be generalize if this type of
            %tasks need to be performed frequently.
            ParadigmName='EMOID';
            %            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            InfoOrganizer=PATHD_getInfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            %            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU1.nii'),scanIDs,'UniformOutput', false);
            scanIDs=InfoOrganizer.AnalysisInfo.ScanIDs;
            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'stats_arf_aCompCor','wmaskSPM.nii'),scanIDs,'UniformOutput', false);
            scanIDs0=scanIDs;
            maskImgs0=maskImgs;
            maskImgs={};
            scanIDs={};
            for i=1:length(scanIDs0)
                mask=maskImgs0{i};
                if(~exist(mask,'file'))
                    continue;
                end
                scanIDs{end+1}=scanIDs0{i};
                maskImgs{end+1}=maskImgs0{i};
            end
            
            V=spm_vol(maskImgs{1});
            A=spm_read_vols(V);
            A=A~=0;
            sizes=zeros(1,length(maskImgs));
            sizes(1)=nnz(A);
            As=cell(1,length(scanIDs));
            len=length(scanIDs);
            As{1}=A;
            for s=2:length(maskImgs)
                A0=spm_read_vols(spm_vol(maskImgs{s}));
                A0=A0~=0;
                sizes(s)=nnz(A0);
                A=A+A0;
                As{s}=A0;
            end
            %             A=A/length(scanIDs);
            %             name='maskSPM-aCompCor1_';
            %             for i=80:2:100
            %                 B=A>i*0.01;
            %                 V.fname = fullfile(outdir,[name num2str(i) '.nii']);
            %                 spm_write_vol(V,B);
            %             end
            overlapping=ones(1,len*(len-1)/2);
            JaccardIndex=ones(1,len*(len-1)/2);
            noPairs=len*(len-1)/2;
            pairIndexes=ones(noPairs,2);
            n=0;
            for i=1:len-1
                Ai=As{i};
                nnzi=sizes(i);
                for j=i+1:len
                    n=n+1;
                    Aj=As{j};
                    nnzj=sizes(j);
                    Aij=Ai.*Aj;
                    nnzij=nnz(Aij);
                    overlapping(n)=nnzij/min(nnzi,nnzj);
                    JaccardIndex(n)=nnzij/nnz(Ai+Aj);
                    pairIndexes(n,:)=[i j];
                end
            end
            [OFS, PIS]=sort(overlapping);
            
            figure();
            hist(JaccardIndex);
            title(gca,'Pairwise Jaccard Index','FontSize', 24);
            set(gca,'FontSize',18);
            
            figure();
            hist(overlapping);
            title(gca,'Pairwise Overlapping Fraction','FontSize', 24);
            set(gca,'FontSize',18);
            %%
            
            B=A==length(scanIDs);
            C=A~=0;
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            V.fname = fullfile(outdir,'maskSPM_Intersection.nii');
            spm_write_vol(V,B);
            V.fname = fullfile(outdir,'maskSPM_Union.nii');
            spm_write_vol(V,C);
            V.fname = fullfile(outdir,'maskSPM_Sum.nii');
            spm_write_vol(V,A);
            
            fname='pairwiseOverlapping_SPM.cvs';
            columnNames={'scanID1' 'scanID2' 'size1' 'size2' 'OverlappingFraction' 'JaccardIndex'};
            rowNames={};
            table=cell(noPairs,length(columnNames));
            for n0=1:noPairs
                n=PIS(n0);
                i=pairIndexes(n,1);
                j=pairIndexes(n,2);
                if(n0==2||n0==(noPairs-3))
                    O=CommonMethods.extractMaskOverlapping(As{i},As{j});
                    CommonMethods.displayVolAsImage(O,3,1);
                    v1=V;
                    fnamev=['overlapping_' scanIDs{i} '_' scanIDs{j} '.nii'];
                    v1.fname=fullfile(outdir,fnamev);
                    v1.dim=size(O);
                    spm_write_vol(v1,O);
                    figFile=fullfile(outdir, CommonMethods.setFileExtension(fnamev,'tif'));
                    print ('-dtiff', figFile);
                end
                table(n0,:)={scanIDs(i) scanIDs(j) num2str(sizes(i)) num2str(sizes(j)) num2str(overlapping(n)) num2str(JaccardIndex(n))};
            end
            CommonMethods.exportTable_CSV(fullfile(outdir,fname),columnNames,rowNames,table);
        end
        
        function sizes=compairAmygdalaMasks()
            %this function build a group mask file by excluding the
            %aCompCor noise mask from the SPM mask.
            %the hard coded parameters should be generalize if this type of
            %tasks need to be performed frequently.
            ParadigmName='EMOID';
            %            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            outdir='/home5/data/process/PTST/GroupAnalysis/GroupMasks/EMOID';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            InfoOrganizer=getPTSD_InfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            %            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU1.nii'),scanIDs,'UniformOutput', false);
            scanIDs=InfoOrganizer.AnalysisInfo.ScanIDs;
            
            maskImgsFS=cellfun(@(x) fullfile(dataDir,x,x,ParadigmName,'1','analysis',['r' x '_Amygdala_whole.nii']),scanIDs,'UniformOutput', false);
            maskImgsRN=cellfun(@(x) fullfile(dataDir,x,x,ParadigmName,'temp_gPPI','gPPI__ROI_B_Amygdala_mask.nii'),scanIDs,'UniformOutput', false);
            V=spm_vol(maskImgsFS{1});
            As={};
            Bs={};
            sizesA={};
            sizesB={};
            maskImgsFS0=maskImgsFS;
            maskImgsRN0=maskImgsRN;
            maskImgsFS={};
            maskImgsRN={};
            scanIDs0=scanIDs;
            scanIDs={};
            for s=1:length(scanIDs0)
                fB=maskImgsRN0{s};
                if(~exist(fB,'file'))
                    continue;
                end
                scanIDs{end+1}=scanIDs0{s};
                maskImgsFS{end+1}=maskImgsFS0{s};
                maskImgsRN{end+1}=maskImgsRN0{s};
            end
            
            %          ResliceImages_job(maskImgsFS(1),maskImgsRN)
            
            for s=1:length(scanIDs)
                fA=maskImgsFS{s};
                fB=maskImgsRN{s};
                cprintf('blue','%s\n',['ScanID: ' scanIDs{s}]);
                A=spm_read_vols(spm_vol(fA));
                B=spm_read_vols(spm_vol(fB));
                %               B=CommonMethods.resliceVolume(B,size(A));
                A=A~=0;
                B=B~=0;
                As{end+1}=A;
                Bs{end+1}=B;
                sizesA{end+1}=nnz(A);
                sizesB{end+1}=nnz(B);
            end
            
            len=length(As);
            overlapping=ones(1,len)*(-1);
            JaccardIndex=ones(1,len)*(-1);
            
            for i=1:len
                A=As{i};
                nnzA=sizesA{i};
                B=Bs{i};
                nnzB=sizesB{i};
                AB=A.*B;
                nnzAB=nnz(AB);
                overlapping(i)=nnzAB/min(nnzA,nnzB);
                JaccardIndex(i)=nnzAB/nnz(A+B);
            end
            [OFS, PIS]=sort(overlapping);
            
            figure();
            hist(JaccardIndex);
            title(gca,'Jaccard Index','FontSize', 24);
            set(gca,'FontSize',18);
            
            figure();
            hist(overlapping);
            title(gca,'Overlapping Fraction','FontSize', 24);
            set(gca,'FontSize',18);
            %%
            
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            
            fname='Overlapping_Amygdala_whole_RNandFS.cvs';
            columnNames={'scanID1' 'size_FS' 'size_RN' 'OverlappingFraction' 'JaccardIndex'};
            rowNames={};
            table=cell(len,length(columnNames));
            for n0=1:len
                n=PIS(n0);
                if(n0==2||n0==(len-3))
                    O=CommonMethods.extractMaskOverlapping(As{n},Bs{n});
                    CommonMethods.displayVolAsImage(O,3,1);
                    v1=V;
                    fnamev=['overlapping_' scanIDs{i} '_FS_RN.nii'];
                    v1.fname=fullfile(outdir,fnamev);
                    v1.dim=size(O);
                    spm_write_vol(v1,O);
                    figFile=fullfile(outdir, CommonMethods.setFileExtension(fnamev,'tif'));
                    print ('-dtiff', figFile);
                end
                table(n0,:)={scanIDs{n} num2str(sizesA{n}) num2str(sizesB{n}) num2str(overlapping(n)) num2str(JaccardIndex(n))};
            end
            CommonMethods.exportTable_CSV(fullfile(outdir,fname),columnNames,rowNames,table);
        end
        
        function sizes=compairAmygdalaMasks_ExtractROIs()
            %this function build a group mask file by excluding the
            %aCompCor noise mask from the SPM mask.
            %the hard coded parameters should be generalize if this type of
            %tasks need to be performed frequently.
            ParadigmName='EMOID';
            %            outdir='/home2/data/process/PATHD/GroupAnalysis/GroupMasks/EMOID';
            outdir='/home5/data/process/PTSD/GroupAnalysis/GroupMasks/EMOID/ROIs';
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            InfoOrganizer=getPTSD_InfoOrganizer(ParadigmName);
            %             groups=PATHD_Grouping.buildGroups(ParadigmName);
            %             scanIDs=groups{1}.ScanIDs;
            dataDir=InfoOrganizer.AnalysisInfo.dataDir;
            %            maskImgs=cellfun(@(x) fullfile(dataDir,x(1:5),x,ParadigmName,'1','analysis','wmaskU1.nii'),scanIDs,'UniformOutput', false);
            scanIDs=InfoOrganizer.AnalysisInfo.ScanIDs;
            
            maskImgsFS=cellfun(@(x) fullfile(dataDir,x,x,ParadigmName,'1','analysis',['r' x '_Amygdala_whole.nii']),scanIDs,'UniformOutput', false);
            maskImgsRN=cellfun(@(x) fullfile(dataDir,x,x,ParadigmName,'temp_gPPI','gPPI__ROI_B_Amygdala_mask.nii'),scanIDs,'UniformOutput', false);
            dataImgs=cellfun(@(x) fullfile(dataDir,x,x,ParadigmName,'1','analysis','arf.nii'),scanIDs,'UniformOutput', false);
            
            maskImgsFS0=maskImgsFS;
            maskImgsRN0=maskImgsRN;
            dataImgs0=dataImgs;
            dataImgs={};
            maskImgsFS={};
            maskImgsRN={};
            scanIDs0=scanIDs;
            scanIDs={};
            for s=1:length(scanIDs0)
                fB=maskImgsRN0{s};
                if(~exist(fB,'file'))
                    continue;
                end
                scanIDs{end+1}=scanIDs0{s};
                maskImgsFS{end+1}=maskImgsFS0{s};
                maskImgsRN{end+1}=maskImgsRN0{s};
                dataImgs{end+1}=dataImgs0{s};
            end
            
            %             roi_imgs={'/home2/data/process/PATHD/GroupAnalysis/clusterImgs/LateralVentricle_from_spmT_0001.img'};
            %             line=[data_imgs{i} ', ' num2str(stats.max_act(i,1)) ', ' num2str(stats.min_act(i,1)) ', ' num2str(stats.mean_act(i,1)) ', ' num2str(stats.median_act(i,1)) ];
            
            
            len=length(maskImgsFS);
            fname='Amygdala_ROIs';
            columnNames={'ScanID'};
            
            rowNames={};
            noVols=171;
            noCons=12;
            
            for i=1:noCons
                columnNames{end+1}=['con' num2str(i) '_meanFS'];
                columnNames{end+1}=['con' num2str(i) '_meanRN'];
            end
            table=cell(noVols,length(columnNames));
            
            for i=1:len
                scanID=scanIDs{i};
                cprintf('blue','%s\n',scanIDs{i});
                rois={char(maskImgsFS{i}) char(maskImgsRN{i})};
                
                %                 datas={char(dataImgs{i})};
                dir=fullfile('/home/tle6/data/process/PTSD/',scanID, scanID, 'EMOID', 'stats_arf_aCompCor');
                datas=arrayfun(@(x) fullfile(dir,['con_' sprintf('%04d.img',x)]), 1:noCons, 'UniformOutput',false);
                stats = CommonMethods.statExtraction(char(rois),char(datas));
                
                table{i,1}=scanID;
                for j=1:noCons
                    table{i,2*j}=num2str(stats.mean_act(j,1));
                    table{i,2*j+1}=num2str(stats.mean_act(j,2));
                end
                %
                %                 for j=1:noVols
                %                     table{j,2*(i-1)+1}=num2str(stats.mean_act(j,1));
                %                     table{j,2*i}=num2str(stats.mean_act(j,2));
                %                 end
            end
            CommonMethods.exportTable_CSV(fullfile(outdir,fname),columnNames,rowNames,table);
            
        end
        
        function resliceImages(imagesToReslice, ReferenceImage, interp)
            %imagesToSlice will be resliced to the same dimension as the
            %Reference Image reslicedImages =
            template=fullfile(CommonEnvMethods.getScriptPath(),'Jobtemplates','ResliceImages_trilinear.mat');
            load (template);
            if(~exist('interp', 'var'))
                %interp: the interpolation option parameter: 0 for the
                %nearest neighbor, 1 for trilnear, 2, 3 ... are for the 2, 3 ... degree of
                %B spline.
                interp=1;
            end
            
            if(ischar(imagesToReslice))
                imagesToReslice={imagesToReslice}
            end
            
            matlabbatch{1,1}.spm.spatial.coreg.write.ref={ReferenceImage};
            matlabbatch{1,1}.spm.spatial.coreg.write.source=imagesToReslice;
            matlabbatch{1,1}.spm.spatial.coreg.write.roptions.interp=interp;
            try
                spm_jobman('run', matlabbatch);
                fprintf('%s\n','DONE');
            catch err
                fprintf('%s\n',['ERROR with reslicing job: ' err.message]);
            end
        end
               
        function changeVoxelSize(img, voxelSize, interp)
            if ~exist('interp','var')
                interp=1;
            end
            voxsiz = voxelSize; % new voxel size {mm}
            V = img;
            V = spm_vol(V);
            for i=1:numel(V)
                bb        = spm_get_bbox(V(i));
                VV(1:2)   = V(i);
                VV(1).mat = spm_matrix([bb(1,:) 0 0 0 voxsiz])*spm_matrix([-1 -1 -1]);
                VV(1).dim = ceil(VV(1).mat \ [bb(2,:) 1]' - 0.1)';
                VV(1).dim = VV(1).dim(1:3);
                spm_reslice(VV,struct('mean',false,'which',1,'interp',interp)); % 1 for linear
            end            
        end
        
        function resliceImages_save(imagesToReslice, ReferenceImage, interp, newImageNames)
            %imagesToSlice will be resliced to the same dimension as the
            %Reference Image reslicedImages =
            
            if(ischar(imagesToReslice))
                imagesToReslice={imagesToReslice};
            end
            
            if(ischar(newImageNames))
                newImageNames={newImageNames};
            end
            
            CommonMethods.resliceImages(imagesToReslice, ReferenceImage, interp);
            if(exist('newImageNames','var'))
                if(length(imagesToReslice)==length(newImageNames))
                    rnames=CommonMethods.getPrefixedFiles(imagesToReslice, 'r');
                    for i=1:length(rnames)
                        rname=rnames{i};
                        V=spm_vol(rname);
                        vols=spm_read_vols(V);
                        V.fname=newImageNames{i};
                        spm_write_vol(V,vols);
                        delete(rname);
                    end
                else
                    error('the number of images to reslice and new image names have to be the same');
                end
            end
        end
        
        function [V, origin] = extractMaskOverlapping(A,B)
            %this function extract the submatrix that encloses A+B~=0
            C=A+B*2;
            dim=size(C);
            indexes=find(C);
            [x,y,z]=Ind2sub(dim,indexes);
            xn=min(x);
            xx=max(x);
            yn=min(y);
            yx=max(y);
            zn=min(z);
            zx=max(z);
            origin=[xn yn zn];
            V=C(xn:xx,yn:yx,zn:zx);
        end
        
        function displayVolAsImage(V, axis, increament)
            %this function displays V as a series of images along axis, as well
            %as displays the tessellated images.
            dim=size(V,axis);
            AxisLabel='XYZ';
            label=AxisLabel(axis)
            len=floor(dim/increament);
            ms=cell(1,len);
            n=0;
            for i=1:increament:dim
                n=n+1;
                switch axis
                    case 1
                        m=V(i,:,:);
                    case 2
                        m=V(:,i,:);
                    case 3
                        m=V(:,:,i);
                end
                %                 figure();
                %                 title([label ' = ' num2str(i)]);
                %                 imagesc(m);
                ms{n}=m;
            end
            M=CommonMethods.tessellateMatrices(ms,2,-2);
            figure();
            imagesc(M);
        end
        
        function range=getMatValueRange(m)
            range=CommonMethods.getRange(m);
        end
        
        
        function M = tessellateMatrices(ms,spacer,base)
            len=length(ms);
            sq=sqrt(len);
            w=ceil(sq);
            h=floor(sq);
            if(w*h<len)
                h=w;
            end
            
            dim=size(ms{1});
            w0=dim(2);
            h0=dim(1);
            M=ones(h*(h0+spacer)-spacer,w*(w0+spacer)-spacer)*base;
            
            for i=1:len
                row=floor((i-1)/w)+1;
                col=i-(row-1)*w;
                x0=(col-1)*(w0+spacer);
                y0=(row-1)*(h0+spacer);
                m=ms{i};
                M(y0+1:y0+h0,x0+1:x0+w0)=m(:,:);
            end
        end
        
        function paths = getExistingFilesOnly(paths0)
            %returns existing files in paths0
            paths={};
            for i=1:length(paths0)
                if(exist(paths0{i},'file'))
                    paths{end+1}=paths0{i};
                end
            end
        end
        
        function showFileNamesInSPM()
            %load SPM and display the file names
            load('SPM.mat');
            files={SPM.xY.VY.fname};
            celldisp(files);
        end
        
        function qacell = calDVARS(procfile,maskfile,aCompCor_file,motionfile,noVols, repetitionTime, bandpass_cutoffs)
            %this function compute DVARS of the image "procfile" and return
            %the results in qacell.
            %the calculation of DVARS is based on voles in "maskfile".
            %Only those voxels in "maskfile" are included in the
            %calculation of DVARS
            %there are 6 version of DVARS computed here, according to the
            %images and the way of normalization.
            %the images are original "procfile", after perform regression
            %only, and after perform regression and band pass filtering.
            
            %the normalization is either per frame or by the mode of the
            %image.
            
            qacell={};
            
            %%% Create Mask
            if ~exist(maskfile, 'file')
                cprintf([0.6 0.3 0.1], 'no mask file\n');
                %                createMask(procfile, 1, ['mask_' procPrefix '.img']);
                return;
            end
            
            %%% Calculate DVARS
            %         if exist(maskfile, 'file') && (~exist(pctdvars_file, 'file') || rerunAnalysis)
            
            mask                = spm_read_vols(spm_vol(maskfile));
            
            motion              = load( motionfile);
            motion_deriv       	= [0 0 0 0 0 0; diff(motion)];
            
            aCompCor            = load(aCompCor_file);
            aCompCor            = aCompCor.aCompCor_stats.pca(:,1:aCompCor.aCompCor_stats.numComp_MonteCarlo);
            
            
            %get the data in shape:
            data                = reshape(spm_read_vols(spm_vol(procfile)), numel(mask), noVols)';
            
            
            %find the voxels in the mask:
            brainvox_idx        = find(mask);
            
            %keep original data
            t=data(:, brainvox_idx);;
            dim=size(t);
            t=reshape(t,[1,dim(1)*dim(2)]);
            m=mode(t);
            
            data_no=1000*(data/m);
            data_n=data_no;
            %std iqr
            %build the design matrix
            X                   = [ones(noVols,1) motion motion_deriv aCompCor];
            
            %run the regression:
            betas              	= arrayfun(@(v) regress(data(:, v), X), 1:size(data,2), 'UniformOutput', false);
            betas               = [betas{:}];
            
            betas_n              	= arrayfun(@(v) regress(data_n(:, v), X), 1:size(data_n,2), 'UniformOutput', false);
            betas_n               = [betas_n{:}];
            %calculate residuals:
            data                = data - (X(:, 2:end) * betas(2:end, :));
            data_n                = data_n - (X(:, 2:end) * betas_n(2:end, :));
            
            %Previously, filter was applied as well:
            fdata= conn_filter(repetitionTime, bandpass_cutoffs, data);
            fdata_n= conn_filter(repetitionTime, bandpass_cutoffs, data_n);
            fdata_no= conn_filter(repetitionTime, bandpass_cutoffs, data_no);
            
            
            %reduce data to the voxels in the mask:
            data = data(:, brainvox_idx);
            fdata= fdata(:, brainvox_idx);
            fdata_n= fdata_n(:, brainvox_idx);
            fdata_no= fdata_no(:, brainvox_idx);
            data_n = data_n(:, brainvox_idx);
            data_no = data_no(:, brainvox_idx);
            
            %DVARS of regressed data, the temporal derivative of a voxel is
            %normalized by the intensity of the voxel.
            pctdiff=diff(data);
            %the time series of DVARS
            pctdvars            = sqrt(mean(((pctdiff ./ data(1:end-1,:))) .^ 2, 2));
            qacell{end+1}=pctdvars;
            noOutliers= nnz(pctdvars > 0.005);
            mdDvars=median(pctdvars);
            %the median of DVARS
            qacell{end+1} = mdDvars;
            %the number of outliers (the frames with DVARS > 0.005)
            qacell{end+1}=noOutliers;
            
            
            %DVARS of regressed and filtered data, the temporal derivative of a voxel is
            %normalized by the intensity of the voxel.
            pctdvars_f = sqrt(mean(((diff(fdata) ./ fdata(1:end-1,:))) .^ 2, 2));
            qacell{end+1}=pctdvars_f;
            noOutliers_f= nnz(pctdvars_f > 0.005);
            mdDvars_f=median(pctdvars_f);
            qacell{end+1} = mdDvars_f;
            qacell{end+1}=noOutliers_f;
            
            %DVARS of normalized, regressed data. The normalization is
            %performed by scaling the entire image so that the mode value
            %the original data is 1000.
            pctdvars_n = sqrt(mean(((diff(data_n))) .^ 2, 2));
            qacell{end+1}=pctdvars_n;
            noOutliers_n= nnz(pctdvars_n > 5);
            mdDvars_n=median(pctdvars_n);
            qacell{end+1} = mdDvars_n;
            qacell{end+1}=noOutliers_n;
            
            %DVARS of normalized, regressed, and filted data.
            pctdvars_fn = sqrt(mean(((diff(fdata_n) )) .^ 2, 2));
            qacell{end+1}=pctdvars_fn;
            noOutliers_fn= nnz(pctdvars_fn > 5);
            mdDvars_fn=median(pctdvars_fn);
            qacell{end+1} = mdDvars_fn;
            qacell{end+1}=noOutliers_fn;
            
            %DVARS of normalized original data.
            pctdvars_no = sqrt(mean(((diff(data_no) )) .^ 2, 2));
            qacell{end+1}=pctdvars_no;
            noOutliers_no= nnz(pctdvars_no > 5);
            mdDvars_no=median(pctdvars_no);
            qacell{end+1} = mdDvars_no;
            qacell{end+1}=noOutliers_no;
            
            %DVARS of normalized and filtered data.
            pctdvars_nof = sqrt(mean(((diff(fdata_no) )) .^ 2, 2));
            qacell{end+1}=pctdvars_nof;
            noOutliers_nof= nnz(pctdvars_nof > 5);
            mdDvars_nof=median(pctdvars_nof);
            qacell{end+1} = mdDvars_nof;
            qacell{end+1}=noOutliers_nof;
        end
        
        
        function fdata = FCAProcessing(data, regressor, noVols, repetitionTime, bandpass_cutoffs)
            %this function performs functional connectivity processing to
            %the data passed in. data is a 4d matrix, the 4th dimension is
            %time.
            
            data                = reshape(spm_read_vols(spm_vol(procfile)), numel(mask), noVols)';
            X                   = [ones(noVols,1) regressor];
            
            %run the regression:
            betas              	= arrayfun(@(v) regress(data(:, v), X), 1:size(data,2), 'UniformOutput', false);
            betas               = [betas{:}];
            
            %calculate residuals:
            data                = data - (X(:, 2:end) * betas(2:end, :));
            
            %Previously, filter was applied as well:
            fdata= conn_filter(repetitionTime, bandpass_cutoffs, data);
        end
        
        function pars=calMovementPars(motionfile)
            rp                 = dlmread(motionfile); %motion file should be the realignment parameter text file
            rpdiff              = diff(rp);
            
            %%% Calculate framewise displacement - in all rotation cases, assume 50mm
            %%% Distance from the origin to the cortex, as Power used
            
            %%% 1.  As in Power et al 2012: FD
            %%% 2.  As in van Dijk et al 2012: VDt and VDr
            
            [FD VDt VDr]        = deal(NaN(size(rpdiff,1), 1));
            pars={};
            for t = 1:size(rpdiff,1)
                
                phi             = abs(rpdiff(t,4));
                theta            = abs(rpdiff(t,5));
                psi             = abs(rpdiff(t,6));
                
                %%% Framewise Displacement
                FD(t)           = sum([abs(rpdiff(t,1:3)) 50*phi 50*theta 50* psi]);
                %%% Framewise Translation
                VDt(t)          = sqrt(sum([rpdiff(t,1:3).^2]));
                
                %%% Framewise rotation (uses Euler's angle formula from Van Dijk 2012 (p%433)
                VDr(t)          = 50*acos((cos(phi)*cos(theta) + cos(phi)*cos(psi) + cos(theta)*cos(psi) + sin(phi)*sin(psi)*sin(theta) - 1)/2);
            end
            pars{end+1}=FD;
            medianFD=median(FD);
            pars{end+1}=medianFD;
            mv=FD>0.5;
            noOutliers_FD=nnz(mv);
            pars{end+1}=noOutliers_FD;
            pars{end+1}=VDt;
            pars{end+1}=VDr;
        end
        
        function str = array2str_CSV(x)
            %this function converts a numeric 1D array into a comma
            %separated string.
            str='';
            for i=1:length(x);
                if(i>1)
                    str=[str ', '];
                end
                str=[str num2str(x(i))];
            end
        end
        
        function str = vec2str(v,conector)
            str=num2str(v(1));
            if(~exist('conector','var'))
                conector='-';
            end
            for i=2:length(v)
                str=[str conector num2str(v(i))];
            end
        end
        
        function str = cell2str_CSV(x,indexes)
            %this function converts elements specified by indexes
            %into a comma separated string.
            if(~exist('indexes','var'))
                indexes=1:length(x);
            end
            str='';
            for i=1:length(indexes);
                if(i>1)
                    str=[str ', '];
                end
                el=x{indexes(i)};
                if(isnumeric(el))
                    el=num2str(el);
                else
                    el=char(el);
                end
                str=[str el];
            end
        end
        
        function scratching(~)
            %this is a function ment to be modified all the time, serving
            %inter place of command window with the convenience of better
            %editting.
            a=Importdata('/home2/data/process/PATHD/GroupAnalysis/Rest_FA77/QA/Saved/QA_Motion_DVARS.csv');
            a=a.data;
            
        end
        
        function createGrayMatterMasks()
            datadir='/home/JinT/matlab_projects/fMRI/data';
            outdir=fullfile(datadir, 'masks');
            if(~exist(outdir,'file'))
                mkdir(outdir);
            end
            %spm probability map
            %PATHD baseline rest mask
            maskfile=fullfile('/home2/data/process/PATHD/GroupAnalysis/GroupMasks','Time1_mask_50_P1198Excluded.nii');
            spmgrayProb='/usr/local/spm8/apriori/grey.nii';
            spmwhiteProb='/usr/local/spm8/apriori/white.nii';
            spmcsfProb='/usr/local/spm8/apriori/csf.nii';
            volg=spm_read_vols(spm_vol(spmgrayProb));
            volw=spm_read_vols(spm_vol(spmwhiteProb));
            volc=spm_read_vols(spm_vol(spmcsfProb));
            
            volgWFU=spm_read_vols(spm_vol('/home/JinT/matlab_projects/fMRI/data/masks/rwfuGray.nii'));
            volgm=spm_read_vols(spm_vol(maskfile));
            
            volgm=(volgm+volgWFU)>1.9;
            
            V=spm_vol(spmgrayProb);
            volg=volg > 0.4;
            volw=volw > 0.99;
            volc=volc > 0.99;
            
            volgwc=(volg-volw-volc)>0;
            
            V.fname=fullfile(outdir,'spmGray40.nii');
            spm_write_vol(V,volg);
            
            V.fname=fullfile(outdir,'spmWhite99.nii');
            spm_write_vol(V,volw);
            
            V.fname=fullfile(outdir,'spmCSF99.nii');
            spm_write_vol(V,volc);
            
            V.fname=fullfile(outdir,'spmGray40-WMCSF99.nii');
            spm_write_vol(V,volgwc);
            
            V=spm_vol('/home/JinT/matlab_projects/fMRI/data/masks/rwfuGray.nii');
            V.fname=fullfile('/home2/data/process/PATHD/GroupAnalysis/GroupMasks','Time1_mask_50_P1198Excluded_Intersect_wfuGray.nii');
            spm_write_vol(V,volgm);
        end
        
        function col = getDataColumn(columnNames, columnName, data)
            %this function returns a column vector of data whose name
            %(contained in columnNames) matches columnName.
            col=[];
            dim=size(data);
            if(length(columnNames)~=dim(2))
                cprintf('red','%s\n','The number of column names are not consistent with the number of columns of the data');
                return;
            end
            fi=find(ismember(columnNames,columnName));
            if(length(fi)~=1)
                cprintf('red','%s\n','The column name must be a unique element in columnNames');
                return;
            end
            index=fi(1);
            col=data(:,index);
        end
        
        %%%%%%
        %Taihao 7/31/2017
        %The following four functions are from xjView. Placed as member
        %functions for the convenience of calling from other functions.
        %%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% p2t
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function t = p2t(p, df, TF)
            if ~iscell(p)
                if upper(TF)=='T' | upper(TF)=='S'
                    t = spm_invTcdf(1-p,df);
                elseif upper(TF) == 'F'
                    t = spm_invFcdf(1-p,df);
                end
            else
                for ii=1:length(p)
                    t{ii} = p2t(p{ii},df{ii},TF{ii});
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% s2t, s is defined as -log10(p)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function t = s2t(s, df, TF)
            
            if ~iscell(s)
                p = 10^(-s);
                if upper(TF)=='T' | upper(TF)=='S'
                    t = spm_invTcdf(1-p,df);
                elseif upper(TF) == 'F'
                    t = spm_invFcdf(1-p,df);
                end
            else
                for ii=1:length(p)
                    t{ii} = s2t(s{ii},df{ii},TF{ii});
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% t2p
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function p = t2p(t, df, TF)
            if ~iscell(t)
                if upper(TF)=='T' | upper(TF)=='S'
                    p = 1-spm_Tcdf(t,df);
                elseif upper(TF) == 'F'
                    p = 1-spm_Fcdf(t,df);
                end
            else
                for ii=1:length(t)
                    p{ii} = t2p(t{ii},df{ii},TF{ii});
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% t2s, s is defined as -log10(p)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function s = t2s(t, df, TF)
            if ~iscell(t)
                if upper(TF)=='T' | upper(TF)=='S'
                    p = 1-spm_Tcdf(t,df);
                elseif upper(TF) == 'F'
                    p = 1-spm_Fcdf(t,df);
                else
                    p = 0.1;
                end
                s = -log10(p);
            else
                s = {10};
                for ii=1:length(t)
                    s{ii} = t2s(t{ii},df{ii},TF{ii});
                end
            end
        end
        
        function cmap = getDefaultColorMap(N)
            %this function returns N rgb colors. it predefined 12 clors.
            s=255;
            i=1;
            cmap=zeros(N,3);
            %/black
            i=i+1;
            if(i>=N)
                return;
            end
            %/red
            c=[s 0 0];
            cmap(i,:)=c(:);
            i=i+1;
            if(i>=N)
                return;
            end
            c=[0 0 s];%blue
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[s 0 s];%magenta
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[0 s s];%cyan
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[0 s 0];%lime
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[10 0 50];%purple
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[128 0 0];%Maroon
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[0 128 0];%green
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[s 165 0];%orange
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[s 20 147];%deep pink
            cmap(i,:)=c(:);
            
            i=i+1;
            if(i>=N)
                return;
            end
            c=[138 43 226];%blue violet
            cmap(i,:)=c(:);
            
            for j=i+1:N
                c=randi(s,[1, 3]);
                cmap(j,:)=c(:);
            end
        end
        
        function customSubplot()
            % - Define dummy data: 11 time series.
            t       = 0 : 0.1 : 10 ;
            data    = 2 * repmat( sin(t).', 1,11 ) + rand( length(t), 11 ) ;
            nSeries = size( data, 2 ) ;
            
            % - Build figure.
            figure(2) ;  clf ;
            set( gcf, 'Color', [255 255 255]/255, 'Unit', 'Normalized', ...
                'Position', [0.1,0.1,0.8,0.8] ) ;
            
            % - Compute #rows/cols, dimensions, and positions of lower-left corners.
            nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
            rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
            colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
            rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
            
            % - Build subplots axes and plot data.
            for dId = 1 : nSeries
                rowId = ceil( dId / nCol ) ;
                colId = dId - (rowId - 1) * nCol ;
                axes( 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
                plot( t, data(:,dId), 'b' ) ;
                grid on ;
                xlabel( '\theta(t) [rad]' ) ;  ylabel( 'Anomaly [m]' ) ;
                title( sprintf( 'Time series %d', dId )) ;
            end
            
            % - Build title axes and title.
            axes( 'Position', [0, 0.95, 1, 0.05] ) ;
            set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
            text( 0.5, 0, 'My Nice Title', 'FontSize', 14', 'FontWeight', 'Bold', ...
                'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
        end
        
        function moveCurrentAxis(mv)
            %this function moves the position of the current axis by mv
            %(normalized).
            position=get(gca,'Position');
            for i=1:2
                position(i)=position(i)+mv(i);
            end
            set(gca,'Position', position);
        end
        
        function markers = getLineMarkers(N)
            %            b     blue          .     point              -     solid
            %            g     green         o     circle             :     dotted
            %            r     red           x     x-mark             -.    dashdot
            %            c     cyan          +     plus               --    dashed
            %            m     magenta       *     star             (none)  no line
            %            y     yellow        s     square
            %            k     black         d     diamond
            %            w     white         v     triangle (down)
            %                                ^     triangle (up)
            %                                <     triangle (left)
            %                                >     triangle (right)
            %                                p     pentagram
            %                                h     hexagram
            markers0 = {'o' 's' 'd' 'v' '*' 'x' '.' 'p' 'h' '^' '<' '>' '+'};
            if(~exist('N','var'))
                markers=markers0;
                return;
            end
            markers=cell(1,N);
            len=length(markers0);
            m=floor(N/len);
            r=floor(N-m*len+0.5);
            markers=repmat(markers0,1,m);
            markers=horzcat(markers, markers0(1:r));
        end
        
        function range = getRange(m)
            %returns a structure containing min, max and max-min of a
            %matrix
            dim=size(m);
            n=min(m);
            x=max(m);
            for i=2:length(dim)
                n=min(n);
                x=max(x);
            end
            range.min=n;
            range.max=x;
            range.range=x-n;
        end
        
        function vec = mat2vec(m)
            vec=reshape(m,[1,numel(m)]);
        end
        
        function summary = getSummary(m)
            %returns summary structure of the matrix m
            num=numel(m);
            m=CommonMethods.mat2vec(m);
            summary.range=CommonMethods.getRange(m);
            %standard deviation
            summary.sd=std(m);
            %standard erro of mean
            summary.sem=summary.sd/sqrt(num);
            summary.mean=mean(m);
            summary.median=median(m);
            summary.mad=mad(m,1);
        end
        
        function cdHome()
            cd '/home/JinT/matlab_projects/fMRI/';
        end
        
        function c = insertIntoCellArray(c,elm, index)
            %this function inserts an element into a 1D cellarray pointed
            %to by index
            dim=size(c);
            if(length(dim)>2||(dim(1)>1&&dim(2)>1))
                cprintf('r','%s\n','the cell array for the function insertIntoCellArray must be 1D');
                return;
            end
            if(dim(1)>1)
                c = (CommonMethods.insertIntoCellArray_1Row(c',elm,index))'
                return;
            end
            c = CommonMethods.insertIntoCellArray_1Row(c,elm,index);
        end
        
        function c = insertIntoCellArray_1Row(c,elm, index)
            %is function inserts an element into one row cellarray
            c=[c(1:index-1) elm c(index:end)];
        end
        
        function exportTable_CSV(filePath,columnNames,rowNames,data,del)
            %this function export a matrax into CSV file.
            %columnNames should be a cell array whose length should be
            %equalt to or one more than the number of columns in data. It
            %will be ignored, it does not meet the requirement.
            %rowName should be a cellarray whose length equals the number
            %of rows in data matrix. rowName will be ignored if the length
            %of this cell array does not meet this requirement.
            dim=size(data);
            
            if(length(dim)~=2)
                cprintf('r','%s\n','data for exportTable_CSV must be a 2D array');
                return;
            end
            
            if(~exist('del','var'))
                del=', ';
            end
            
            if(length(rowNames)~=dim(1))
                rowNames = {};
            end
            
            if(length(columnNames)==dim(2))
                if(~isempty(rowNames))
                    columnNames=CommonMethods.insertIntoCellArray(columnNames,' ',1);
                end
            else
                if((length(columnNames)-dim(2))~=1)
                    columnNames={};
                end
            end
            
            ext=CommonMethods.getFileExtension(filePath);
            if(~strcmpi(ext,'csv')&&~strcmpi(ext,'tsv'))
                filePath=[filePath '.csv'];
            end
            fid = fopen(filePath,'wt');
            
            if(~isempty(columnNames))
                fprintf(fid,'%s\n',CommonMethods.cell2str(columnNames,del));
            end
            
            for i=1:dim(1)
                if(iscell(data))
                    line=CommonMethods.cell2line(data(i,:),del);
                else
                    line=CommonMethods.array2str_CSV(data(i,:));
                end
                if(~isempty(rowNames))
                    line=[rowNames{i} ', ' line];
                end
                fprintf(fid,'%s\n', line);
            end
            fclose(fid);
        end
        
        function c = getColumnCells(m)
            %this function returns a cell array of each column in the
            %matrix m. m is assumed to be a 2D mumeric array.
            dim=size(m);
            c=mat2cell(m,dim(1),ones(1,dim(2)));
        end
        
        function m = getColumns(m0, indexes)
            %this function returns a matrix of columns of m specified by
            %indexes
            dim=size(m0);
            m=zeros(dim(1),length(indexes));
            for i=1:length(indexes)
                m(:,i)=m0(:,indexes{i});
            end
        end
        
        function [rows, cols]=getSubplotArrangement(numPlots)
            sq=sqrt(numPlots);
            rows=floor(sq);
            cols=ceil(sq);
            if(cols*rows<numPlots)
                rows=rows+1;
            end
        end
        
        function displayDataColumnComparison(fh,data,columnNames,fontSize, fontName)
            if(~exist('fontSize','var'))
                fontSize=14;
            end
            if(~exist('fontName','var'))
                fontName='Arial';
            end
            if(~exist('columnNames','var'))
                columnNames=arrayfun(@(x) ['Column' num2str(x)], 1:size(data,2), 'UniformOutput', false);
            end
            
            
            [Rho PVal]=corr(data,'Type','Pearson');
            
            %             dim=size(data);
            %             cols=dim(2);
            %             rows=cols;
            %             num=1;
            len=size(data,2);
            [rows, cols]=CommonMethods.getSubplotArrangement(len*(len-1)/2);
            num=1;
            
            for r=1:len
                for c=1:len
                    
                    if(c>=r)
                        continue;
                    end
                    subplot(rows,cols,num);
                    
                    x=data(:,c);
                    y=data(:,r);
                    line(x,y,'Color',[0 0 1],'LineStyle','none','Marker', 'o', 'MarkerSize', 3);
                    b=regress(y,x);
                    t1=['b= ' num2str(b,2)];
                    t2=['r= ' num2str(Rho(r,c),3)];
                    t3=['p= ' num2str(PVal(r,c),2)];
                    
                    xlabel(gca, columnNames{c},'FontName', fontName, 'FontSize', fontSize);
                    %                     if(r==rows)
                    %                         xlabel(gca, columnNames{c},'FontName', fontName, 'FontSize', fontSize);
                    %                     else
                    %                         set(gca,'xticklabel',[]);
                    %                     end
                    
                    ylabel(gca, columnNames{r},'FontName', fontName, 'FontSize', fontSize);
                    %                     if(c==1)
                    %                         ylabel(gca, columnNames{r},'FontName', fontName, 'FontSize', fontSize);
                    %                     else
                    %                         set(gca,'yticklabel',[]);
                    %                     end
                    lsline;
                    xRange=xlim;
                    yRange=ylim;
                    position=get(gca,'Position');
                    
                    %                   text(xRange(1),yRange(2),{t1,t2,t3});
                    title([t1 ' ' t2 ' ' t3], 'fontsize', 14);
                    num=num+1;
                    
                    set(get(gca,'xlabel'),'fontsize', 14);
                    set(get(gca,'ylabel'),'fontsize', 14);
                end
            end
        end
        
        function cleanDirs(InfoOrganizer)
            subjectNodes=InfoOrganizer.AnalysisInfo.SubjectNodes;
            if(true) return;end;%Comment this line out to clearn dirs
            for i=1:length(subjectNodes)
                %        for i=1:1
                if ~fMRIAnalysisInfoOrganizer.completedScans(subjectNodes{i},{'a' 'b'})
                    continue;
                end
                scanNodes=subjectNodes{i}.ScanNodes;
                ParadigmName=InfoOrganizer.AnalysisInfo.ParadigmName;
                for j=1:length(scanNodes)
                    rundirs=scanNodes{j}.RunDirs;
                    for k=1:length(rundirs)
                        CommonMethods.cleanFunctionalDir(rundirs(k));
                    end
                    anatDir=strrep(rundirs{1},ParadigmName,'T1');
                    CommonMethods.cleanAnatDir(anatDir);
                end
            end
        end
        
        function cleanFunctionalDir(dir0)
            dir=dir0{1};
            if(exist(dir,'dir'))
                cprintf('blue', '\nDir : %s \n', dir);
                c = strsplit(dir, filesep);
                %to makesure not deleting the original data
                if(strcmp(c{length(c)},'analysis'))
                    cd (dir);
                    cprintf('blue', 'Dir : %s \n', dir);
                    system(['rm *.*']);
                    system(['cp ' fullfile(strrep(dir, 'analysis', 'nifti'), 'f.nii') ' ' dir]);
                    system(['cp ' fullfile(strrep(dir, 'analysis', 'nifti'), 'f.bxh') ' ' dir]);
                end
            end
        end
        
        function cleanAnatDir(dir0)
            if(exist(dir0,'dir'))
                cprintf('blue', '\nDir : %s \n', dir0);
                c = strsplit(dir0, filesep);
                %to makesure not deleting the original data
                if(strcmp(c{length(c)},'analysis'))
                    cd (dir0);
                    cprintf('blue', 'Dir : %s \n', dir0);
                    
                    allFiles=dir(dir0);
                    names={allFiles.name};
                    for i=1:length(names)
                        name=names{i};
                        if~isempty(find(ismember({'f.bxh', 'f.nii', 'mf.nii', 'c1f.nii', 'c2f.nii', 'c3f.nii'},name)))
                            continue;
                        end
                        if(length(name)<5)
                            continue;
                        end
                        system(['rm ' name]);
                    end
                    
                end
            end
        end
        
        function compareFiles()
            dir1='/home2/data/process/PATHD/P1110/P1110a/EMOID/1/analysis/savedMasks';
            dir2='/home2/data/process/PATHD/P1110/P1110a/EMOID/1/analysis';
            name1='swmaskU1.nii';
            name2='swmaskU1.nii';
            file1=fullfile(dir1,name1);
            file2=fullfile(dir2,name2);
            %           [~, ~] = system(['gunzip ' file2 '.gz']);
            visdiff(file1,file2);
        end
        function data = fillColumn(data,column)
            %this function fills empty cells of the column of data using
            %the string from the previous non-empty cell.
            dim=size(data);
            rows=dim(1);
            st0=data{1,column};
            for r=1:rows
                st=data{r,column};
                if(strcmp(st,''))
                    data{r,column}=st0;
                else
                    st0=st;
                end
            end
        end
        function row = getRow(data, columns, contents)
            %this function returns the first row on which the cells in
            %columns match the centents
            dim=size(data);
            rows=dim(1);
            for r=1:rows
                match=true;
                for i=1:length(columns)
                    c=columns{i};
                    if(~strcmp(contents{i},data{r,c}))
                        match=false;
                        break;
                    end
                    if(match)
                        row=r;
                        return;
                    end
                end
            end
        end
        
        function displayIndexedTimeSeries(data, signalNames, uniformYScales)
            if(~exist('uniformYScales', 'var'))
                uniformYScales=false;
            end
            if(~exist('signalNames','var'))
                signalNames=arrayfun(@(x)num2str(x),1:size(data,2),'UniformOutput', false);
            end
            
            dim=size(data);
            cols=dim(2);
            N=cols*(cols-1)/2;
            num=1;
            indexes=(1:dim(1));
            for c=1:cols
                for c2=c+1:cols
                    if(c==c2)
                        continue;
                    end
                    subplot(N,1,num);
                    y1=data(:,c);
                    y2=data(:,c2);
                    
                    [hAx,line1,line2]=plotyy(indexes,y1,indexes,y2);
                    ylabel(hAx(1),signalNames{c});
                    ylabel(hAx(2),signalNames{c2});
                    set(line1,'linewidth',2);
                    set(line2,'linewidth',2);
                    set(line1, 'displayname', signalNames{c});
                    set(line2, 'displayname', signalNames{c2});
                    if(uniformYScales)
                        lim1=[min(y1) max(y1)];
                        lim2=[min(y2) max(y2)];
                        range1=lim1(2)-lim1(1);
                        range2=lim2(2)-lim2(1);
                        if(range1>range2)
                            CommonMethods.copyYscale(hAx(1),hAx(2));
                        else
                            CommonMethods.copyYscale(hAx(2),hAx(1));
                        end
                    end
                    if(num==N)
                        xlabel('Time (index)');
                    end
                    num=num+1;
                    set(get(gca,'xlabel'),'fontsize', 14);
                    set(get(gca,'ylabel'),'fontsize', 14);
                    
                    set(get(hAx(2),'xlabel'),'fontsize', 14);
                    set(get(hAx(2),'ylabel'),'fontsize', 14);
                end
            end
        end
        
        function copyYscale(hAx1, hAx2)
            %this function set y scale of hAx2 the same range and ticks as
            %hAx1, but keep the original center
            ylim1=get(hAx1,'ylim');
            ylim2=get(hAx2,'ylim');
            m1=mean(ylim1);
            m2=mean(ylim2);
            ylim2=ylim1-m1+m2;
            
            yTick1=get(hAx1,'ytick');
            yTick2=yTick1-m1+m2;
            set(hAx2,'ylim',ylim2,'yTick',yTick2);
        end
        
        function str=cell2str_csv(strs,del)
            if ~exist('del','var')
                del=',';
            end
            
            len=0;
            for i=1:length(strs)
                len=len+length(strs{i})+2;
            end
            len=len-2;
            
            p=0;
            for i=1:length(strs)
                st=strs{i};
                for j=1:length(st);
                    p=p+1;
                    str(1,p)=st(j);
                end
                if(i==length(strs))
                    break;
                end
                p=p+1;
                str(1,p)=del;
                p=p+1;
                str(1,p)=' ';
            end
        end
        
        function labels = getLabels(names0,labels0,names)
            %this function returns a cell array of strings (labels) chosen
            %from labels0 at positions at which the contents of names0
            %match those in names
            indexes = CommonMethods.findIndexes(names0,names);
            labels=cellfun(@(x) labels0(x), indexes, 'UniformOutput', false);
        end
        function index = findIndex(strs, str)
            fi=find(ismember(strs,str));
            if(isempty(fi))
                index=0;
                return;
            end
            index=fi(1);
        end
        function indexes = findIndexes(strs0, strs)
            %this function returns a cell arrays of indexes at whitch the
            %contents of strs0 are contained in strs
            indexes={};
            for i=1:length(strs)
                indexes = horzcat(indexes, find(ismember(strs0,strs(i))));
            end
        end
        
        function markVerticalLines_Figure(fh,X,cm)
            ahs=findall(fh,'type','axes');
            for i=1:length(ahs)
                CommonMethods.markVerticalLines_SingleAxis(ahs(i),X,cm);
            end
        end
        
        function markVerticalLines_SingleAxis(ah,X,cm)
            %this function mark each x value in X using a verticle line
            %(colored by cm) across the entire y range.
            if(~exist('cm','var'))
                cm=[1 0 0];
            end
            
            axes(ah);
            yr=ylim();
            
            for i=1:length(X)
                x=X(i);
                line([x x],yr,'color',cm);
            end
        end
        function refreshAxes(fh)
            ahs=findall(fh,'type','axes');
            for i=1:length(ahs)
                axes(ahs(i));
            end
        end
        function [vol,file]=loadVols()
            %this function use spm_select to select images and return the
            %volumes and the file names.
            file0=spm_select([0 Inf],'image','Select an image file to load the volume');
            if(isempty(file0))
                vol={};
                file={};
                return;
            end
            dim=size(file0);
            len=dim(1);
            vol=cell(1,len);
            file=cell(1,len);
            for i=1:len
                f=file0(i,:);
                file{i}=f;
                vol{i}={spm_read_vols(spm_vol(f))};
            end
        end
        function o = loadFile()
            %this function use spm_select to select images and return the
            %volumes and the file names.
            file=spm_select(1,'mat','Select the file to be loaded in the work space');
            if(isempty(file))
                return;
            end
            dim=size(file);
            len=dim(1);
            for i=1:len
                f=file(i,:);
                o=load(f);
            end
        end
        
        function copyNifitiFiles(scanID,home)
            %this is to copy every files in nifti directory of home to
            %analysis directory for scanID.
            mysql('open','localhost','root','denali','MRI');
            home=[filesep home];
            
            if nargin ~= 2
                %                images                                                                  = mysql('select * from analyzeData;');
                return;
            else
                images                                                                  = mysql(['select * from analyzeData where scanID=''' scanID ''';']);
            end
            
            for i = 1:size(images,1)
                
                temp                                                                    = images(i);
                scaninfo                                                                = mysql(['select * from scans where scanID=''' deblank(temp.scanID) ''';']);
                
                if strcmp(temp.studyName,'Prodrome')
                    scanfolder                                                          = [strrep(scaninfo.scandate,'-',''), '_', deblank(temp.scanID)];
                else
                    scanfolder                                                          = scanID;
                end
                
                niftiPath                                                             = fullfile(home, 'data', 'process', deblank(temp.studyName), deblank(scaninfo(1).subjectID), scanfolder, deblank(temp.paradigmName), num2str(temp.instance), 'nifti');
                path2=strrep(niftiPath,'nifti','analysis');
                if(exist(path2,'file'))
                    cmd=['rm ' path2 filesep '*'];
                    cprintf('blue','%s\n',cmd);
                    system(cmd);
                else
                    mkdir(path2);
                end
                cd (niftiPath);
                cmd=['cp ' '* ' path2];
                cprintf('blue','%s\n',cmd);
                system (cmd);
                %         system(['cp ' fullfile(strrep(rundir{r}, 'analysis', 'nifti'), 'f.nii') ' ' rundir{r}]);
            end
            mysql('close');
        end
        
        function dist = pdist2(p1,p2)
            len1=size(p1,1);
            len2 = size(p2,1);
            dist = zeros(len1,len2);
            for i=1:len1
                for j=1:len2
                    dist(i,j)=CommonMethods.dist(p1(i,:),p2(j,:));
                end
            end
        end
        
        function str = Ind2str(sub, sep)
            len=length(sub);
            str = num2str(sub(1));
            for i=2:len
                str=[str sep num2str(sub(i))];
            end
        end
        
        function dist = dist(p,p0)
            dp=p-p0;
            dp2=dp.*dp;
            dist=sqrt(sum(dp2));
        end
        
        function p=getNormalizedPosition_Axes(axis, p0)
            %            axes(axis);
            xr=xlim;
            yr=ylim;
            x=CommonMethods.LinearInterpolation(0,xr(1),1,xr(2),p0(1));
            y=CommonMethods.LinearInterpolation(0,yr(1),1,yr(2),p0(2));
            p=[x y];
        end
        
        function openShowsrs3_new()
            anat    = readmr('/home2/data/process/PATHD/P0001/P0001a/EMOID/1/analysis/rmT1.nii');
            %            HC      = readmr('/home2/data/process/PATHD/P0001/P0001a/EMOID/1/analysis/arf.nii', 'NOPROGRESSBAR');
            HC      = readmr('/home2/data/process/PATHD/P0001/P0001a/Rest_FA77/1/analysis/arf.nii', 'NOPROGRESSBAR');
            rangeHC=[min(min(min(HC.data))) max(max(max(HC.data)))];
            Tmap_HC = readmr(['/home2/data/process/PATHD/P0001/P0001a/EMOID/stats_arf_aCompCor/spmT_0005.img'], 'NOPROGRESSBAR');
            con_HC = readmr(['/home2/data/process/PATHD/P0001/P0001a/EMOID/stats_arf_aCompCor/con_0005.img'], 'NOPROGRESSBAR');
            rangeTmap_HC=[min(min(min(Tmap_HC.data))) max(max(max(Tmap_HC.data)))];
            %             neg_Tmap_HC      = Tmap_HC;
            %             neg_Tmap_HC.data = neg_Tmap_HC.data*-1;neg_Tmap_HC      = Tmap_HC;
            %             neg_Tmap_HC.data = neg_Tmap_HC.data*-1;
            
            %             h1 = showsrs3_new(...
            %                 anat,        struct('cmapLim', [],         'cmap', gray(256)),...
            %                 Tmap_HC,     struct('cmapLim', rangeTmap_HC, 'cmap', redgrad(256),  'transLevel',0.5, 'showTimePlot',0),...
            %                 HC,          struct('cmapLim', rangeHC, 'cmap', redpos(256),   'transLevel',0,   'showTimePlot',1));
            
            showsrs3(con_HC,[]);
            % 'cmapLim'      controls the min and max values for the overlays.
            % 'cmap'         is the colormap used for the overlay.
            % 'transLevel'   is the transparency level for the overlay
            % 'showTimePlot' 0=don't show timeseries for this overlay, 1=show overlay.
        end
        
        function [desMats, condNames] = exactDesignMatrices(matfile)
            load(matfile);
            noRuns=length(SPM.Sess);
            desMats=cell(1,noRuns);
            condNames=cell(1,noRuns);
            noVols=length(SPM.Sess(1).row);
            %           figure();
            for r = 1:noRuns
                noConds=length(SPM.Sess(r).U);
                model_desMat                                                        = SPM.xX.X((r-1)*noVols+1:r*noVols, find(cellfun(@(x) ~isempty(strfind(x, ['Sn(' num2str(r) ')'])) & ~isempty(strfind(x, 'bf(1)')), {SPM.Vbeta.descrip})));
                desMats{r}=model_desMat;
                condNames{r}=arrayfun(@(x) char(SPM.Sess(r).U(x).name{1}),1:noConds,'UniformOutput', false);
                %                subplot(2,3,r);
                %                plot(model_desMat);
            end
            condMat=CommonMethods.getCondition(desMats{3},condNames{3},condNames{3}{3});
        end
        
        function condMat = getCondition(desMat,condNames,condName)
            indexes=find(ismember(condNames,condName));
            num=length(indexes);
            condMat = desMat(:,indexes);
        end
        
        function displayTimeSeries(v,disp)
            persistent displayer;
            %           displayer=disp;
            if(exist('disp','var'))
                displayer=disp;
            end
            
            if(isa(displayer,'TimeSeriesDisplayer'))
                displayer.display(v);
            end
        end
        
        function display4D()
            %            showsrs3_new();
            InfoOrganizer=PATHD_getInfoOrganizer('EMOID');
            InfoOrganizerR=PATHD_getInfoOrganizer('Rest_FA77');
            subjectID='P1110';
            scanID=[subjectID 'a'];
            InfoOrganizer.selectSubjectNodes({subjectID});
            InfoOrganizerR.selectSubjectNodes({subjectID});
            rundirs=InfoOrganizer.getRunDirs(InfoOrganizer.AnalysisInfo.SubjectNodes{1});
            rundirsR=InfoOrganizerR.getRunDirs(InfoOrganizerR.AnalysisInfo.SubjectNodes{1});
            files=cellfun(@(x) fullfile(x, 'arf.nii'), rundirs,'UniformOutput',false);
            files=horzcat(files,cellfun(@(x) fullfile(x, 'arf.nii'), rundirsR,'UniformOutput',false));
            
            taskfolder='stats_arf_aCompCor';
            datadir=InfoOrganizer.AnalysisInfo.dataDir;
            statdir=fullfile(datadir,subjectID,scanID,InfoOrganizer.AnalysisInfo.ParadigmName,taskfolder);
            matfile=fullfile(statdir,'SPM.mat');
            [desMats, condNames]=CommonMethods.exactDesignMatrices(matfile);
            cond='emotion';
            yy={};
            
            for i=1:size(desMats,2)
                names=cellfun(@(x) x(1:(length(x)-2)), condNames{i}, 'UniformOutput',false);
                yy{end+1}=CommonMethods.getCondition(desMats{i},names,cond);
            end
            for i=1:length(rundirsR)
                yy{end+1}=[];
            end
            
            statsdir='/home2/data/process/PATHD/P1110/P1110a/EMOID/stats_arf_aCompCor';
            cname=[cond 'T1'];
            img=CommonMethods.getConImg(statsdir, cname);
            
            displayer = TimeSeriesDisplayer(files,yy,statsdir,cond);
            CommonMethods.displayTimeSeries([23 23 23],displayer);
            showsrs3_new(readmr(img),[]);
        end
        
        function img = getConImg(statsdir,cname)
            load(fullfile(statsdir,'SPM.mat'));
            cnames=arrayfun(@(x) SPM.xCon(x).name, 1:length(SPM.xCon),'UniformOutput', false);
            fi=find(ismember(cnames,cname));
            if(length(fi)~=1)
                img='';
                cprintf('red','%s\n', ['the contrast name ' cname ' is not found in SPM.xCon']);
                return;
            end
            img=fullfile(statsdir, sprintf('con_%04d.img', fi(1)));
        end
        
        function setEqualXScales(Axes)
            len=length(Axes);
            range=get(Axes{1},'xlim');
            for a=2:len
                xs=get(Axes{a},'xlim');
                range(1)=min(range(1),xs(1));
                range(2)=max(range(2),xs(2));
            end
            for a=1:length(Axes)
                set(Axes{a},'xlim',range);
            end
        end
        
        function setEqualYScales(Axes)
            len=length(Axes);
            range=get(Axes{1},'ylim');
            for a=2:len
                xs=get(Axes{a},'ylim');
                range(1)=min(range(1),xs(1));
                range(2)=max(range(2),xs(2));
            end
            for a=1:length(Axes)
                set(Axes{a},'ylim',range);
            end
        end
        
        
        function setEqualYScales1(Axes)
            len=length(Axes);
            axes(Axes{1});
            range=ylim;
            for a=2:len
                axes(Axes{2});
                xs=ylim;
                range(1)=min(range(1),xs(1));
                range(2)=max(range(2),xs(2));
            end
            for a=1:length(Axes)
                set(Axes{a},'ylim',range);
            end
        end
        
        function dispRange = getDispRange(y,range)
            %return the min and max, max-min=range to display y in the
            %middle.
            w=range/2;
            m=(min(y)+max(y));
            dispRange=[m-w m+w];
        end
        
        function rv = getResidualVolume(imgfile,spmfile,runNo)
            r=runNo;
            load(spmfile);
            vol=spm_read_vols(spm_vol(imgfile));
            model_desMat = SPM.xX.X((r-1)*noVols+1:r*noVols, find(cellfun(@(x) ~isempty(strfind(x, ['Sn(' num2str(r) ')'])) & ~isempty(strfind(x, 'bf(1)')), {SPM.Vbeta.descrip})));
            desMat=SPM.xX.X
        end
        
        function condNames = getCondNames_FristLevelModel(SPM,runNo)
            load('/home2/data/process/PATHD/P1110/P1110a/EMOID/stats_arf_aCompCor/SPM.mat');
        end
        
        function run = sleepUntil(stateFun, args, increament)
            pause('on');
            cprintf('blue','pause until ''stateFun'' become true. (%s)\n', datestr(clock));
            while(~stateFun(args))
                pause(increament);
            end
            run=true;
            cprintf('blue','''stateFun'' become true. (%s)\n', datestr(clock));
            return;
        end
        
        function allExist = allFilesExist(files)
            allExist=~any(cellfun(@(x) ~exist(x,'file'), files));
        end
        
        function existingFIles=selectExistingFiles(files)
            existingFIles=files(find(cellfun(@(x) exist(x,'file'),files)));
        end
        
        function existingFIles=selectNonExistingFiles(files)
            existingFIles=files(find(cellfun(@(x) ~exist(x,'file'),files)));
        end
        
        function std = mad2std(mad)
            std=1.4826*mad;
        end
        
        function past=PastTheTime(time)
            %this function test whether current time is past the time
            %specified by "time". time=[yyyy mm dd hour min sec ....].
            %"time" can have any number of elements.
            cur_time=clock();
            len=length(time);
            past=false;
            for i=1:len
                if(cur_time(i) == time(i))
                    continue;%check the next element
                elseif(cur_time(i) > time(i))%passed the time
                    past=true;
                    return;
                else %not yet
                    return;
                end
            end
        end
        
        function pval = getGaussianP_MAD (vec0)
            %this function compute the p value of each element of vec (a
            %vector) in a gaussian distribution whose mean and standard
            %deviation is median and 1.4826*mad.
            
            mad0=mad(vec0,1);
            if(mad0==0)
                mad0=1;
            end
            std0=CommonMethods.mad2std(mad0);
            median0=median(vec0);
            vec=(vec0-median0)/std0;
            pval=normcdf(vec);
            indexes=find(pval>0.5);
            pval(indexes)=1-pval(indexes);
        end
        
        function displayFiles(files)
            eFiles=files(find(cellfun(@(x) exist(x,'file'),files)));
            nFiles=files(find(cellfun(@(x) ~exist(x,'file'),files)));
            nNum=length(nFiles);
            eNum=length(eFiles);
            cprintf('blue','List of existing files (%s) : (name and modification date)\n',num2str(eNum));
            for i=1:length(eFiles)
                file=eFiles{i};
                list=dir(file);
                cprintf('blue','%s %s:    %s\n', num2str(i), file, list.date);
            end
            cprintf('blue','\n\n List of non-exesiting files (%s):\n',num2str(nNum));
            for i=1:length(nFiles)
                file=nFiles{i};
                cprintf('blue','%s %s (file does not exist)\n', num2str(i), file);
            end
        end
        
        function pval = mat2pval_col(data)
            %this function compute pValue for each element in the matrix
            %data and stores at the corresponding element in the matrix
            %pval.
            %Each element (i,j) of pval is the pvalue of the element of data
            %assuming the j-th column is from a normal distribution with
            %mean and std is median and 1.4826 of the MAD of the column,
            %respectively.
            len=size(data,2);
            pval=-1*ones(size(data));
            cprintf('blue','mat2pval_col: %s)\n',datestr(clock));
            for i=1:len
                vec=data(:,i);
                p=CommonMethods.getGaussianP_MAD(vec);
                pval(:,i)=p;
                if(mod(i,1000)==0)
                    cprintf('blue','mat2pval_col: %s of %s\n',num2str(i),num2str(len));
                    cprintf('blue','mat2pval_col: %s)\n',datestr(clock));
                end
            end
        end
        
        function pval = mat2pval(data,dim)
            %this function compute pValue for each element in the matrix
            %data and stores at the corresponding element in the matrix
            %pval.
            %if dim==1, pval is computed from each row, other wise pval is
            %computed from each column.
            if(~exist('dim','var')) dim=2; end;
            if(dim==1)
                data=data';
                pval=CommonMethods.mat2pval_col(data);
                pval=pval';
            else
                pval=CommonMethods.mat2pval_col(data);
            end
        end
        
        function strs = removeStrings(strs0, strsR)
            %strs0 and strsR are cellarray of str.
            %this function remove cells in strs0 that are memmbers of
            %strsR and return the remain cells as a new cell array.
            strs={};
            for i=1:length(strs0)
                str=strs0{i};
                if(any(ismember(strsR,str)))
                    continue;
                end
                strs{end+1}=str;
            end
        end
        
        function volIndex = mni2ind_spm(st,coor,volIndex)
            %this function converts mmi coordinates into the matrix indexes
            if(~exist('volIndex','var'))
                volIndex=1;
            end
            V=st.vols{volIndex};
            mat=V.mat;
            invm=inv(mat);
            if(size(coor,1)==1) coor=coor';end;
            coor=cat(1,coor,[1]);
            volIndex=round(inv(mat)*coor);
            volIndex=volIndex(1:3)';
            volIndex(volIndex<1)=1;
            volIndex(volIndex>V.dim)=V.dim(volIndex>V.dim);
        end
        
        function coor = ind2mni_spm(st, ind, volIndex)
            %this function converts maxtrix indexes into mmi coor
            if(~exist('volIndex','var'))
                volIndex=1;
            end
            V=st.vols{volIndex};
            mat=V.mat;
            ind=[ind(1) ind(2) ind(3) 1];
            coor = mat*ind';
            coor = coor(1:3);
        end
        
        function coor = ind2mni_spmMat(mat, ind)
            %this function converts maxtrix indexes into mmi coor
            ind=[ind(1) ind(2) ind(3) 1];
            coor = mat*ind';
            coor = coor(1:3);
        end
        
        function coor=ind2coor_mat(ind, mat)
            %ind is a mby3 matrix. each row is three voxel indexes
            %voxel indexes to 3d coordinates
            if ~size(ind,2)==3
                ind=ind';
            end
            ind=[ind ones(size(ind,1),1)];
            coor=(mat*ind')';
            coor=coor(:,1:3);
        end
        
        function ind=coor2ind_mat(coor, mat)
            %this function converts the world coordinates into the voxel
            %indexes
            %coor is a mby3 matrix. each row is three voxel indexes
            
            if isempty(coor)
                ind = [];
                return;
            end
            if(size(coor,2)~=3)
                coor=coor';
            end
            ind = [coor(:,1) coor(:,2) coor(:,3) ones(size(coor,1),1)]*(inv(mat))';
            ind = ind(:,1:3);
            ind = round(ind);
        end
        
        function [st, R4o] = resetOrigin_spm(st0,coor)
            if(~exist('coor','var')) coor=st0.centre';end;
            if(size(coor,1)==1) coor=coor';end;
            st=st0;
            for i=1:length(st0.vols)
                V=st0.vols{i};
                if isempty(V)
                    continue;
                end
                mat=V.mat;
                inds=CommonMethods.coor2ind_mat(coor, mat);
                coor=CommonMethods.ind2coor_mat(inds, mat);
                R4=eye(4);
                R4(:,4)=R4(:,4)-[coor 0]'
                mat=R4*mat;
%                mat(1:3,4)=mat(1:3,4)-coor';
                st.vols{i}.mat=mat;
                if(i==1)
                    R4o=R4;
                end
            end
            st.centre=[0 0 0]';
        end
        
        
        function [st, R4o] = toYAxis_spm(st0,coor)
            %this function replaces the matrix of each volume in st
            %such that the coordinate of the point coor will be on Y axis
            %when direction==-1, the point is aligned to the negative axis
                if(~exist('coor','var')) coor=st0.centre';end;
                
            for i=1:length(st0.vols)
                
                V=st0.vols{i};
                if(isempty(V))
                    continue;
                end
                mat=V.mat;
                inds=CommonMethods.coor2ind_mat(coor, mat);
                coor=CommonMethods.ind2coor_mat(inds, mat);
                x=coor(1);
                y=coor(2);
                z=coor(3);
                
                v1=[0 -1 0]';
                v2=coor';
                R=CommonMethods.computeRotationMatrix_align2vec(v1, v2);
                R4=eye(4);
                R4(1:3,1:3)=R;
                mat=R4*mat;
                if(i==1)
                    R4o=R4;
                end
                st0.vols{i}.mat=mat;
            end
            st=st0;
        end
        function [st, R4o] = toYZPlane_spm(st0,coor)
            %this function replaces the matrix of each volume in st
            %such that the coordinate of the point coor will be on Y axis
            %when direction==-1, the point is aligned to the negative axis
                if(~exist('coor','var')) coor=st0.centre';end;
                
            for i=1:length(st0.vols)
                
                V=st0.vols{i};
                if(isempty(V))
                    continue;
                end
                mat=V.mat;
                inds=CommonMethods.coor2ind_mat(coor, mat);
                coor=CommonMethods.ind2coor_mat(inds, mat);
                x=coor(1);
                y=coor(2);
                z=coor(3);
                
                v1=[0 0 1]';
                v2=coor';
                R=CommonMethods.computeRotationMatrix_align2vec(v1, v2);
                R4=eye(4);
                R4(1:3,1:3)=R;
                mat=R4*mat;
                
                if(i==1)
                    R4o=R4;
                end
               st0.vols{i}.mat=mat;
            end
            st=st0;
        end
        
        function R=computeRotationMatrix_align2vec(v1, v2)
            %this function computes a rotation matrix that makes v2
            %parallel to v1. 
            %R: the rotation matrix;
            %v1, v2 are to 3D column vectors;
            %from https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
            vu = @(v) v/norm(v);
            GG = @(A,B) [ dot(A,B) -norm(cross(A,B)) 0; norm(cross(A,B)) dot(A,B)  0; 0              0           1];
            FFi = @(A,B) [ A (B-dot(A,B)*A)/norm(B-dot(A,B)*A) cross(B,A) ];
            UU = @(Fi,G) Fi*G*inv(Fi);
            a=vu(v2);
            b=vu(v1);
            R=UU(FFi(a,b), GG(a,b));
        end
        
        function st = toYZPlane_spm0(st0, coor)
            %this function replaces the matrix of the first volume in st
            %such that the coordinate of the point coor will be on YZ plane
            if(~exist('coor','var')) coor=st0.centre;end;
            x=coor(1);
            y=coor(2);
            z=coor(3);
            
            for i=1:length(st0.vols)
                
                V=st0.vols{i};
                if isempty(V)
                    continue;
                end
                mat=V.mat;
                %rotation around y axis
                z1=sqrt(x*x+z*z);
                st=-x/z1;
                ct=z/z1;
                Ry=[ct 0 st; 0 1 0; -st 0 ct];
                
                
                R0=mat(1:3,1:3);
                T0=mat(1:3,4);
                
                T1=Ry*T0;
                R1=Ry*R0;
                
                mat(1:3,1:3)=R1;
                mat(1:3,4)=T1;
                
                V.mat=mat;
                
                st0.vols{i}.mat=mat;
            end
        end
        
        function st=fromSPMExt(action,args)
            st=false;
            switch action
                case 'Reposition'
                    if(length(args)==2)
                        if(strcmp(args{2},'fromSPMExt'))
                            st=true;
                        end
                    end
            end
        end
        
        function spm_setZoom(index)
            %for spm*: index should be between 1 and 8, for the zooming
            %option 'BBox(Y> ...' , 'BBox (nonzero)', '10mm', '20mm',
            %'40mm', '80mm', '160mm', and 'full volume', respectively.
            [zl rl] = spm_orthviews('ZoomMenu');
            cz=index;
            spm_orthviews('Zoom',zl(cz),rl(cz));
        end
        
        function spm_setView(st,coor)
            %this function set the view center of spm to coor
            spm_orthviews('Reposition',coor,'fromSPMExt');
            [zl rl] = spm_orthviews('ZoomMenu');
            if(isfield(st,'zoomer'))
                cz = numel(zl)-get(st.zoomer,'Value')+1;
                if(cz<7)
                    CommonMethods.spm_setZoom(7);
                    CommonMethods.spm_setZoom(cz);
                end
            end
        end
        
        function out = spm_newImageLoaded(st,spmExt)
            out=false;
            if(isempty(spmExt))
                out=true;
                return;
            end;
            
            guiData=guidata(spmExt);
            spmStruct=guiData.spmStruct;
            if(isempty(spmStruct))
                out=true;
                return;
            end
            axis1=guiData.spmStruct.vols{1}.ax{1};
            V=st.vols{1};
            if(~isempty(V))
                if(isfield(V,'ax'))
                    axis2=st.vols{1}.ax{1};
                    if(axis1.ax~=axis2.ax) out=true; end;
                else
                    out=true;
                end
            else
                out=true;
            end
        end
        
        function out = isFigureHandle(h)
            if(isempty(h))
                out=false;
                return;
            end
            out=ishghandle(h);
            %            fig = h-rand;%commented out on 4/15/2016
            %             if isscalar(h) &&  ishghandle(h)
            %                 fig = getParentFigure(h);
            %             end
            %            out=fig==h;
        end
        
        function ind = getCurrentVolumeIndex(st)
            %This function determines the current volume index by comparing
            %the axis handles of the tree views of each volume with the
            %current axis.
            %It will return -1 if there is no volumes in st or the current
            %axis is not one of the orthogonal views of the volumes.
            ind=-1;
            if(isempty(st))
                return;
            end
            fig=st.fig;
            figure(fig);
            axc=gca;
            vols=st.vols;
            for i=1:length(vols)
                vol=vols{i};
                if(~isfield(vol,'ax'))
                    continue;
                end
                ax=vol.ax;
                for j=1:length(ax)
                    if(axc==ax{j}.ax)
                        ind=i;
                        return;
                    end
                end
            end
        end
        
        function spm_setCurrentVolumeIndex(st, ind)
            vol=st.vols{ind};
            if(~isempty(vol))
                axc=vol.ax;
                if(~isempty(axc))
                    axes(axc{1}.ax);
                end
            end
        end
        
        function [imgs, fnames]=getOpenImages(st)
            imgs={};
            fnames={};
            if(~isfield(st,'vols'))
                return;
            end
            vols=st.vols;
            for i=1:length(vols)
                vol=vols{i};
                if(~isempty(vol))
                    [path,name,ext]=fileparts(vol.fname);
                    imgs{end+1}=name;
                    fnames{end+1}=vol.fname;
                end
            end
        end
        
        function buildHistogram()
            %            fname='/home5/data/process/Prodrome_sf_VBM/Group/Outliers_AC_PC_Aligned/pValues.csv';
            fname='/home5/data/process/Prodrome_sf_VBM/Group/Outliers/pValues_withoutModulation_15311.csv';
            data=ImportData(fname);
            data=data.data;
            %            testData=data.textdata;
            RHO=data(:,1);
            figure;
            hist(RHO,25);
        end
        
        function valid = ValidAuxiliaryFile(file1, file2)
            valid=false;
            if(~exist(file1,'file')||~exist(file2,'file'))
                return;
            end
            f1=dir(file1);
            f2=dir(file2);
            valid=f2.datenum > f1.datenum;
        end
        
        function peak = findMaxNeighbor(vol, pos0, sign)
            persistent strel;
            dim=size(vol);
            if(length(dim)~=3)
                peak=pos0;
                return;
            end
            if(isempty(strel))
                strel=cell(1,26);
                num=0;
                for x=-1:1
                    for y=-1:1
                        for z=-1:1
                            if(x==0&&y==0&&z==0)
                                continue;
                            end;
                            num=num+1;
                            strel{num}=[x y z];
                        end
                    end
                end
                %put [0 0 0] to the first position;
                temp=strel{14};
                strel{14}=strel{1};
                strel{1}=temp;
            end
            
            peak=pos0;
            valP=sign*vol(peak(1), peak(2), peak(3));
            vals=zeros(1,26);
            for i=1:26
                pos=pos0+strel{i};
                if(any(pos<=0)||any(pos>dim))
                    continue;
                end
                val=sign*vol(pos(1),pos(2),pos(3));
                vals(i)=val;
                if(val>valP)
                    peak=pos;
                    valP=val;
                end
            end
        end
        
        function [Maxima,MaxPos,Minima,MinPos]=localExtrema_clusterCenters(vol, connectivity)
            if(exist('connectivity','var'))
                [~,MaxPos0,~,MinPos0]=CommonMethods.localExtrema(vol,connectivity);
                %MaxPos0 is a mx3 matrix
            else
                [~,MaxPos0,~,MinPos0]=CommonMethods.localExtrema(vol);
            end
            dim=size(vol);
            MaxPos1=CommonMethods.findClusterCenters(MaxPos0')';
            idx=sub2ind(dim,MaxPos1(:,1), MaxPos1(:,2), MaxPos1(:,3));
            Maxima=vol(idx);
            
            MinPos1=CommonMethods.findClusterCenters(MinPos0')';
            idx=sub2ind(dim,MinPos1(:,1), MinPos1(:,2), MinPos1(:,3));
            Minima=vol(idx);
            
            [Maxima, idx]=sort(Maxima, 'descend');
            MaxPos=MaxPos1(idx,:);
            
            [Minima, idx]=sort(Minima, 'ascend');
            MinPos=MinPos1(idx,:);
        end
        
        function rank=getRanking(Vs,val)
            %find the ranking of v in the vector Vs
            dV=abs(Vs-val);
            [~,rank]=min(dV);
        end
        function [Maxima,MaxPos,Minima,MinPos]=localExtrema(vol,connectivity)
            dim=size(vol);
            nanPos=isnan(vol);
            vol(nanPos)=0;
            if(exist('connectivity','var'))
                lx=imregionalmax(vol,connectivity);
                ln=imregionalmax(-vol,connectivity);
            else
                lx=imregionalmax(vol);
                ln=imregionalmax(-vol);
            end
            if(nnz(nanPos)>0)
                lx(nanPos)=0;
                ln(nanPos)=0;
            end
            
            idxx=find(lx);
            idxn=find(ln);
            Maxima0=vol(idxx);
            Minima0=vol(idxn);
            
            [Maxima, idx]=sort(Maxima0, 'descend');
            idxx=idxx(idx);
            MaxPos=CommonMethods.ind2sub_expanded(dim,idxx);
            
            [Minima, idx]=sort(Minima0, 'ascend');
            idxn=idxn(idx);
            MinPos=CommonMethods.ind2sub_expanded(dim,idxn);
            %MinPos and MaxPos are mx3 matrixes, m is the number of voxels
        end
        
        function subs=ind2sub_expanded(dim, idx)
            len=length(dim);
            
            comnd='[sub1';
            for i=2:len
                comnd=[comnd ', sub' num2str(i)];
            end
            comnd=[comnd '] = ind2sub(dim, idx);'];
            eval(comnd);
            subs=zeros(length(idx),len);
            
            for i=1:len
                eval(['subs(:,i)= sub' num2str(i) ';']);
            end
        end
        
        function [peak, poss] = findPeak(vol, pos, sign)
            peak=CommonMethods.findMaxNeighbor(vol,pos,sign);
            poss=pos;
            while(~all(peak==pos))
                pos=peak;
                poss=[poss;pos];
                peak=CommonMethods.findMaxNeighbor(vol,pos,sign);
            end
        end
        
        function spm_check_registration(imgs)
            if(ischar(imgs))
                imgs={imgs};
            end
            imgs=cellfun(@(x) spm_vol(x),imgs,'UniformOutput',false);
            imgsm=[];
            for i=1:length(imgs)
                imgsm=[imgsm;imgs{i}];
            end
            spm_check_registration(imgsm);
        end
        
        function nums = str2nums_decimal(str)
            len=length(str);
            spliter=';';
            for i=1:len
                if(~CommonMethods.isPartOfNum_decimal(str,i,len))
                    str(i)=spliter;
                end
            end
            nums=CommonMethods.str2nums(str, spliter);
        end
        
        function line=cell2line(strs,conn)
            %connect cell array of str with conn
            line='';
            if(isempty(strs))
                return;
            end
            if ~exist('conn','var')
                conn='_';
            end
            line=strs{1};
            for i=2:length(strs)
                line=[line conn strs{i}];
            end
        end       
        function groups = str2numGroups_decimal(str, spliter)
            strs=strsplit(str, spliter);
            groups=cellfun(@(x)CommonMethods.str2nums_decimal(x),strs, 'UniformOutput',false);
            groups=groups(arrayfun(@(x)~isempty(x), groups));
        end
        
        function [nums, numerical, strs] = str2nums(str, spliter)
            %this function extract numbers in a string
            if(~ischar(str))
                nums=[];
                numerical=[];
                strs={};
                return;
            end
            if(iscell(spliter))
                strs=strSplit_CellArray(str, spliter);
            else
                strs=strsplit(str, spliter);
            end
            numerical=cellfun(@(x) CommonMethods.isnum(x), strs);
            numstrs=strs(numerical);
            nums=cellfun(@(x) str2num(x), numstrs);
            strs=strs(~numerical);
        end
        
        function IS = isnum(str)
            %this function checks whether str is a number
            IS=true;
            try
                num=str2num(str);
            catch
                IS=false;
                return;
            end
            if(isempty(num))
                IS=false;
            end
        end
        
        function properties = getChildrenProperties(h, name)
            children = get(h, 'Children');
            properties=arrayfun(@(x) CommonMethods.getProperty(x, name), children, 'UniformOutput', false);
        end
        
        function prop = getProperty(h, name)
            prop = [];
            if(CommonMethods.isFigureHandle(h))
                props= get(h);
                if(isfield(props, name))
                    prop = get(h, name);
                end
            end
        end
        
        function IS=isdigit(c)
            IS=false;
            if(length(IS)==1)
                IS=(c>='0')&&(c<='9');
            end
        end
        
        function IS = isPartOfNum_decimal(str,p, len)
            IS=false;
            c=str(p);
            if(CommonMethods.isdigit(c))
                IS=true;
            elseif(strcmp(c,'-')||strcmp(c,'.'))
                if(p<len)
                    IS=CommonMethods.isdigit(str(p+1));
                end
            end
        end
        
        function [signals, pc, V] = getPCs(data)
            [m, n]   = size(data);
            %m is the number of measurements (data points, trials), and n is the number of the
            %measurement types (variables, dimension)
            %data: the input data
            %signals: mXn matrix of projected data
            %pc: each column is a PC, the eighen vector of the covarience
            %matrix
            %V: nX1 matrix of Variance
            ave=mean(data);
            ave=repmat(ave, m, 1);
            y=(data-ave);
            
            [u, s, pc]= svd(y);
            s=diag(s)/sqrt((m-1));
            V=s.*s;
            signals=y*pc;
        end
        
        function V = spm_setImageType(V0, type)
            V=V0;
            switch (type)
                case 'unit8'
                    V.dt(1) = 2;
                case 'int16'
                    V.dt(1) =4;
                case 'int32'
                    V.dt(1)=8;
                case 'float32'
                    V.dt(1)=16;
                case 'foat64'
                    V.dt(1)=64;
                case 'int8'
                    V.dt(1)=256;
                case 'uint16'
                    V.dt(1)=512;
                case 'uint32'
                    V.dt(1)=768;
            end
        end
        
        function matching = spm_IsMatchingMask(V, Vm)
            matching=CommonMethods.spm_MatchingHeaders(V,Vm);
        end
        
        function matching = spm_MatchingHeaders(V1, V2)
            %this function checks whether the dim and mat of the two
            %heather match.
            matching=true;
            if(any(V1.dim~=V2.dim))
                matching=false;
                return;
            end
            tol=0.0000000001;
            if(any(any(abs(V1.mat-V2.mat))>tol))
                matching=false;
            end
        end
        
        function IS=duplicatedImgs(fname1,fname2)
            IS=false;
            if ~exist(fname1,'file')||~exist(fname2,'file')
                return;
            end
            V1=spm_vol(fname1);
            V2=spm_vol(fname2);
            if(~CommonMethods.equalMatrices(V1.mat,V2.mat))
                return;
            end
            if(~CommonMethods.equalMatrices(V1.dim,V2.dim))
                return;
            end
            vols1=spm_read_vols(V1);
            vols2=spm_read_vols(V2);
            n1=isnan(vols1);
            n2=isnan(vols2);
            if nnz(n1-n2)
                return;
            end
            if ~nnz(vols1(~n1)-vols2(~n2))
                IS=true;
            end
        end
        
        function eq = equalMatrices(M1, M2)
            %this function test whether two matrices are the same size and
            %the equal matrix elements.
            eq=true;
            dim1=size(M1);
            dim2=size(M2);
            if(length(dim1)~=length(dim2))
                eq=false;
                return;
            end
            if(any(dim1-dim2))
                eq=false;
                return;
            end
            if(any(any(M1-M2)))
                eq=false;
            end
        end
        
        function eq = sameDim(M1, M2)
            %this function test whether the dimensions of the two matrices
            %are the same
            eq=CommonMethods.equalMatrices(size(M1), size(M2));
        end
        
        function removeFile(dir, fileName)
            %this function removes a file after escape special characters
            %inthe file name.
            fname=fullfile(dir,fileName);
            charsToEscape={'(' ')'};
            if(exist(fname,'file'))
                fname=CommonMethods.escapeCharachters(fname, charsToEscape);
                system(['rm ' fname]);
            end
        end
        
        function newStr=escapeCharachters(str, charsToEscape)
            %this function escapes the characters specified in
            %charsToEscape from the string 'str' and return the newStr
            newStr=[];
            for i=1:length(str)
                s=str(i);
                if(any(ismember(charsToEscape,s)))
                    s=['\' s];
                end
                newStr=[newStr s];
            end
        end
        
        function [coeff, score, latent] = spm_princomp(y)
            
            %this is the part of spm_regions that computes "eigenvariate"
            %of the data matrix.
            [m n]   = size(y);
            
            temp=m;
            m=n;
            n=temp;
            
            if m > n
                [v s v] = svd(y'*y);
                s       = diag(s);
                v       = v(:,1);
                u       = y*v/sqrt(s(1));
            else
                [u s u] = svd(y*y');
                s       = diag(s);
                u       = u(:,1);
                v       = y'*u/sqrt(s(1));
            end
            d       = sign(sum(v));
            u       = u*d;
            v       = v*d;
            Y       = u*sqrt(s(1)/n);
            coeff=v;
            score=Y;
            latent=s;
        end
        
        function [app, st] = getSourceApp()
            %this function detect the type of appliction that the current
            %figure is associated with.
            st=[];
            h=gcf;
            app=[];
            name=get(h,'name');
            
            if(length(name)<3)
                return;
            end
            
            if(strcmpi(name,'volume'))
                return;
            end
            
            if(strcmpi(name(1:3),'spm'))
                if(strcmpi(name,'spm_extension'))
                    guiData=guidata(h);
                    app=guiData.interactingApp;
                else
                    strs=strsplit(name, ' ');
                    st=strs{1};
                    len=length(st);
                    if(len>3)
                        st=st(4:len);
                        if(CommonMethods.isnum(st))
                            app='SPM';
                        end
                    end
                end
            elseif(length(name)>=7)
                if(strcmpi(name(1:7),'xJView:'))
                    app='xjView';
                end
            end
        end
        
        function IS = startsWith(str1, str2)
            %this function tests whether str1 starts with str2
            IS=strfind(str1,str2)==1;
        end
        
        function info=updateXJViewInfo(h, info)
            name=get(h,'name');
            len=length(name);
            
            if(len>=7)
                if(strcmpi(name(1:7),'xjView:'))
                    info.name=name;
                    guiData=guidata(h);
                    info.sectionViewTargetFile=guiData.sectionViewTargetFile;
                    imageFileName=guiData.imageFileName;
                    info.imageFileName=imageFileName;
                    info.currentxyz=guiData.currentxyz;
                    info.handles=guiData;
                    info.xjViewFigureHandle=h;;
                else
                    return;
                end
            else
                return;
            end
            
            imageFileName=info.imageFileName;
            numImgs=length(imageFileName);
            imgs0=arrayfun(@(x) num2str(x), 1:numImgs, 'UniformOutput', false);
            if(isfield(info,'vols'))
                vols0=info.vols;
                if(length(vols0)==numImgs)
                    vols0=info.vols;
                    imgs0=cellfun(@(x) x.fname, vols0, 'UniformOutput',false);
                end
            end
            
            vols=cell(1,numImgs);
            for i=1:numImgs
                img=imageFileName{i};
                if(CommonMethods.startsWith(img,imgs0{i}))
                    V=vols0{i};
                else
                    V=spm_vol(img);
                end
                if(~isfield(V,'vol'))
                    V.vol=spm_read_vols(V);
                    V.MaxStr={'Local Maxima List not Built'};
                    V.MinStr={'Local Minima List not Built'};
                end
                
                mni=guiData.currentmni;
                mni=mni{i};
                renew=false;
                if(~isfield(V,'clusterHandler'))
                    renew=true;
                elseif(V.clusterHandler.getSize()~=size(mni,2))
                    renew=true;
                end
                
                if(renew)
                    V.clusterHandler=ClusterHandler(mni', V.mat, V.vol);
                end
                vols{i}=V;
            end
            %           info.vols=cell2mat(vols);
            info.vols=vols;
        end
        
        function st = mni2str(mni)
            st=['(' sprintf('%4.0f', mni(1)) sprintf('%4.0f', mni(2)) sprintf('%4.0f', mni(3)) ')'];
        end
        
        function plotSummary(hAxis,y0)
            %this function plots bar graph from y0. y0 should be a 2D matrix or a cell array of 2D matrix. 
             %each column of y is assumed to be a set of data. It computes
             %the mean and sem from each set of data and plots them.
           if(isempty(hAxis))
                hAxis=gca;
            end
            
            if(~iscell(y0))
                y0={y0};
            end
            m=[];
            sem=[];
            for i=1:length(y0)
                y=y0{i};
                m=[m mean(y)];
                len=size(y,1);
                sem=[sem std(y)/sqrt(len)];
            end
            
            axes(hAxis);
            hbar=bar(m);
            set(hbar,'facecolor',[0 0 1]);
            c=get(hbar,'facecolor');
            hold(hAxis,'on');
            
            herr=errorbar(m,sem);
            set(herr,'LineStyle','none');
            set(herr,'LineWidth',3);
            set(herr,'color',c);
            hold(hAxis,'off');
            set(hAxis,'fontsize',14);
            if(mean(m)<0)
                %               set(hAxis,'ydir','reverse');
            end
        end
        
        function boxplot_Cellarray(hAxis,y0)
            %each column of y is assumed to be a set of data
            if(isempty(hAxis))
                hAxis=gca;
            end
            
            if(~iscell(y0))
                boxplot(y0);
                return;
            end
            
            numEls=0;
            for i=1:length(y0)
                numEls=numEls+numel(y0{i});
            end
            
            data=zeros(1,numEls);
            group=zeros(1,numEls);
            
            numEls=0;
            for i=1:length(y0);
                iI=numEls+1;
                numEls=numEls+numel(y0{i});
                group(iI:numEls)=i;
                data(iI:numEls)=y0{i};
            end
            boxplot(data,group);
        end
        
        function cntr = findCenter(coors)
            %coors are 3D coordinate, 3xm.
            %finding the point in coors that is closest to the center of
            %gravity of coors.
            gc=mean(coors,2);%center of the gravity
            %dists=arrayfun(@(i)CommonMethods.getDist(gc,coors(:,i)),1:size(coors,2));
            dists1=coors(1,:)-gc(1);
            dists2=coors(2,:)-gc(2);
            dists3=coors(3,:)-gc(3);
            
            dists=dists1.*dists1+dists2.*dists2+dists3.*dists3;
            %             dn=CommonMethods.getDist(gc,coors(:,1));
            %             for i=2:size(coors,2)
            %                 dist=CommonMethods.getDist(gc,coors(:,i));
            %                 if(dist<dn)
            %                     dn=dist;
            %                     ci=i;
            %                 end
            %             end
            [~,pn]=min(dists);
            cntr=coors(:,pn);
        end
        function cluster_coors=findClusters(coors, connectivity)
            %This function groups coors into clusters according
            %coors are 3D coordinate, 3xm, m is the number of coordinates.
            %coors: voxel indexes as 3 by m matrix
            transposed=false;
            if(size(coors,1)~=3&& size(coors,2)==3)
                coors=coors';
                transposed=true;
            end
            if ~exist('connectivity','var')
                connectivity=26;
            end
            cl=spm_clusters(coors,connectivity);%cluster label of each coordinate
            clusters=unique(cl);
            cluster_coors=cell(1,length(clusters));
            %the index indicating the cluster of the voxels belong to
            lenc0=max(clusters);%the number of clusters
            %counting the cluster sizes
            sizes0=zeros(1,lenc0);
            indexes0=cell(1,lenc0);
            empty=zeros(1,lenc0);
            for i=1:length(clusters);
                indc=clusters(i);
                inds=find(cl==indc);
                cluster_coors{i}=coors(:,inds);
            end
        end
        function cluster_coors=getClusters_mask(volsm)
            %returns cell array of volume indexes. Each cell is for one
            %cluster.
            dim=size(volsm);
            inds=find(volsm);
            [x,y,z]=ind2sub(dim,inds);
            coors=[x y z]';%coors for spm_clusters need to be 3xn matrix, n is the number of voxels
            cl=spm_clusters(coors);
            clusters=unique(cl);
            
            cluster_coors=cell(1,length(clusters));
            %the index indicating the cluster of the voxels belong to
            lenc0=max(clusters);%the number of clusters
            %counting the cluster sizes
            empty=zeros(1,lenc0);
            for i=1:length(clusters);
                indc=clusters(i);
                inds=find(cl==indc);
                cluster_coors{i}=coors(:,inds);
                n=length(inds);
                if(n==0)
                    empty(i)=1;
                    continue;
                end
            end
            cluster_coors=cluster_coors(~empty);
        end
        function sumStructs=getClusterSummary(cluster_coors,vols,vRange)
            %vRange: this function also compute the number of points(coors)
            %whose values are in the range specified by "vRange"
            if(~exist('vRange','var'))
                vRange=[];
            end
            num=length(cluster_coors);
            sumStructs=cell(1,num);
            for i=1:num
                coors=cluster_coors{i};
                if(size(coors,1)~=3)
                    coors=coors';
                end
                len=length(coors);
                v=arrayfun(@(x)vols(coors(x,1),coors(x,2),coors(x,3)),1:len);
                m=min(10,len);
                p=arrayfun(@(x)1/x,m:-1:1);
                sumStruct.quant=quantile(v,p);
                sumStruct.mean=mean(v);
                sumStruct.median=median(v);
                sumStruct.mad=mad(v);
                sumStruct.std=std(v);
                sumStruct.max=max(v);
                sumStruct.min=min(x);
                sumStruct.num=len;
                [~,~,volume,~,edgelength] = minboundbox(coors(:,1),coors(:,2),coors(:,3),'v',3);
                box.v=volume;
                box.edgelength=edgelength;
                sumStruct.box=box;
                %[rotmat,cornerpoints,volume,surface,edgelength] = minboundbox(x,y,z,metric,level)
                if(~isempty(vRange))
                    sumStruct.pointsInRange=nnz((v<=vRange(1).*v<vRange(2))>0);
                else
                    sumStruct.pointsInRange=nan;
                end
                sumStructs{i}=sumStruct;
            end
        end
        function [centers, sizes, cluster_coors] = findClusterCenters_mask(mask)
            centers={};
            sizes=[];
            cluster_coors={};
            dim=size(mask);
            inds=find(mask>0);
            if(isempty(inds))
                return;
            end
            [x,y,z]=ind2sub(dim,inds);
            coors=[x y z]';
            [centers, sizes, cluster_coors] = CommonMethods.findClusterCenters(coors);
        end
        function [centers, sizes, cluster_coors] = findClusterCenters(coors)
            %This function groups coors into clusters according to 18
            %connectivity scheme
            %coors are 3D coordinate, 3xm, m is the number of coordinates.
            %coors: voxel indexes as 3 by m matrix
            %centers: cluster centers
            %cluster_inds: one d indexes of coors belongs to each cluster
            transposed=false;
            if(size(coors,1)~=3&& size(coors,2)==3)
                coors=coors';
                transposed=true;
            end
            
            cl=spm_clusters(coors);%cluster label of each coordinate
            clusters=unique(cl);
            cluster_coors=cell(1,length(clusters));
            %the index indicating the cluster of the voxels belong to
            lenc0=max(clusters);%the number of clusters
            %counting the cluster sizes
            sizes0=zeros(1,lenc0);
            indexes0=cell(1,lenc0);
            centers0=zeros(3,lenc0);
            empty=zeros(1,lenc0);
            for i=1:length(clusters);
                indc=clusters(i);
                inds=find(cl==indc);
                cluster_coors{i}=coors(:,inds);
                n=length(inds);
                sizes0(i)=n;
                indexes0{i}=inds;
                if(n==0)
                    empty(i)=1;
                    continue;
                end
                center=CommonMethods.findCenter(coors(:,inds));
                centers0(:,i)=center;
            end
            centers=centers0(:,~empty);
            sizes=sizes0(~empty);
            if(transposed)
                centers=centers';
                sizes=sizes';
                cluster_coors=cellfun(@(x)x', cluster_coors, 'UniformOutput', false);
            end
        end
        
        function svdComparison()
            clear;
            %           matfile='/home2/data/process/PATHD/P1085/P1085a/EMOID/stats_sarf_aCompCor_ARTMotion7_ARTReg/VOI_Amygdala_1.mat';
            %           matfile='/home2/data/process/PATHD/P1067/P1067a/EMOWM/stats_ar
            %           f_aCompCor/VOI_Amygdala_1.mat';
            matfile='/home/tle6/data/process/PTSD/PP001/PP001/EMOID/stats_arf_aCompCor/VOI_Amygdala_1.mat';
            %             matfile='/home2/data/process/PATHD/P1085/P1085a/EMOID/stats_sarf_aCompCor_ARTMotion7_ARTReg/VOI_VOI_spm_1.mat';
            load(matfile);
            data0=xY.y;
            
            dim=size(data0);
            iI=1;
            iF=dim(1);
            %             iI=10;
            %             iF=dim(1)-10;
            data=data0(iI:iF,:);
            
            [coeff, score] = princomp(data);
            d=sign(sum(coeff(:,1)));
            PC_matlab=d*score(:,1);
            
            ave=mean(data);
            ave=repmat(ave,size(data,1),1);
            [coeff_spm, score_spm] = CommonMethods.spm_princomp(data-ave);
            d=sign(sum(coeff_spm(:,1)));
            PC_spm=d*score_spm(:,1);
            
            m=mean(data,2);
            x=max(data')';
            PCs=[PC_matlab Y(iI:iF) PC_spm];
            signalNames={'Matlab' 'spm' 'spm\_demeaned'};
            %             PCs=[PC_matlab Y(iI:iF) m];
            %             signalNames={'Matlab' 'spm' 'mean'};
            %             PCs=[PC_matlab Y(iI:iF) PC_spm m x];
            %             signalNames={'Matlab' 'spm' 'spmSection' 'max' 'mean'};
            %             PCs=[PC_matlab Y(iI:iF)];
            %             signalNames={'Matlab' 'spm'};
            figure();
            CommonMethods.displayIndexedTimeSeries(PCs, signalNames,false);
            h=figure();
            CommonMethods.displayDataColumnComparison(h,PCs,signalNames);
            
            PCs=[coeff(:,1) coeff_spm(:,1)];
            norms=Vnorm(PCs,1);
            %              norms=repmat(norms,size(PCs));
            %              PCs=PCs./norms;
            signalNames={'Matlab' 'spm'};
            figure();
            CommonMethods.displayIndexedTimeSeries(PCs, signalNames,false);
            h=figure();
            CommonMethods.displayDataColumnComparison(h,PCs,signalNames);
            figure();
            plot(PCs);
        end
        
        function datam = centeringData(data, direction)
            if(~exist('direction', 'var'))
                direction=1; %subtracting the columnmeans
            end
            dim=size(data);
            ave=zeros(size(dim));
            switch (direction)
                case 1
                    m=mean(data);
                    ave=repmat(m,dim(1),1);
                case 2
                    m=mean(data,2);
                    ave=repmat(m,1,dim(2));
            end
            datam=data-ave;
        end
        
        function svdComparisonResPond()
            clear;
            datadir='/home/tle6/data/process/PTSD/PP001/PP001/EMOID/stats_arf_aCompCor';
            
            %matfile is the file saved by spm_regions.m
            matfile=fullfile(datadir,'VOI_Amygdala_1.mat');
            load(matfile);
            data=xY.y;
            
            [coeff, score] = princomp(data);
            d=sign(sum(coeff(:,1)));
            PC_matlab=d*score(:,1);
            Y2=xY.u*10*rand(1)-mean(xY.u)*rand(1);
            figure();
            CommonMethods.displayDataColumnComparison(gcf, [PC_matlab xY.u Y2], ...
                {'PC\_matlab', 'xY.u', 'xY.u*a + b'});
            
            Var_matlab=var(PC_matlab);
            Var_Y=var(Y);
            
            
            data=CommonMethods.centeringData(data);
            [coeff_spm, score_spm, latent_spm] = CommonMethods.spm_princomp(data);
            d=sign(sum(coeff_spm(:,1)));
            PC_spm=d*score_spm(:,1);
            Var_spm=var(PC_spm);
            
            PCs=[PC_matlab./sqrt(dim(1)) Y PC_spm];
            %            signalNames={'Matlab' 'Mean\_time\_series' 'spm\_Homogeneous'};
            signalNames={'matlab/sqrt(#columns)' 'spm' 'spm\_demeaned'};
            
            figure();
            CommonMethods.displayIndexedTimeSeries(PCs, signalNames,false);
            h=figure();
            CommonMethods.displayDataColumnComparison(h,PCs,signalNames);
            
            %test princomp.m with and without centering
            a=rand(100,50);
            [coeffa, scorea]=princomp(a);
            [coeffb, scoreb]=princomp_withoutCentering(a);
            
            figure();
            plot([scorea(:,1) scoreb(:,1)], 'linewidth', 2);
            h=legend('PC Score: Centering', 'PC Score: No Centering');
            set(h,'fontsize', 16);
            figure();
            CommonMethods.displayDataColumnComparison(gcf, [scorea(:,1) scoreb(:,1)], ...
                {'PC Score: Centering', 'PC Score: No Centering'});
        end
        
        function ind = mni2ind(coor, mat, fraction)
            %turning mni coordinate into volume indexes, assuming a single
            %voxel, as a column or a row
            if(length(coor)==4)
                coor=coor(1:3);
            end
            dim=size(coor);
            if(dim(1) < dim(2))
                coor=([coor 1])';
            else
                coor=[coor; 1];
            end
            ixyz = mat \ coor;
            if(~exist('fraction','var'))
                ind=round(ixyz(1:3)');
            else
                if(~fraction)
                    ind=round(ixyz(1:3)');
                else
                    ind=ixyz(1:3)';
                end
            end
        end
        
        function ind = mni2ind_cols(coor, mat,noRound)
            %turning mni coordinate into volume indexes, assuming a
            %three-row matrix, each column is a 3D coordinate.
            %noRount: do not round the ind to integers
            if(~exist('noRound', 'var'))
                noRound=true;
            end
            ind=[];
            if(isempty(coor))
                return;
            end
            if(size(coor,1)~=3)
                coor=coor';
            end
            coor=[coor;ones(1,size(coor,2))];
            ixyz = mat \ coor;
            if(noRound)
                ind=round(ixyz(1:3,:));
            else
                ind=ixyz(1:3,:);
            end
        end
        
        function dist = getDist(p1, p2)
            dist=0;
            p=p2-p1;
            for i=1:length(p1)
                dist=dist+p(i)*p(i);
            end
            dist=sqrt(dist);
        end
        
        function spacing = getVoxelSpacing()
            %            file='/home2/data/process/PATHD/GroupAnalysis/Rest_FA77/Connectivity/HEP(T2_T1)/Amygdala3d_d0_b_TD/general/spmT_0001.img';
            file='/home4/Projects/fBIRN3_RESTdata/ConnectivityGroup/HC-SZ_flexANOVA_Group_Site_meGroup_meSite/falphNW_swlarf_ttestZ/NW_allIEclust3_p01FWE.img';
            V=spm_vol(file);
            mat=V(1).mat;
            coor=zeros(4);
            coor(:,1)=mat*[0 0 0 1]';
            coor(:,2)=mat*[1 0 0 1]';
            coor(:,3)=mat*[0 1 0 1]';
            coor(:,4)=mat*[0 0 1 1]';
            
            spacing=arrayfun(@(x) CommonMethods.getDist(coor(1:3,1), coor(1:3,x)), 2:4);
        end
        
        function volIndex = getValidMatrixIndexes(dim, volIndex)
            %replacing out of limit indexes with the boundary values
            volIndex(volIndex>dim)=dim(volIndex>dim);
            lows=ones(1,length(volIndex));
            volIndex(volIndex<1)=lows(volIndex<1);
        end
        
        function out = spm_loadMat(varargin)
            %load parameters in a SPM.mat file
            if(~isempty(varargin))
                SPM=varargin{1};
                if(ischar(SPM))
                    path=fileparts(SPM);
                    a=load(SPM);
                    SPM=a.SPM;
                    SPM.filepath=path;
                end
            end
            
            scans=SPM.xY.P;
            %scans in the first level modeling are stored as a column array
            if(~iscell(scans))
                dim = size(scans);
                
                scan=scans(1,:);
                try
                    V=spm_vol(scan);
                    di=1;
                catch ME
                    scan=scans(:,1)
                    try
                        V=spm_vol(scan)
                        di=2;
                    catch ME2
                        rethrow(ME2)
                    end
                end
                
                temp=scans;
                scans=cell(1,dim(di));
                for i=1:dim(di)
                    if(di==1)
                        scans{i}=temp(i,:);
                    else
                        scans{i}=temp(:,i);
                    end
                end
            end
            
            if(~isfield(SPM, 'vols'))
                temp=cellfun(@(x) spm_read_vols(spm_vol(x)), scans, 'UniformOutput', false);
                len=length(temp);
                dim=[size(temp{1}) len];
                vols=zeros(dim);
                for i=1:len
                    vols(:,:,:,i)=temp{i};
                end
                SPM.vols=vols;
                V=spm_vol(scans{1});
                SPM.mat=V.mat;
            end
            %SPM.vols is turned into a 4D volume.
            
            dim=size(SPM.vols);
            dim=dim(1:3);
            out=SPM;
        end
        
        function snames = spm_getShortScanNames(fnames)
            %makes shorter version of scan names
            len=length(fnames);
            paths=cell(1,len);
            names=cell(1,len);
            exts=cell(1,len);
            for s = 1:length(fnames)
                if(iscell(fnames))
                    [path, name, ext] = fileparts(fnames{s});
                else
                    [path, name, ext] = fileparts(fnames(s,:));
                end
                paths{s}=path;
                names{s}=name;
                exts{s}=ext;
            end
            
            ext=exts{1};
            if(~isempty(strfind(ext,';')))%scans are different frames of the same image
                snames=names;
                for i=1:len
                    snames{i}=[names{i} exts{i}];
                end
            else
                snames=fnames;
            end
        end
        
        function combineFigures()
            % First, create 4  figures with four different graphs (each with  a
            
            % colorbar):
            
            figure(1)
            
            surf(peaks(10))
            
            colorbar
            
            figure(2)
            
            mesh(peaks(10))
            
            colorbar
            
            figure(3)
            
            contour(peaks(10))
            
            colorbar
            
            figure(4)
            
            pcolor(peaks(10))
            
            colorbar
            
            % Now create destination graph
            
            figure(5)
            
            ax = zeros(4,1);
            
            for i = 1:4
                
                ax(i)=subplot(4,1,i);
                
            end
            
            % Now copy contents of each figure over to destination figure
            
            % Modify position of each axes as it is transferred
            
            for i = 1:4
                
                figure(i)
                
                h = get(gcf,'Children');
                
                newh = copyobj(h,5)
                
                for j = 1:length(newh)
                    
                    posnewh = get(newh(j),'Position');
                    
                    possub  = get(ax(i),'Position');
                    
                    set(newh(j),'Position', [posnewh(1) possub(2) posnewh(3) possub(4)])
                    
                end
                
                delete(ax(i));
                
            end
            
            figure(5)
            
        end
        
        function arrangeSubplots(hF, dim, hFunc)
            %create subplot according to dim, assigns 'position' to the
            %existing axis.
            if(~exist('hFunc','var'))
                hFunc='';
            end
            cols=dim(2);
            rows=dim(1);
            subso=CommonMethods.getSubplotHandlers(hF);
            nSubso=length(subso);
            subs=cell(1,cols*rows);
            if(~exist('hFunc','var'))
                hFunc=[];
            end
            
            if(~CommonMethods.isFigureHandle(hF))
                return;
            end
            
            axesN=zeros(rows,cols);
            hN=figure();
            for i=1:nSubso
                if(i<=rows*cols)
                    axesN(i)=subplot(rows,cols,i);
                    position=get(axesN(i),'position');
                    set(subso{i}.hAxs,'position',position,'ButtonDownFcn', hFunc);
                    subs{i}=subso{i};;
                    subs{i}.position=position;
                    subs{i}.rows=dim(1);
                    subs{i}.cols=dim(2);
                    subs{i}.rank=i;
                else
                    delete(subso{i}.hAxs)
                end
            end
            close (hN);
            
            
            
            num=dim(1)*dim(2);
            for i=(nSubso+1):num
                ax=subplot(rows, cols, i);
                set(ax,'ButtonDownFcn',hFunc);
                subs{i}=SubplotHandler();
                subs{i}.hFig=hF;
                subs{i}.hAxs=ax;
                subs{i}.position=get(ax,'position');
                subs{i}.ButtonDownFcn=hFunc;
                subs{i}.rows=dim(1);
                subs{i}.cols=dim(2);
                subs{i}.rank=i;
            end
            
            guiData=guidata(hF);
            guiData.subHandlers=subs;
            guiData.subDim=dim;
            guidata(hF, guiData);
        end
        
        function makeSubplots(hFig, subs)
            Axs=findall(hFig,'type','axes');
            for a=1:length(Axs)
                delete(Axs(a));
            end
            figure(hFig);
            
            subsN=cell(1,length(subs));
            for s=1:length(subs)
                sub=subs{s};
                subN=sub.clone();
                subN.hFig=hFig;
                hAxs=sub.hAxs;
                hAxsN=[];
                for a=1:length(hAxs)
                    axis=axes('position', subN.position, 'ButtonDownFcn', subN.ButtonDownFcn);
                    hAxsN=[hAxsN axis];
                end
                subN.hAxs=hAxsN;
                subN.hFig=hFig;
                subN.plot();
                subsN{s}=subN;
            end
            guiData=guidata(hFig);
            guiData.subHandlers=subsN;
            guidata(hFig, guiData);
        end
        
        function app = getCreatorApp(hFig)
            %find the label of the app that clreated the figure. The label
            %is stored as a field in guidata by the creating app.
            app=[];
            if(ishandle(hFig))
                guiData=guidata(hFig);
                if(isfield(guiData,'creatorApp'))
                    app=guiData.creatorApp;
                end
            end
        end
        
        function subs= getSubplots(hFig)
            %gets a cell array of structure describing all axes in hFig.
            %The the axes are created using subplot, then the structures
            %will also get the row and collumn of the axes.
            if(ishandle(hFig))
                subplots=findall(hFig,'Type','Axes');
                
                len=length(subplots);
                left=zeros(1,len);
                bottom=zeros(1,len);
                size=zeros(1,len);
                titles=cell(1,len);
                
                for i=1:len
                    pos=get(subplots(i), 'position');
                    left(i)=pos(1);
                    bottom(i)=pos(2);
                    size(i)=min(pos(3),pos(4));
                    title=get(subplots(i),'title');
                    titles{i}=get(title,'string');
                end
                rows=unique(bottom);
                cols=unique(left);
                
                subs=cell(1,len);
                num=0;
                noRows=length(rows);
                noCols=length(cols);
                for r=1:length(rows)
                    for c=1:length(cols)
                        for i=1:len
                            ax=subplots(i);
                            pos=get(ax, 'position');
                            if(pos(2)==rows(noRows-r+1)&&pos(1)==cols(c))
                                num=num+1;
                                sub=struct;
                                sub.ax=ax;
                                sub.row=r;
                                sub.col=c;
                                sub.title=titles{i};
                                sub.YALocation=get(ax, 'YAxisLocation');
                                sub.description=[num2str(r) '-' num2str(c) '-' sub.YALocation];
                                sub.rows=noRows;
                                sub.cols=noCols;
                                subs{num}=sub;
                            end
                        end
                    end
                end
            end
        end
        
        function label = makeFigureLabel(labelo, len)
            %this function replace '_' with '\_', and get the first len
            %letters followed by '...'
            %             if(~isempty(strfind(label0,'_')))
            %                 label=strrep(label0, '_', '\_');
            %             else
            %                 label=label0;
            %             end
            if(iscell(labelo))
                label=cellfun(@(x) CommonMethods.makeFigureLabel(x), labelo, 'UniformOutput', false);
                return;
            end
            
            label=labelo;
            inds=strfind(label,'_');
            num=0;
            for i=1:length(inds)
                ind=inds(i);
                replace=true;
                if(ind > 1)
                    if(strcmp(labelo(ind-1),'\'))
                        replace=false;
                    end
                end
                if(replace)
                    label=[label(1:ind+num-1) '\' label(ind+num:end)];
                    num=num+1;
                end
            end
            if(exist('len','var'))
                if(len>0&&len<(length(label)-3))
                    label=[label(1:len+num) '...'];
                end
            end
        end
        
        function axs=getAllAxes(hFig)
            %
            if(~exist('hFig', 'var'))
                hFig=gcf;
            end
            axs=findall(hFig,'type', 'axes');
        end
        
        function IS = IsOverlapping(position1, position2)
            %
            if(ishandle(position1))
                %some old way of callthing this function
                position1=get(position1,'position');
            end
            if(ishandle(position2))
                %some old way of callthing this function
                position2=get(position2,'position');
            end
            
            IS=CommonMethods.all_R(position1==position2);
        end
        
        function IS = IsOverlappingPositions(pos1, pos2)
            IS=CommonMethods.all_R(pos1==pos2);
        end
        
        function sub = getOverlappingSubplot(pos, subs)
            sub={subs{find(cellfun(@(x) CommonMethods.IsOverlappingPositions(pos,x.position),subs))}};
            if(~isempty(sub))
                sub=sub{1};
            else
                sub=[];
            end
        end
        
        function displayAxes(hA)
            %This function displays an axes hA in a new figure
            hF=figure
            hA1=subplot(1,1,1);
            p=get(hA1,'position');
            delete(hA1);
            copyobj(hA,hF);
            set(gca,'position',p);
        end
        
        function axs = getOverlappingAxes(hFig, position)
            %
            if(ishandle(position))
                %some old way of callthing this function
                position=get(position,'position');
            end
            axs=CommonMethods.getAllAxes(hFig);
            axs=axs(find(arrayfun(@(x) CommonMethods.IsOverlapping(position,x), axs)));
        end
        
        function purgeOverlappingAxes(hFig, position)
            %
            if(ishandle(position))
                %some old way of callthing this function
                position=get(position,'position');
            end
            axs=CommonMethods.getOverlappingAxes(hFig,position);
            for i=1:length(axs)
                axis=axs(i);
                set(axis,'xticklabel',{[]});
                set(axis,'xtick',[])
                set(axis,'yticklabel',{[]});
                set(axis,'ytick',[])
                delete(axis);
            end
        end
        
        function IS = all_R(m)
            %return a single 1 or zero regardless of the size of m
            IS=all(m);
            while(~isscalar(IS))
                IS=all(IS);
            end
        end
        
        function IS = any_R(m)
            %return a single 1 or zero regardless of the size of m
            IS=any(m);
            while(~isscalar(IS))
                IS=any(IS);
            end
        end
        
        function unifyYScales(Axs)
            for i=1:length(Axs)
                ylim(Axs(i),'auto');
            end
            refresh(get(Axs(1),'parent'));
            ylims=arrayfun(@(x) get(x,'ylim'), Axs,'UniformOutput', false);
            ylims=cell2mat(ylims);
            nYlims=[min((min(ylims))) max(max(ylims))];
            for i=1:length(Axs)
                axes(Axs(i));
                ylim(nYlims);
            end
        end
        
        function unifyXScales(Axs)
            for i=1:length(Axs)
                ylim(Axs(i),'auto');
            end
            refresh(get(Axs(1),'parent'));
            xlims=arrayfun(@(x) get(x,'xlim'), Axs,'UniformOutput', false);
            xlims=cell2mat(xlims);
            nXlims=[min((min(xlims))) max(max(xlims))];
            for i=1:length(Axs)
                axes(Axs(i));
                xlim(nXlims);
            end
        end
        
        function sphs = getSubplotHandlers(h)
            sphs=[];
            if(ishandle(h))
                guiData=guidata(h);
                if(isfield(guiData,'subHandlers'))
                    sphs=guiData.subHandlers;
                end
            end
        end
        
        function pos=findstrs(str, patterns)
            pos=cellfun(@(x)findstr(x,str), patterns, 'UniformOutput', false);
        end
        
        %         function nstr=strrep_cell(st, ostrs, nstrs)
        %             for i=1:length(ostrs)
        %                 st=strrep(st, ostrs{i}, nstrs{i});
        %             end
        %             nstr=st;
        %         end
        
        function code = getInequalityCode(st1, st2)
            %standard ineqality lb st1 sym st2 hb,
            if(strcmp(st1,'<')&&strcmp(st2,'<'))
                code=1;
            elseif(strcmp(st1,'<=')&&strcmp(st2,'<'))
                code=2;
            elseif(strcmp(st1,'<')&&strcmp(st2,'<='))
                code=3;
            elseif(strcmp(st1,'<=')&&strcmp(st2,'<='))
                code=4;
            else
                error('Wrong type of ineqality operator(s)');
            end
        end
        
        function valid=validInequality(lb,hb,num,code)
            switch code
                case 1
                    valid=num>lb&&num<hb;
                case 2
                    valid=num>=lb&&num<hb;
                case 3
                    valid=num>lb&&num<=hb;
                case 4
                    valid=num>=lb&&num<=hb;
                otherwise
                    error('invalide inequality code of CommonMethods.validInequality(lb,hb,num,code)');
            end
        end
        
        function [lb, st1, sym, st2, hb] = parseInequality(ieq)
            %this function break an inequality into lower bound (lb),
            %symbol, and higher bound. st is '<' or '<=', gt is '>' or '>='
            %the output are parameter for a standard ineqality
            %lb st1 sym st2 hb, st1 and st2 are '<' or '<='
            
            patterns={'<=' '<' '>=' '>'}; % the order of '<=' and '<" is important
            poss = CommonMethods.findstrs(ieq, patterns);
            
            if(~isempty([poss{1} poss{2}]))
                st=true;%smaller than, '<' or '<=';
            elseif(~isempty([poss{3} poss{4}]))
                st=false;
            else
                error('Invalid inequality');
            end
            
            ieq=CommonMethods.strrep_cell(ieq, patterns, {';' ';' ';' ';'});
            
            [nums,numerical,strs]=CommonMethods.str2nums(';',ieq);
            len=length(nums);
            
            lb=-inf;
            hb=inf;
            st1=[];
            st2=[];
            sym=trim(strs{1});
            
            switch len
                case 1 %there is only one number, eg 1 < x or x > 1
                    if(st)%e.g. '1<x, x<1 ...
                        if(numerical(1)) %the number is on the left side, e.g. 1<x or 1<=x
                            lb=nums(1);
                            if(~isempty(poss{1}))
                                st1='<=';
                            else
                                st1='<';
                            end
                        else %the number is on the right side, e.g. x< 1 or x<=1
                            hb=nums(1);
                            if(~isempty(poss{1}))
                                st2='<=';
                            else
                                st2='<';
                            end
                        end
                    else
                        if(numerical(1)) %e.g. 1>x or 1>=x
                            hb=nums(1);
                            if(~isempty(poss{3}))
                                st2='<=';
                            else
                                st2='<';
                            end
                        else %e.g. x>1 or x>=1
                            lb=nums(1);
                            if(~isempty(poss{3}))
                                st1='<=';
                            else
                                st1='<';
                            end
                        end
                    end
                case 2 %there are two numbers, eg 1 < x < 2 or 2 > x > 1
                    if(st) % e.g. 1 < x < 2
                        len1=length(poss{1});
                        lb=nums(1);
                        hb=nums(2);
                        switch len1
                            case 0
                                st1='<';
                                st2='<';
                            case 1
                                pos1=poss{1};
                                pos2=poss{2};
                                if(pos1(1)<=pos2(1))
                                    st1='<=';
                                    st2='<'
                                else
                                    st1='<';
                                    st2='<=';
                                end
                            case 2
                                st1='<=';
                                st2='<=';
                        end
                    else % 2 > x > 1
                        len1=length(poss{3});
                        lb=nums(2);
                        hb=nums(1);
                        switch len1
                            case 0
                                st1='<';
                                st2='<';
                            case 1
                                pos3=poss{3};
                                pos4=poss{4};
                                if(pos3(1)<=pos4(1))
                                    st2='<=';
                                    st1='<'
                                else
                                    st2='<';
                                    st1='<=';
                                end
                            case 2
                                st1='<=';
                                st2='<=';
                        end
                    end
            end
            if(lb==-inf)
                st1='<';
            end
            if(hb==inf)
                st2='<';
            end
        end
        
        function hc = getCreatorHandle(hFig)
            hc=[];
            if(CommonMethods.isFigureHandle(hFig))
                guiData=guidata(hFig);
                if(isfield(guiData,'creatorHandle'))
                    hc=guiData.creatorHandle;
                end
            end
        end
        
        function IS = strcmp_cell(st1, st2)
            IS=false;
            if(length(st1)==length(st2))
                if(size(st1,1)==size(st2,2))
                    st1=st1';
                end
                IS=all(strcmp(st1,st2));
            end
        end
        
        function IS = isAssigned(v)
            IS=true;
            if(~exist('v', 'var'))
                return;
            end
            
            IS=~isempty(v);
        end
        
        function name = getUserName ()
            if isunix()
                name = getenv('USER');
            else
                name = getenv('username');
            end
        end
        
        function home = getHome()
            home=getenv('HOME');
        end
        
        function [mcenter, morphDist]= getMorphCenter(volIndexes, dim0)
            %find the morphological center.
            %morphDist is the morphologocial distance of the voxels pointed
            %by volIndexes.
            len=length(volIndexes);
            coors=zeros(3,len);
            [coors(1,:),coors(2,:),coors(3,:)]=ind2sub(dim0,volIndexes);
            offset=zeros(3,1);
            maxs=zeros(3,1);
            
            for i=1:3
                offset(i)=min(coors(i,:))-1;
                maxs(i)=max(coors(i,:));
            end
            
            dim=(maxs-offset)';
            %subscripts and indexes of a 3d volume with size dim.
            coorst=coors-repmat(offset,1,len);
            volIndexest=sub2ind(dim, coorst(1,:),coorst(2,:),coorst(3,:));
            
            vol=zeros(dim);
            
            vol(volIndexest)=1;
            
            md=Shape_Handler(vol,true,false);
            morphDist=md.vol_MorphDist(volIndexest);
            
            mcenter=md.morphCenter+offset;
        end
        
        function [mcenter, vMorphDist]= getMorphCenter_vol(vol)
            md=Shape_Handler(vol,true,false);
            mcenter=md.morphCenter;
            vMorphDist=md.vol_MorphDist;
        end
        
        function fname = getCaseMatchingFileName(fname0)
            %this function find the existing file name if fname0 is a
            %correctfile name in a case insensitive os.
            fname=[];
            parts=strsplit(fname0, filesep);
            dir=filesep;
            len=length(parts);
            for i=1:len-1
                name0=parts{i};
                if(isempty(name0))
                    continue;
                end
                name=CommonMethods.getCaseMatchingFile(dir, name0);
                if(isempty(name))
                    return;
                end
                dir=fullfile(dir,name);
            end
            name0=parts(len);
            name=CommonMethods.getCaseMatchingFile(dir, name0);
            if(isempty(name))
                return;
            end
            fname=fullfile(dir,name);
        end
        
        function file = getCaseMatchingFile(path, fname0)
            lst=dir(path);
            file=[];
            for i=1:length(lst)
                name=lst(i).name;
                if(strcmpi(fname0, name))
                    file=name;
                    return;
                end
            end
        end
        
        function var = getBaseWorkspaceVariable(varName)
            %this function get the variable value from the base workspace
            if(evalin('base', ['exist(''' varName ''', ''var'')']))
                var=evalin('base',varName);
            else
                var=[];
            end
        end
        
        function assignToBaseWrokspace(name, value)
            assignin('base', name, value);
        end
        
        function varList = getFunctionVariables(Formular)
            id1=strfind(Formular,'(');
            id2=strfind(Formular,')');
            varList=strsplit(Formular(id1+1:id2-1), ',');
        end
        
        function IS = IsANewerFile(file1, file2)
            %true if file2 is newer
            if(~exist(file1,'file'))
                cprintf('red','%s\n', ['File1 ' file1 ' does not exist. Func: IsANewerFile']);
                IS=false;
                return;
            end
            if(~exist(file1,'file')||~exist(file2,'file'))
                cprintf('red','%s\n', ['File2 ' file2 ' does not exist. Func: IsANewerFile']);
                IS=false;
                return;
            end
            s1=dir(file1);
            s2=dir(file2);
            IS=s2.datenum>s1.datenum;
        end
        
        function IS = IsAnOlderFile(file, year, month, day)
            if(~exist(file,'file'))
                cprintf('red','%s\n', ['File ' file 'does not exist. Func: IsAnOlderFile']);
            end
            s=dir(file);
            IS=s.datenum<datenum(year,month,day);
        end
        
        function b = strstartswith(s, pat)
            %STRSTARTSWITH Determines whether a string starts with a specified pattern
            %
            %   b = strstartswith(s, pat);
            %       returns whether the string s starts with a sub-string pat.
            %
            
            %   History
            %   -------
            %       - Created by Dahua Lin, on Oct 9, 2008
            %
            
            %% main
            
            sl = length(s);
            pl = length(pat);
            
            b = (sl >= pl && strcmp(s(1:pl), pat)) || isempty(pat);
        end
        
        function [betas, columns]=getBetaNames_SPM(SPM)
            bnames=SPM.xX.name;
            
            betas={};
            cols=[];
            for i=1:length(bnames)
                name=bnames{i};
                parts=strsplit(name, ' ');
                name=parts{2};
                len=length(name);
                if(len<=5)
                    continue;
                end
                if(strcmp(name(len-5:len),'*bf(1)'))
                    beta=name(1:len-6);
                    betas{end+1}=beta;
                    cols=[cols i];
                end
            end
            betau=unique(betas);
            columns=cellfun(@(x) cols(find(ismember(betas,x))), betau,'UniformOutput', false);
            betas=betau;
        end
        
        function cropNativeRois(rois, maskfname, functionalFname)
            len=length(rois);
            for i=1:len
                roi=rois{i};
                [path, ~, ~]=fileparts(roi);
                mask=fullfile(path, [maskfname '.img']);
                if(~exist(mask,'file'))
                    createMask(functionalFname, mask);
                end
                volm=spm_read_vols(spm_vol(mask));
                Vr=spm_vol(roi);
                volr=spm_read_vols(Vr);
                volr(volm==0)=0;
                spm_write_vol(Vr,volr);
            end
        end
        
        function smoothImage(imgs, FWHM)
            %imgs: a cell array of images
            load(fullfile(CommonEnvMethods.getScriptPath(),'Jobtemplates', 'Smooth.mat'));
            if(ischar(imgs))
                imgs={imgs};
            end
            if(~exist('FWHM', 'var'))
                FWHM=[6 6 6];
            end
            matlabbatch{1}.spm.spatial.smooth.data=imgs;
            matlabbatch{1}.spm.spatial.smooth.fwhm=FWHM;
            spm_jobman('run', matlabbatch);
        end
        
        function pfnames=getPrefixedFiles(fnames, prefix)
            if(ischar(fnames))
                [path, name, ext]=fileparts(fnames);
                pfnames=fullfile(path,[prefix name ext]);
            else
                [paths, names, exts]=cellfun(@(x)fileparts(x),fnames, 'UniformOutput', false);
                pfnames=arrayfun(@(x)fullfile(paths{x}, [prefix names{x} exts{x}]), 1:length(fnames), 'UniformOutput', false);
            end
        end
        
        function fname=getSuffixedFilename(fname, suffix)
            [path, name, ext]=fileparts(fname);
            fname=fullfile(path, [name suffix ext]);
        end
        
        function [newfiles, dateStr] = getNewFiles( root, year, month, day )
            cutoff=datenum(year,month,day);
            addpath('/home/JinT/matlab_projects/fMRI/OperationalClasses/utilities');
            newfiles={};
            dateStr={};
            addNewFilesRecursively(root);
            function addNewFilesRecursively(directory)
                
                files=dir(directory);
                cprintf('blue','\n%s\n\n', ['Checking' directory ': ' num2str(length(files)-2) ' files.']);
                for i=1:length(files)
                    name=files(i).name;
                    if(strcmp(name,'.')||strcmp(name,'..'))
                        continue;
                    end
                    if(isdir(fullfile(directory,name)))
                        addNewFilesRecursively(fullfile(directory,name));
                    else
                        prop=dir(fullfile(directory, name));
                        if(length(prop)>1)
                            prop=prop;
                        end
                        if(prop.datenum>cutoff)
                            newfiles{end+1}=[fullfile(directory, name) '_' prop.date];
                            dateStr{end+1}=prop.date;
                        end
                    end
                end
            end
            %savepath('/home/JinT/pathdef.m');
        end
        
        function [idx, dateStr] = getNewerFileIndexes( files, year, month, day )
            cutoff=datenum(year,month,day);
            addpath('/home/JinT/matlab_projects/fMRI/OperationalClasses/utilities');
            props=cellfun(@(x) dir(x), files, 'UniformOutput', false);
            idx=find(cellfun(@(x) x.datenum>cutoff, props));
            props=props(idx);
            dateStr=cellfun(@(x) x.date, props,'UniformOutput', false);
        end
        
        function datenums=getFileDatenums(files)
            props=cellfun(@(x)dir(x),files,'UniformOutput',false);
            datenums=cellfun(@(x)x.datenum,props);
        end
        
        function dn=getFileDatenum(f)
            if(exist(f,'file'))
                prop=dir(f);
                dn=prop.datenum;
            else
                dn=nan;
            end
        end
        function st=getDatestr_fname(option)
            if ~exist('option','var')
                option='none'
            end
            st=datestr(clock());
            st=strrep(st,' ','_');
            st=strrep(st,':','-');
            if ~strcmp(option,'second')
                strs=strsplit(st,'_');
                st=strs{1};
            end
        end
        function [idx, dateStr] = getOlderFileIndexes( files, year, month, day )
            if ~exist('year','var')
                year=[];
            end
            if ~exist('month','var')
                month=1;
            end
            if ~exist('day','var')
                day=1;
            end
            
            if isempty(year)
                cutoff=now;
            elseif isempty(month)
                cutoff=datenum(year,month,day);
            end
            
            props=cellfun(@(x) dir(x), files, 'UniformOutput', false);
            idx=find(cellfun(@(x) x.datenum<cutoff, props));
            props=props(idx);
            dateStr=cellfun(@(x) x.date, props,'UniformOutput', false);
        end
        
        function IS=isEarlierFile(file, year, month, day)
            prop=dir(file);
            IS=prop.datenum<datenum(year,month,day);
        end
        
        function IS=IsNewerFile(file1, file2)
            if(~exist(file1,'file')||~exist(file2,'file'))
                IS=false;
                return;
            end
            prop1=dir(file1);
            prop2=dir(file2);
            IS=prop1.datenum>prop2.datenum;
        end
        
        %         function strs =getCommonStrings(strs1, strs2)
        %             %this function returns the common string elements of the two
        %             %string cell arrays strs1 and strs2
        %             strs=strs2(cellfun(@(x)any(ismember(strs1, x)), strs2))
        %         end
        
        function statst = seedStatExtraction(roi_imgs, data_imgs, fcpPars)
            %modified from \home\vpalzes\scripts\seedStatExtraction.m
            %1. to do functional connectivity processing before computing the
            %stats. 2. make it capable of extracting multiple seed ios
            %per data images. 3. the output is a cell array
            % fcpPars: functional connectivity processing parameters
            
            %--------------------------------------------------------------------------
            % seedStatExtraction
            %
            % Extracts various measures of a seed:
            %   - number of voxels
            %   - mean time series
            %   - summed total activity time series
            %   - coefficient alpha
            %
            % You will only be able to extract stats from one seed at a time with this
            % script, sorry!
            %
            % USAGE:
            % stats = seedStatExtraction(seed_imgs, data_imgs)
            %
            % INPUT:
            % seed_imgs:     cell array of reverse-normalized seed image filenames
            % data_imgs:    cell array of data image filenames for which the means and
            %               medians are extracted
            %
            % OUTPUT:
            % stats.mean:     matrix of mean activations
            % stats.sum:      matrix of summed activations
            % stats.n:        number of voxels in seed
            % stats.total_var:    variance of total activity time series
            % stats.time_var:     variance of each timepoint across all voxels
            % stats.vox_var:      variance of each voxel across all timepoints
            % stats.v_roi:        SPM V structure for the ROIs
            % stats.v_data:       SPM V structure for the data images
            %
            % NOTES:
            % Since our connectivity analyses now mainly use non-normalized,
            % non-smoothed data, you will need to reverse normalize your ROI before
            % running this script.
            %----------------------------------------------------------------------
            
            % Produce an error if the seed_imgs and data_imgs are not cell arrays
            if ~iscell(roi_imgs)
                error('Input variable "seed_imgs" must be a cell array!');
            end
            if ~iscell(data_imgs)
                error('Input variable "data_imgs" must be a cell array!');
            end
            
            % Produce an error if the number of seed_imgs and data_imgs do not match,
            % because each subject's data img should have a corresponding
            % reverse-normalized ROI img
            if length(roi_imgs)~=length(data_imgs)
                error('The size of "seed_imgs" and "data_imgs" should be the same!');
            end
            
            stats.int = [];
            
            % Get number of data imgs
            numSubs = length(data_imgs);
            
            % Get number of data points
            volume = spm_vol(data_imgs{1});
            numVols = length(volume);
            dim=volume(1).dim;
            numVoxels=dim(1)*dim(2)*dim(3);
            
            statst=cell(1,numSubs);
            % Initialize some variables
            
            for s = 1:numSubs
                
                cprintf('blue', 'Extracting %s of %s\n', num2str(s), num2str(numSubs));
                clear data vox_sum vox_var vox_var_sum sum_var k
                
                seed_imgs=roi_imgs{s};
                % Load the volumes
                
                numSeeds=length(seed_imgs);
                dataVol = spm_vol(data_imgs{s});
                % Load the data img
                dataVolY = spm_read_vols(dataVol);
                dataVolY_reshape = reshape(dataVolY,[numVoxels numVols]);
                stats.mean = nan(numVols, numSeeds);
                stats.sum = nan(numVols, numSeeds);
                stats.n0 = nan(numSeeds,1);
                stats.n = nan(numSeeds,1);
                stats.coeffalpha0 = nan(numSeeds,1);
                stats.coeffalpha = nan(numSeeds,1);
                
                for r=1:numSeeds
                    seedVol = spm_vol(seed_imgs{r});
                    
                    % Get seed name
                    [pathstr, seedName, seedExt] = fileparts(seedVol.fname);
                    seedName = [seedName seedExt];
                    
                    %%% xyz mm of non-zero points in seed img
                    %                    [Y, XYZ] = spm_read_vols(seedVol);
                    Y = spm_read_vols(seedVol);
                    idx = find(Y);
                    stats.n0(r)=length(idx);
                    %                   XYZ = [XYZ(:, idx); ones(1, length(idx))];
                    
                    % If the seed dimensions differ from the data img, then create error
                    % (NOTE: THIS IS BECAUSE YOU REVERSE NORMALIZED YOUR SEED, RIGHT?)
                    if any(dataVol(1).dim~=seedVol.dim)
                        error(['Volume dimensions of the seed differ from the data volume, remember you need to reverse normalize first!'...
                            '\nSeed = %s\nDataImg = %s\n'], seedVol.fname, dataVol.fname);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % coefficient alpha calculation
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    data=arrayfun(@(x)dataVolY_reshape(idx,x), 1:numVols,'UniformOutput',false);
                    data=cell2mat(data)';
                    data=CommonMethods.excludeNanColumns(data);
                    
                    stats.n(r)=size(data,2);%number of valid voxels
                    
                    data0=data;
                    data=CommonMethods.functionalConnectivityProcessing(data0, fcpPars{s});
                    
                    vox_sum0 = sum(data0, 2); %sum of all voxel at each time point.
                    vox_var0 = var(data0,0,1); % variance of each voxel, size should be the number of voxels
                    vox_var_sum0 = sum(vox_var0,2); % summed variance of each voxel
                    sum_var0 = var(vox_sum0,0,1); % variance of summed totals across voxels
                    
                    vox_sum = sum(data, 2); %sum of all voxel at each time point.
                    vox_var = var(data,0,1); % variance of each voxel, size should be the number of voxels
                    vox_var_sum = sum(vox_var,2); % summed variance of each voxel
                    sum_var = var(vox_sum,0,1); % variance of summed totals across voxels
                    
                    k = stats.n(r); % total number of voxels in seed
                    stats.coeffalpha0(r) = (k/(k-1.0)) * ((sum_var0 - vox_var_sum0) / sum_var0); % coefficient alpha
                    stats.coeffalpha(r) = (k/(k-1.0)) * ((sum_var - vox_var_sum) / sum_var); % coefficient alpha
                end
                statst{s}=stats;
            end
        end
        
        function statst = seedStatExtraction_testing(roi_imgs, data_imgs, fcpPars)
            %modified from \home\vpalzes\scripts\seedStatExtraction.m
            %1. to do functional connectivity processing before computing the
            %stats. 2. make it capable of extracting multiple seed ios
            %per data images. 3. the output is a cell array
            % fcpPars: functional connectivity processing parameters
            
            %--------------------------------------------------------------------------
            % seedStatExtraction
            %
            % Extracts various measures of a seed:
            %   - number of voxels
            %   - mean time series
            %   - summed total activity time series
            %   - coefficient alpha
            %
            % You will only be able to extract stats from one seed at a time with this
            % script, sorry!
            %
            % USAGE:
            % stats = seedStatExtraction(seed_imgs, data_imgs)
            %
            % INPUT:
            % seed_imgs:     cell array of reverse-normalized seed image filenames
            % data_imgs:    cell array of data image filenames for which the means and
            %               medians are extracted
            %
            % OUTPUT:
            % stats.mean:     matrix of mean activations
            % stats.sum:      matrix of summed activations
            % stats.n:        number of voxels in seed
            % stats.total_var:    variance of total activity time series
            % stats.time_var:     variance of each timepoint across all voxels
            % stats.vox_var:      variance of each voxel across all timepoints
            % stats.v_roi:        SPM V structure for the ROIs
            % stats.v_data:       SPM V structure for the data images
            %
            % NOTES:
            % Since our connectivity analyses now mainly use non-normalized,
            % non-smoothed data, you will need to reverse normalize your ROI before
            % running this script.
            %----------------------------------------------------------------------
            
            % Produce an error if the seed_imgs and data_imgs are not cell arrays
            if ~iscell(roi_imgs)
                error('Input variable "seed_imgs" must be a cell array!');
            end
            if ~iscell(data_imgs)
                error('Input variable "data_imgs" must be a cell array!');
            end
            
            % Produce an error if the number of seed_imgs and data_imgs do not match,
            % because each subject's data img should have a corresponding
            % reverse-normalized ROI img
            if length(roi_imgs)~=length(data_imgs)
                error('The size of "seed_imgs" and "data_imgs" should be the same!');
            end
            
            stats.int = [];
            
            % Get number of data imgs
            numSubs = length(data_imgs);
            
            % Get number of data points
            volume = spm_vol(data_imgs{1});
            numVols = length(volume);
            dim=volume(1).dim;
            numVoxels=dim(1)*dim(2)*dim(3);
            
            statst=cell(1,numSubs);
            % Initialize some variables
            numSeeds=4;
            for s = 1:numSubs
                
                cprintf('blue', 'Extracting %s of %s\n', num2str(s), num2str(numSubs));
                clear data vox_sum vox_var vox_var_sum sum_var k
                
                seed_imgs=roi_imgs{s};
                % Load the volumes
                
                dataVol = spm_vol(data_imgs{s});
                % Load the data img
                dataVolY = spm_read_vols(dataVol);
                dataVolY_reshape = reshape(dataVolY,[numVoxels numVols]);
                stats.mean = nan(numVols, numSeeds);
                stats.sum = nan(numVols, numSeeds);
                stats.n0 = nan(numSeeds,1);
                stats.n = nan(numSeeds,1);
                stats.coeffalpha0 = nan(numSeeds,1);
                stats.coeffalpha = nan(numSeeds,1);
                
                for r=1:numSeeds
                    seedVol = spm_vol(seed_imgs{1});
                    
                    % Get seed name
                    [pathstr, seedName, seedExt] = fileparts(seedVol.fname);
                    seedName = [seedName seedExt];
                    
                    %%% xyz mm of non-zero points in seed img
                    %                    [Y, XYZ] = spm_read_vols(seedVol);
                    Y = spm_read_vols(seedVol);
                    idx = find(Y);
                    stats.n0(r)=length(idx);
                    %                   XYZ = [XYZ(:, idx); ones(1, length(idx))];
                    
                    % If the seed dimensions differ from the data img, then create error
                    % (NOTE: THIS IS BECAUSE YOU REVERSE NORMALIZED YOUR SEED, RIGHT?)
                    if any(dataVol(1).dim~=seedVol.dim)
                        error(['Volume dimensions of the seed differ from the data volume, remember you need to reverse normalize first!'...
                            '\nSeed = %s\nDataImg = %s\n'], seedVol.fname, dataVol.fname);
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % coefficient alpha calculation
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    data=arrayfun(@(x)dataVolY_reshape(idx,x), 1:numVols,'UniformOutput',false);
                    data=cell2mat(data)';
                    data=CommonMethods.excludeNanColumns(data);
                    
                    stats.n(r)=size(data,2);%number of valid voxels
                    
                    cols=size(data,2);
                    rows=size(data,1);
                    
                    switch r
                        case 1
                            %data=data;
                        case 2
                            data=rand(size(data));
                        case 3
                            data=rand(size(data));
                            for c=2:cols
                                data(:,c)=data(:,1);
                            end
                        case 4
                            if(mod(cols,2)==1)
                                cols=cols+1;
                            end
                            data=rand(rows,cols);
                            for c=2:cols
                                sign=1;
                                if(mod(c,2)==0)
                                    sign=-1;
                                end
                                data(:,c)=sign*data(:,1);
                            end
                    end
                    
                    data0=data;
                    data=CommonMethods.functionalConnectivityProcessing(data0, fcpPars{s});
                    
                    vox_sum0 = sum(data0, 2); %sum of all voxel at each time point.
                    vox_var0 = var(data0,0,1); % variance of each voxel, size should be the number of voxels
                    vox_var_sum0 = sum(vox_var0,2); % summed variance of each voxel
                    sum_var0 = var(vox_sum0,0,1); % variance of summed totals across voxels
                    
                    vox_sum = sum(data, 2); %sum of all voxel at each time point.
                    vox_var = var(data,0,1); % variance of each voxel, size should be the number of voxels
                    vox_var_sum = sum(vox_var,2); % summed variance of each voxel
                    sum_var = var(vox_sum,0,1); % variance of summed totals across voxels
                    
                    k = stats.n(r); % total number of voxels in seed
                    stats.coeffalpha0(r) = (k/(k-1.0)) * ((sum_var0 - vox_var_sum0) / sum_var0); % coefficient alpha
                    stats.coeffalpha(r) = (k/(k-1.0)) * ((sum_var - vox_var_sum) / sum_var); % coefficient alpha
                end
                statst{s}=stats;
            end
        end
        
        function validdata=excludeNanColumns(data)
            %this function removes all collumns containing nans
            cols=size(data,2);
            validCols=find(arrayfun(@(x)~any(isnan(data(:,x))),1:cols));
            validdata=data(:,validCols);
        end
        
        function data=functionalConnectivityProcessing(data, pars)
            %This function performs functional connectivity processing
            %based on pars.
            bandpass_filter=pars.bandpass_filter;
            X = pars.desMat;
            repetitionTime=pars.repetitionTime;
            iX = pinv(X); %k by n
            y                                                  	                = data; %n by m matrix, m: number of voxels and n: the length of time series
            b                                                                   = iX*y; % each column is a beta vector for the corresponding voxel
            y                                                                   = y - X*b;
            data                                                                   = conn_filter(repetitionTime, bandpass_filter, y);
        end
        
        function test()
            fname='/home4/Projects/fBIRN3_RESTdata/ConnectivityGroup_RN_aCompCor_r2z/coefficientsOfAlpha/CoefficentOfAlpha_Thalamus.csv';
            datat=importdata(fname);
            data=datat.data;
            figure();
            subplot(2,2,1);
            hist(data(:,2),20);
            xlim([0.65 1]);
            ylim([0 120]);
            title('Thulamus (arf)', 'fontsize', 14);
            set(gca,'fontSize', 16);
            set(get(gca,'xlabel'),'string', 'Coeff. of Alph', 'fontsize', 14)
            set(get(gca,'ylabel'),'string', 'Counts', 'fontsize', 14)
            
            subplot(2,2,2);
            hist(data(:,3),20);
            xlim([0.65 1]);
            ylim([0 120]);
            title('Thulamus (aCompCor)', 'fontsize', 14);
            set(gca,'fontSize', 16);
            set(get(gca,'xlabel'),'string', 'Coeff. of Alph', 'fontsize', 14)
            set(get(gca,'ylabel'),'string', 'Counts', 'fontsize', 14)
            
            counts1=CommonMethods.getCumulativeDistributionCounts(data(:,2),0.5,1.,0.02);
            counts2=CommonMethods.getCumulativeDistributionCounts(data(:,3),0.5,1.,0.02);
            
            fname='/home4/Projects/fBIRN3_RESTdata/ConnectivityGroup_RN_aCompCor_ART_r2z/coefficientsOfAlpha/CoefficentOfAlpha_Thalamus.csv';
            
            datat=importdata(fname);
            data=datat.data;
            subplot(2,2,3);
            hist(data(:,3),20);
            xlim([0.65 1]);
            ylim([0 120]);
            title('Thulamus (aCompCor\_ART)', 'fontsize', 14);
            set(gca,'fontSize', 16);
            set(get(gca,'xlabel'),'string', 'Coeff. of Alph', 'fontsize', 14)
            set(get(gca,'ylabel'),'string', 'Counts', 'fontsize', 14)
            set(gca,'fontSize', 16);
            counts3=data(:,3);
            [counts3,x]=CommonMethods.getCumulativeDistributionCounts(data(:,3),0.5,1.,0.02);
            
            subplot(2,2,4);
            counts=[counts1;counts2;counts3]';
            plot(x',counts);
            set(get(gca,'xlabel'),'string', 'Coeff. of Alph', 'fontsize', 14)
            set(get(gca,'ylabel'),'string', 'Counts', 'fontsize', 14)
            
            
            figure
            subplot(2,1,1)
            plot(data0(:,[1 5]))
            title('Thulamus (arf)', 'fontsize', 14);
            subplot(2,1,2)
            plot(data(:,[1 5]))
            title('Thulamus (aCompCor)', 'fontsize', 14);
            
            %             subplot(2,2,4);
            %             hist(data(:,3),20);
            %             title('coeff. of alpha (fc processed)', 'fontsize', 18);
            %             xlim([0.65 1]);
            %             ylim([0 120]);
            %             set(gca,'fontSize', 18);
            %             X=[1.1 1.4 1.2 1.1];
            %             ex=rand(1,length(X));
            %             Y=[1.4 1.4 1.1];
            %             ey=rand(1,length(Y));
            %             a = [2 11] - 1;
            %
            %             bar((1:numel(X))+a(1), X, 'b');
            %             hold on;
            %             errorbar((1:numel(X))+a(1), X,  ex, 'b', 'LineStyle', 'None');
            %             hold on;
            %             bar((1:numel(Y))+a(2), Y, 'r');
            %             hold on;
            %             errorbar((1:numel(Y))+a(2), Y,  ey, 'r', 'LineStyle', 'None');
            %             hold off;
            %             set(gca,'XTickMode','auto');
            %             legend({'X','Y'})   ;
            %            fopen('test.txt','wt');
        end
        
        function [counts,x]=getCumulativeDistributionCounts(y,x1,x2,dx)
            n=ceil((x2-x1)/dx);
            counts=zeros(1,n);
            x=zeros(1,n);
            for i=1:n
                x(i)=x1+(i-1)*dx;
                counts(i)=length(find(y<=x1+(i-1)*dx));
            end
        end
        
        function str=getConstString(len,char)
            str=arrayfun(@(x) char(1), 1:len);
        end
        
        function str=leftfill(str0,width, char)
            str=str0;
            len=length(str0);
            if(len>width)
                return;
            end
            str=CommonMethods.getConstString(width,char);
            offset=width-len;
            
            for i=1:len;
                str(i+offset)=str0(i);
            end
        end
        
        function setImageType(img, type)
            V=spm_vol(img);
            vols=spm_read_vols(V);
            V.dt(1)=type;
            spm_write_vol(V,vols);
        end
        
        function imgs=splitMaskImage(img)
            %this function is to split each cluster in the mask image as a
            %severate image.
            %coors are 3D coordinate, 3xm, m is the number of coordinates.
            
            V=spm_vol(img);
            vols=spm_read_vols(V);
            idx=find(vols);
            [ind1, ind2, ind3]=ind2sub(V.dim, idx);%ind1, ind2 and ind3 are column vectors
            coors=[ind1';ind2';ind3'];
            cl=spm_clusters(coors);%cluster label of each coordinate
            clusters=unique(cl);
            len=length(clusters);
            
            imgs=cell(1,len);
            [path, name, ext]=fileparts(img);
            imgs=arrayfun(@(x) fullfile(path,[name '_' num2str(x) ext]), 1:len, 'UniformOutput', false);
            volsOne=ones(V.dim);
            volsZero=zeros(V.dim);
            for i=1:len
                idxt=idx(cl==clusters(i));
                volst=volsZero;
                volst(idxt)=volsOne(idxt);
                V.fname=imgs{i};
                spm_write_vol(V,volst);
            end
        end
        
        function img = mergeMaskImages(imgs, img)
            %this function does image addition and save to img
            len=length(imgs);
            V=spm_vol(imgs{1});
            vol=spm_read_vols(V);
            for i=2:len
                vol=vol+spm_read_vols(spm_vol(imgs{i}));
            end
            V.fname=img;
            spm_write_vol(V,vol);
        end
        
        function clearPaths()
            rmpath(genpath('/home/JinT/matlab_projects/fMRI/Scripts/utilities/wfu_bpm_beta'));
        end
        
        function line=strcell2line(strcell,linker)
            if(~exist('linker', 'var'))
                linker=', ';
            end
            line='';
            if(isempty(strcell))
                return;
            end
            line=strcell{1};
            for i=2:length(strcell)
                str=strcell{i};
                if(isempty(str))
                    str=' ';
                end
                line=[line linker str];
            end
        end
        
        function line=vec2line(vec,linker)
            if(~exist('linker', 'var'))
                linker=', ';
            end
            strcell=arrayfun(@(x) num2str(x), vec, 'UniformOutput', false);
            line=CommonMethods.strcell2line(strcell,linker);
        end
        
        function strcell=strrep_cell(strcell0, stro, strn)
            strcell=cellfun(@(x) strrep(x, stro, strn), strcell0, 'UniformOutput', false);
        end
        
        function [idx, str]=getSelection(handle, display)
            idx=[];
            str={};
            options=get(handle, 'String');
            idx=get(handle, 'Value');
            if(~exist('display','var'))
                display=false;
            end
            str=options(idx);
            if(display)
                line=[CommonMethods.vec2line(idx, '_') ': ' CommonMethods.strcell2line(str,': ')];
                disp(line);
            end
        end
        function [imgs,idx, ind]=getImgsValidAtTheVoxel(imgs0, coor)
            %this function select imgs whose value is not a nan at coor
            imgs={};
            idx=[];
            for i=1:length(imgs0)
                img=imgs0{i};
                V=spm_vol(img);
                vols=spm_read_vols(V);
                ind = CommonMethods.mni2ind(coor, V.mat);
                if(any(isnan(vols(ind(1),ind(2),ind(3),:))))
                    continue;
                end
                imgs{end+1}=img;
                idx(end+1)=i;
            end
        end
        function roi_stats_to_csv_maskIntersect(p, setup)
            %--------------------------------------------------------------------------
            % roi_stats_to_csv
            %
            % USAGE:
            %
            % roi_stats_to_csv(p, setup)
            %
            % takes the p structure from roi_stat_extraction and print to a CSV file.
            % The setup structure includes some additional information
            %
            % INPUT:
            %
            %  p: structure from roi_stat_extraction as a cell array:
            %   p{1} for the first contrast
            %   p{2} for the second contrast, etc.
            %
            % p{1}.contrastname
            % The p stucture should include the contrast name in the field contrastname
            % e.g. p{1}.contrastname.
            %
            % setup.filename: output filename (including path)
            %   setup.columns is addition columns to include in the CSV file, e.g.:
            %   setup.columns{1}.vector = subjectIDs;
            %   setup.columns{1}.name   = 'subjectID';
            %   setup.columns{2}.vector = group;
            %   setup.columns{2}.name   = 'group';
            %
            % KWJ, Dec 2009
            % Modified 7/7/15 by KBW to include PC eigenvalue
            %
            % see also: roi_stat_extraction
            
            
            % open output file
            fid = fopen(setup.filename ,'wt');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % print the header
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for extra=1:length(setup.columns)
                fprintf(fid,setup.columns{extra}.name);
                fprintf(fid,',');
            end
            
            if(isfield(setup,'dataSuffix'))
                suffix=setup.dataSuffix;
            else
                suffix={''};
            end
            
            for c=1:size(p,1)
                for i=1:length(p{c}.v_roi)
                    [dummy_pth file dummy_ext]=fileparts(p{c}.v_roi(i).fname);
                    
                    for(s=1:length(suffix))
                        fprintf(fid,[p{c,s}.contrastname '__mean_' suffix{s} '_' file]);
                        fprintf(fid,',');
                        
                        fprintf(fid,[p{c,s}.contrastname '__median_' suffix{s} '_' file]);
                        fprintf(fid,',');
                        
                        fprintf(fid,[p{c,s}.contrastname '__N_' suffix{s} '_' file]);
                        fprintf(fid,',');
                    end
                    
                    fprintf(fid,['ROISize__' file]);
                    fprintf(fid,',');
                    
                    fprintf(fid,['ROISizeO__' file]);
                    fprintf(fid,',');
                end
            end
            fprintf(fid,'\n');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % print the data
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for i=1:size(p{1}.mean_act,1)
                for c=1:size(p,1)
                    for ii=1:size(p{c}.mean_act,2)
                        if ii==1 && c==1
                            for extra=1:length(setup.columns)
                                fprintf(fid,setup.columns{extra}.vector{i});
                                fprintf(fid,',');
                            end
                        end
                        
                        % means
                        for(s=1:length(suffix))
                            fprintf(fid,num2str(p{c,s}.mean_act(i,ii)));
                            fprintf(fid,',');
                            
                            % medians
                            fprintf(fid,num2str(p{c,s}.median_act(i,ii)));
                            fprintf(fid,',');
                            
                            % n
                            fprintf(fid,num2str(p{c,s}.n_act(i,ii)));
                            fprintf(fid,',');
                            
                            % n
                        end
                        fprintf(fid,num2str(p{c}.nRoi_act(i,ii)));
                        fprintf(fid,',');
                        
                        % n
                        fprintf(fid,num2str(p{c}.nRoio_act(i,ii)));
                        fprintf(fid,',');
                    end
                end
                fprintf(fid,'\n');
            end
            
            % close output file
            fclose(fid);
        end
        
        function [gCoeff, dCoeff, p, o, po] = calVoxelICC(img1, img2, mask)
            %computing voxelwise g and d coefficient.
            addpath /home/kwj6/scripts/fbirn/ects/;
            vol1=spm_read_vols(spm_vol(img1));
            vol2=spm_read_vols(spm_vol(img2));
            if(~exist('mask','var'))
                volt=vol1+vol2;
                if(any(isnan(volt)))
                    idx=find(~isnan(volt));
                else
                    V=spm_vol(img1);
                    fname=V.fname;
                    [path,name,ext]=fileparts(fname);
                    V.fname=fullfile(path,[name '_temp' ext]);
                    V.dim=[V.dim 2];
                    volt=zeros(V.dim);
                    volt(:,:,:,1)=vol1;
                    volt(:,:,:,2)=vol2;
                    spm_write_vol(V,volt);
                    mask=fullfile(path,[name '_maskTemp' ext]);
                    cratemask(V.fname,mask);
                    volt=spm_read_vols(spm_vol(mask));
                    idx=find(volt);
                    delete(mask);
                    delete(V.fname);
                end
            else
                volt=spm_read_vols(spm_vol(mask));
                idx=find(volt);
            end
            
            scores=[vol1(idx) vol2(idx)];
            [vcs] = G(scores, 2, NaN, 1);
            p= vcs(1);
            o(x,y,z)= vcs(2);
            po(x,y,z)= vcs(3);
            gCoeff=p/(p+o);
            dCoeff=p/(p+o+po);
        end
        
        function mask=getSPM_brainmask()
            mask='/usr/local/spm8/apriori/brainmask.nii';
        end
        
        function elapsedTimeStr=getElapsedTimeStr(clock1, clock2)
            units={'Year' 'Months' 'Days' 'Hours' 'Minutes' 'Seconds'};
            elapsedTimeStr='';
            eT=clock2-clock1;
            for i=1:length(clock1)
                dt=eT(i);
                if(dt==0)
                    continue;
                end
                elapsedTimeStr=[elapsedTimeStr ' ' num2str(dt) ' ' units{i}];
            end
        end
        
        function line=displayElapsedTime(clock1, clock2)
            etstr=CommonMethods.getElapsedTimeStr(clock1, clock2);
            line=['Started at ' datestr(clock1) ' finished at ' datestr(clock2) ' (Elapsed time: ' etstr ').'];
            cprintf('blue', '%s\n', line);
        end
        
        function outputTimingLog(path, clock1, clock2, filename)
            if(~exist('filename','file'))
                filename='timingLog.txt';
            end
            fid=fopen(fullfile(path,filename),'wt');
            line=CommonMethods.displayElapsedTime(clock1, clock2);
            fprintf(fid,'%s\n', line);
            fclose(fid);
        end
        
        function [strdiff, idxd] = getStringCellarrayComparison(strs1, strs2)
            %This function returns cell arrays of string that are not
            %common to the two arrays
            idx1d=find(cellfun(@(x)~any(ismember(strs2, x)), strs1));
            idx2d=find(cellfun(@(x)~any(ismember(strs1, x)), strs2));
            strdiff={strs1(idx1d) strs2(idx2d)};
            idxd={idx1d idx2d};
        end
        
        function [ca, idxr] = removeCells(ca, idx)
            if(iscell(idx))
                idx=cell2mat(idx);
            end
            %ca = ca(idxrO)
            len1=length(ca);
            len2=length(idx);
            if(len2>len1)
                return;
            end
            
            os=ones(1,len1);
            os(idx)=zeros(1,len2);
            idxr=find(os);
            ca=ca(idxr);
        end
        
        function [ca, idxr] = removeElements(ca, idx)
            if(iscell(ca))
                [ca, idxr]=CommonMethods.removeCells(ca,idx);
            else
                tmp=ones(1,length(ca));
                tmp(idx)=0;
                idxr=find(tmp);
                ca=ca(idxr);
            end
        end
        
        function [ca, idx] = removeStringEles(strs1, strs2)
            %ca=strs1(idx)
            idx=find(cellfun(@(x)~any(ismember(strs2, x)), strs1));
            ca=strs1(idx);
        end
        
        function [commonStrs, idxc1, idxc2] = getCommonStrings(strs1, strs2)
            %This function returns cell arrays of string that are
            %common to the two arrays
            idxc1=find(cellfun(@(x)any(ismember(strs2, x)), strs1));
            idxc2=find(cellfun(@(x)any(ismember(strs2, x)), strs1));
            commonStrs=strs1(idxc1);
        end
        
        function cnames = getContrastNames(SPM)
            cnames=arrayfun(@(x) x.name, SPM.xCon, 'UniformOutput', false);
        end
        
        function [fnames, maskSizes] = getMaskSizes_nan(matfile)
            %/home4/data/process/Prodrome/E13003/20050103_05031311/GO/stats_arf_weighted_by_run/wmask.img
            %/home4/data/process/Prodrome/Group/newGOAnalysis/GO_z_age/Data/stats_arf_weighted_by_run/10trialmin/FirstEpisode/E13003__05031311_ngo-trg.img,1
            %_/home4/data/process/Prodrome/E33012/20070221_07033642/GO/stats_arf_weighted_by_run/wmask.img
            %/home4/data/process/Prodrome/Group/newGOAnalysis/GO_z_age/Data/stats_arf_weighted_by_run/10trialmin/FirstEpisode/E33012__07033642_ngo-trg.img,1
            
            if(~exist('matfile','var'))
                matfile='/home4/data/process/Prodrome/Group/newGOAnalysis/GO_z_age/Analysis/stats_arf_weighted_by_run/10trialmin/ANOVA_fullfact/ngo-trg/SPM.mat';
            end
            load(matfile);
            fnames=SPM.xY.P;
            vols=spm_read_vols(spm_vol(char(fnames)));
            maskSizes=arrayfun(@(x)nnz(~isnan(vols(:,:,:,x))), 1:length(fnames));
            figure
            subplot(2,2,1);
            hist(maskSizes)
            title('histogram and scattered plot of mask sizes (all subjects)');
            subplot(2,2,2);
            scatter(1:length(maskSizes),maskSizes);
            maskSizes=maskSizes(maskSizes>48000);
            subplot(2,2,3);
            hist(maskSizes)
            title('histogram and scattered plot of mask sizes (after remove two subjects)');
            subplot(2,2,4);
            scatter(1:length(maskSizes),maskSizes);
        end
        
        function plotVoxels()
            coor=[-39 -46 -44];
            matfile='/home4/Projects/fBIRN3_RESTdata/ConnectivityGroup_RN_aCompCor_ART/HC-SZ_flexANOVA_Group_Site_meGroup_meSite/Thalamus/SPM.mat';
            load(matfile);
            imgs0=SPM.xY.P;
            imgsn=CommonMethods.getPrefixedFiles(imgs0,'n');
            imgs0=char(imgs0);
            imgsn=char(imgsn);
            V0=spm_vol(imgs0);
            Vn=spm_vol(imgsn);
            V=V0(1);
            ind=CommonMethods.mni2ind(coor, V.mat);
            vols0=spm_read_vols(V0);
            volsn=spm_read_vols(Vn);
            vec0=squeeze(vols0(ind(1),ind(2),ind(3),:));
            vecn=squeeze(volsn(ind(1),ind(2),ind(3),:));
            indn=find(isnan(vecn));
            ind0=find(~isnan(vecn));
            figure();
            subplot(2,1,1);
            set(gca,'FontSize', 18);
            hist(vec0(ind0));
            xlim([-0.5 0.8]);
            title(gca, 'Histogram of Included Voxels');
            set(gca,'FontSize', 18);
            subplot(2,1,2);
            hist(vec0(indn));
            xlim([-0.5 0.8]);
            title(gca,'Histogram of Excluded Voxels');
            V=V;
        end
        
        function copyscripts()
            path1='/home/vpalzes/scripts/FordMerit';
            path2='/home/JinT/matlab_projects/fMRI/Scripts/Prodrome/GoNoGo/gPPI';
            files={...
                'FordMerit_gPPI_reverseNormalize_ROIs.m' ...
                'FordMerit_block_gPPIAnalysis.m' ...
                'FordMerit_event_gPPIAnalysis.m'...
                'FordMerit_block_gPPI_1stLevelContrasts.m' ...
                'FordMerit_GroupAnalysis_Block_gPPI.m' ...
                'FordMerit_MissingVoxel_GroupAnalysis_Block_gPPI.m' ...
                };
            arrayfun(@(x) copyfile(fullfile(path1,files{x}), fullfile(path2, files{x})), 1:length(files));
        end
        
        function h=removedSelection_LstBx(h)
            strs=get(h,'String');
            idx=get(h, 'Value');
            strs=CommonMethods.removeCells(strs,idx);
            set(h,'String',strs);
            set(h,'Value', max(0,length(strs)));
        end
        
        function folder = getParentFolder(folder, levels)
            parts=strsplit(folder, filesep, 'include');
            folder='';
            for i=1:length(parts)-levels
                folder=[folder parts{i}];
            end
        end
        
        function strs=searchStr_LstBx(h, str)
            strs=get(h,'String');
            if(isempty(str))
                return;
            end
            idx=cellfun(@(x)~isempty(strfind(x,str)), strs);
            strs=strs(idx);
        end
        
        function strs=searchStri_LstBx(h, str)
            strs=get(h,'String');
            if(isempty(str))
                return;
            end
            strsl=cellfun(@(x) lower(x), strs, 'UniformOutput', false);
            strl=lower(str);
            idx=cellfun(@(x)~isempty(strfind(x,strl)), strsl);
            strs=strs(idx);
        end
        function MakeSphereROI(img, center, radius)
            
            %--------------------------------------------------------------------------
            %Modified from
            % Name : /home/vpalzes/scripts/FordMerit/FordMerit_PeaktoSphereROI.m by
            % Vanessa Palzes
            % Creation Date : 11/18/2015
            %
            % Purpose : This script will make img a spere roi specified by shpere
            %
            % Inputs:
            %
            % img: the input image
            % shere: a structure with two fields, center and radius (in mm)
            %--------------------------------------------------------------------------
            
            sphere.centre=center;
            sphere.radius=radius;
            
            addpath /home2/software/marsbar_allVERSIONS/marsbar-0.43/
            
            %             spm('defaults', 'fmri')
            %             spm_jobman('initcfg');
            if(~exist('label', 'var'));
                label='';
            end
            % Load marsbar default options .mat file
            load /home/vpalzes/scripts/marsbaroptions.mat
            
            % Set the ROI space for MarsBaR to the contrsat space
            MARS.OPTIONS.spacebase.fname = img;
            mars_options('put');
            
            
            sphere_roi = maroi_sphere(sphere);
            
            % Label the ROI
            %sphere_roi = label(sphere_roi, label);
            
            % Save the .mat structure
            %saveroi(sphere_roi, fullfile(curdir, [name '.mat']));
            
            % Save the ROI mask
            save_as_image(sphere_roi, img);
            
            % Set the ROI space for MarsBaR back to default
            MARS.OPTIONS.spacebase.fname = '/usr/local/spm8/templates/T1.nii';
            mars_options('put');
            rmpath  /home2/software/marsbar_allVERSIONS/marsbar-0.43/
        end
        
        function roisn=makeToSingleVoxelMask(rois,affix)
            %this function will shring masks in rois into a single voxel
            %masks and change the data type into floating point.
            roisn=cell(1,length(rois));
            if~exist('affix','var')
                affix='';
            end
            for i=1:length(rois)
                roi=rois{i};
                [path,name,ext]=fileparts(roi);
                V=spm_vol(roi);
                vols=spm_read_vols(V);
                idx=find(vols>0);
                [x,y,z]=ind2sub(V.dim,idx);
                X=round(mean(x));
                Y=round(mean(y));
                Z=round(mean(z));
                volsn=zeros(V.dim);
                volsn(X, Y, Z)=1;
                V.dt(1)=16;
                V.fname=fullfile(path,[name '_SingleVoxel' affix ext]);
                spm_write_vol(V, volsn);
                roisn{i}=V.fname;
            end
        end
        
        function saveHihgerResolutionImage(img, factor, img_new)
            %this function increase the resolution of the image
            %this function is not verified yet.
            V=spm_vol(img);
            dimo=V.dim;
            volso=spm_read_vols(V);
            dim=dimo.*factor;
            vf=factor(1)*factor(2)*factor(3);
            vols=nan(dim);
            os=ones(factor);
            for x=1:dimo(1)
                for y=1:dimo(2)
                    for z=1:dimo(3)
                        iI=[factor(1)*(x-1)+1 factor(2)*(y-1)+1 factor(3)*(z-1)+1];
                        iF=[factor(1)*x factor(2)*y factor(3)*z];
                        t=volso(x,y,z);
                        vols(iI(1):iF(1), iI(2):iF(2), iI(3):iF(3))=t*os;
                    end
                end
            end
            %            arrayfun(@(x,y,z)vols((factor(1)*(x-1)+1):factor(1)*x, (factor(2)*(y-1)+1):factor(2)*y, (factor(3)*(z-1)+1):factor(3)*z), 1:dimo(1), 1:dimo(2), 1:dimo(3));
            V.fname=img_new;
            V.dim=dim;
            spm_write_vol(V,vols);
        end
        
        function jac=getMaskJaccardIndex(volsA, volsB, cutoff)
            %this function compute the jaccard index of between the two mask
            %images created using common cutoff. volsA and volsB are the
            %voxel values of two images of the same dimensionality.
            volsA=volsA>cutoff;
            volsB=volsB>cutoff;
            volsU=(volsA+volsB)>0;
            volsI=(volsA.*volsB)>0;
            nnzU=nnz(volsU);
            if(nnzU==0)
                jac=1;
                return;
            end
            jac=nnz(volsI)/nnzU;
        end
        
        function jac=getMaskJaccardIndex_unionMask(volsAs, volsBs, cutoff)
            %this function compute the jaccard index of between the two mask
            %images created using common cutoff. volsAs and volsBs are cell arrays containing the
            %voxel values of two group of images of the same dimensionality.
            maskA=volsAs{1}>cutoff;
            for i=2:length(volsAs)
                maskA=maskA+volsAs{i}>cutoff;
            end
            
            maskB=volsBs{1}>cutoff;
            for i=2:length(volsBs)
                maskB=maskB+volsBs{i}>cutoff;
            end
            
            volsU=(maskA+maskB)>0;
            volsI=(maskA.*maskB)>0;
            nnzU=nnz(volsU);
            if(nnzU==0)
                jac=1;
                return;
            end
            jac=nnz(volsI)/nnzU;
        end
        
        function overlayImages(underlayImg, overlayImgs, outfname)
            % By Taihao, adapted from /home/vpalzes/scripts/overlaycoregFSmasks.m by Kate and Vanessa
            %
            % Creation Date : 12/11/2015
            %
            % Purpose : Creates slice overlay .png files.
            %
            % underlayImg: the underlay image
            % overlayImgs: the overlay images
            % outfname: the name of the output png file name
            %
            % The modification made by Taihao:
            % 1. directly input the underlay and overly images as well as
            % the output put file path.
            % 2. addjust the display range based on the bounding box of the
            % underlay image.
            % 3. adjust the display intensity range of the images based on
            % the voxel values of the images instead of fixed ranges.
            %--------------------------------------------------------------------------
            
 %           load /home/kwj6/scripts/core_scripts/blobs_so.mat; % load slover object template
            
            % Clear out extra images in template
            SO = slover;
            a.vol='';
            SO.img=[a a];
            
            % Directory where the .png images will be saved
            if(~exist('outfname','var'))
                outfname='';
            end
            if(~isempty(outfname))
                outdir=fileparts(outfname);
                if ~isdir(fullfile(outdir))
                    mkdir(outdir);
                end
            end
            % Output filename
            filename = outfname;
            
            
            %%%% AXIAL VIEW
            SO.transform    = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
            slices          = 72:-4:-40;
            
            % Change images for overlay
            SO.img(1).vol   = spm_vol(underlayImg);
            range = CommonMethods.getMatValueRange(spm_read_vols(SO.img(1).vol)); % meanf range
            SO.img(1).range = [range.min range.max]; % meanf range
            SO.cbar         = [1];
            
            %            [def_flags.bb, def_flags.vox] = bbvox_from_V(spm_vol(underlayImg))
            %getting the bounding box of the underlay image
            bb = spm_get_bbox(spm_vol(underlayImg));
            
            coors=CommonMethods.getCoorRanges(underlayImg);
            
            coormin=coors(1,:);
            coormax=coors(2,:);
            
            
            %             xrange=[bb(1,1), bb(2,1)];
            %             yrange=[bb(1,2), bb(2,2)];
            %             zrange=[bb(2,3), bb(1,3)];
            
            xrange=[coormin(1) coormax(1)];
            yrange=[coormin(2) coormax(2)];
            zrange=[coormin(3) coormax(3)];
            
            dx=1;
            dy=1;
            dx=dx*sign(xrange(2)-xrange(1));
            dy=dy*sign(yrange(2)-yrange(1));
            dz=4*sign(zrange(2)-zrange(1));
            % the slice ranges of z coordinates.
            SO.slices       = floor(zrange(2)):-dz:ceil(zrange(1));
            % the slie range of x and y coordinates
            SO.slicedef=[floor(xrange(1)) dx ceil(xrange(2)); floor(yrange(1)) dy ceil(yrange(2))];
            
            if(~iscell(overlayImgs))
                overlayImgs={overlayImgs};
            end
            
            cmaps={'winter' 'winter(1)'; 'cool' 'cool(1)'; 'spring' 'spring(1)'};
            
            for i=1:length(overlayImgs)
                V=spm_vol(overlayImgs{i});
                SO.img(i+1).vol=V;
                vols=spm_read_vols(V);
                range=CommonMethods.getMatValueRange(vols);
                SO.img(i+1).range=[0.1 1]*range.max;
                
                if length(unique(vols))>2
                    SO.img(i+1).cmap  = eval(cmaps{i, 1});
                else
                    SO.img(2).cmap  = eval(cmaps{i, 2});
                end
            end
            
            SO.figure = spm_figure('GetWin', 'Graphics');
            set(findobj('Tag', 'Graphics'), 'Visible', 'off');
            
            paint(SO);
            [~,name,~]=fileparts(filename);
            title(name);
            print_fig(SO, filename, 'print -dpng');
            close(SO.figure);
        end
        
        function coors=getCoorRanges(fname)
            %this function returns the range of the coordinates, coors, in
            %the image file, fname
            
            V=spm_vol(fname);
            coors=zeros(8,3);
            cont=0;
            for x=1:2
                indx=1;
                if(x==2)
                    indx=V.dim(1);
                end;
                for y=1:2
                    indy=1;
                    if(y==2)
                        indy=V.dim(2);
                    end;
                    for z=1:3
                        indz=1;
                        if(z==2)
                            indz=V.dim(3);
                        end;
                        coor=V.mat*[indx indy indz 1]';
                        cont=cont+1;
                        coors(cont,:)=coor(1:3);
                    end
                end
            end
            coors=[min(coors); max(coors)];
        end
        
        function [terms, delimiters, positions]=decomposeString(str, delimiters0)
            %This function decompose a string according to the delimiters.
            %ters is a cell array of strings contain the substrings of the
            %origional string, str. delimiters is a cell array of
            %delimiters as stored as they appear in the original string.
            
            % eg. str='23+x*y-z*34*abde-56y' delimiters0={'+' '-' '*'} then
            % term = {'23' 'x' 'y' 'z' '34' 'abde' '56y'} and delimiters={'+' '*' '-' '*' '*' '-'}
            % positions=[3 5 7 9 12 17]
            
            len=length(delimiters0);
            positions0=[];
            syms={};
            for i=1:len
                idx=findstr(str, delimiters0{i});
                if(isempty(idx))
                    continue;
                end
                positions0=[positions0 idx];
                syms=horzcat(syms, arrayfun(@(x) delimiters0{i}, idx, 'UniformOutput', false));
            end
            
            if(isempty(syms))
                terms={str};
                delimiters={};
                positions=[];
                return;
            end
            
            [positions, I]=sort(positions0);
            
            if(positions(1)>1)
                terms={str(1:positions(1)-1)};
            else
                terms={};
            end
            
            len1=length(positions);
            
            sym=syms{I(1)};
            delimiters{1}=sym;
            p0=positions(1)+length(sym);
            for i=2:len1
                p=positions(i);
                terms{end+1}=str(p0:p-1);
                sym=syms{I(i)};
                delimiters{i}=sym;
                p0=p+length(sym);
            end
            
            len2=length(str);
            if(p0<len2)
                terms{end+1}=str(p0:len2);
            end
        end
        
        function [sign, factor]=getSignOfMultiplication(signs)
            %%% this function ditermines the sign of the product of the
            %%% multiplication of several terms whose signs are stored in
            %%% the cellarray, signs ('+" and '-' s).
            
            num=length(find(ismember(signs, '-')));
            if(mod(num,2))%%odd number of '-'
                sign='-';
                factor=-1;
            else
                sign='+';
                factor=1;
            end
        end
        
        function displayMatrixes(fname)
            V=spm_vol(fname);
            for i=1:numel(V)
                fprintf('[%g %g %g %g; %g %g %g %g; %g %g %g %g]  %s\n',...
                    V(i).mat(1:3,:)', V(i).fname);
            end
        end
        
        function sleepUntil_time(year, mon, day, hour, min, sec, msg)
            if(~exist('msg','var'))
                msg = ['Sleep until ' num2str(year) '-', num2str(mon) '-', num2str(day) '-', num2str(hour) '-', num2str(min) '-', num2str(sec)];
            end
            dispstat(msg,'init');
            dispstat(msg,'timestamp');
            %dispstat('Begining the process...','keepthis','timespamp');
            dnum0=datenum(year, mon, day, hour, min, sec);
            tvec=clock();
            dnum=datenum(tvec(1), tvec(2), tvec(3), tvec(4), tvec(5), tvec(6));
            while (dnum<dnum0)
                pause(30);
                tvec=clock();
                dnum=datenum(tvec(1), tvec(2), tvec(3), tvec(4), tvec(5), tvec(6));
                dispstat(msg,'timestamp');
            end
            dispstat('Starting work','keepprev','timestamp');
        end
        
        
        function rsyncFolders(folder1, folder2, logfile)
            if(~exist(folder1,'file'))
                error('folder1 does not exist. ---rsyncFolders');
            end
            if(~exist(fileparts(logfile),'file'))
                error('the directory of the logfile does not exist. ---rsyncFolders');
            end
            if(~exist(folder2,'file'))
                mkdir(folder2);
            end
            command=['rsync -auvzr ' folder1 ' ' folder2 ' >> ' logfile];
            system(command);
        end
        
        function makeHTMLFiles(statsdir, cur_outdir, p, getNewScanIDs, getProcessedScanIDs, writeScanIDs, reUpload2bieeglview)
            matfile = fullfile(statsdir, 'SPM.mat');
            
            if ~exist(statsdir, 'dir') || ~exist(matfile, 'file')
                return;
            end
            
            try
                info=load(fullfile(statsdir, 'SPM.mat'));
                new_scanIDs=getNewScanIDs(info.SPM);
            catch err
                return;
            end
            
            if(~exist('reUpload2bieeglview', 'var'))
                reUpload2bieeglview = false;
            end
            
            % Find Subject List text file of subjects already processed
            processed_scanIDs=getProcessedScanIDs(cur_outdir);
            % For each contrast (+/-)
            for conNo = 1:length(info.SPM.xCon)
                p.con = info.SPM.xCon(conNo).name;
                p.conNo = conNo;
                p.conStat = info.SPM.xCon(conNo).STAT;
                if(strcmp(p.conStat,'T'))
                    continue;
                end
                rmpath(genpath('/usr/local/spm8/'));
                addpath(genpath('/usr/local/spm5/'));
                p.spmmatpth = statsdir;
                p.outputbase = fullfile(cur_outdir, p.con);
                
                % For each p value
                for pv = 1:length(p.pvals)
                    
                    % Make blobs
                    if p.doBlobs
                        
                        p.mode = 'blobs';
                        p.p1 = p.pvals(pv);
                        p.outputbase = fullfile(cur_outdir, ['con' num2str(conNo)], p.correction, num2str(p.p1), p.mode);
                        
                        % Make directory
                        if ~isdir(p.outputbase)
                            mkdir(p.outputbase);
                        end
                        
                        % Output text to screen
                        cprintf('blue','%s\n', p.outputbase);
                        
                        % If need to run
                        if length(dir(fullfile(p.outputbase, '*.jpg'))) ~=3 || reUpload2bieeglview || ~isempty(setxor(new_scanIDs, processed_scanIDs(~cellfun(@isempty, processed_scanIDs))))
                            
                            % For each view
                            for v = 1:length(p.views)
                                p.view = p.views{v};
                                try
                                    if(~strcmp(p.conStat,'T'))
                                        %                                     rmpath(genpath('/usr/local/spm8/'));
                                        %                                     addpath(genpath('/usr/local/spm5/'));
                                        makeHTML_spm5('blobs', p);
                                    else
                                        makeHTML('blobs', p);
                                    end
                                    cprintf('green','View: %s\tDONE\n',p.view);
                                    close all
                                catch err
                                    cprintf('red','ERROR\n');
                                    [~, ~] = system(['rm ' fullfile(p.outputbase, '*.jpg')]);
                                end
                            end
                            % If don't need to run
                        else
                            cprintf('blue','%s\n', p.outputbase);
                            cprintf('green','NO NEW SUBJECTS FOUND\n');
                        end
                    end
                    
                    % Make render
                    if p.doRender
                        
                        p.mode = 'render';
                        p.p1 = p.pvals(pv);
                        p.outputbase = fullfile(cur_outdir, ['con' num2str(conNo)], p.correction, num2str(p.p1), p.mode);
                        
                        % Make directory
                        if ~isdir(p.outputbase)
                            mkdir(p.outputbase);
                        end
                        
                        % Output text to screen
                        cprintf('blue','%s\n', p.outputbase);
                        
                        % If need to run
                        if ~exist(fullfile(p.outputbase, 'render.jpg'), 'file') || reUpload2bieeglview || ~isempty(setxor(new_scanIDs, processed_scanIDs(~cellfun(@isempty, processed_scanIDs))))
                            try
                                if(~strcmp(p.conStat,'T'))
                                    %                                     rmpath(genpath('/usr/local/spm8/'));
                                    %                                     addpath(genpath('/usr/local/spm5/'));
                                    makeHTML_spm5('cluster', p);
                                else
                                    makeHTML('cluster', p);
                                end
                                cprintf('green','DONE\n');
                                close all
                            catch err
                                cprintf('red','ERROR\n');
                            end
                            % If don't need to run
                        else
                            cprintf('green','NO NEW SUBJECTS FOUND\n');
                        end
                    end
                    
                    %%% Extract clusters
                    if p.doClusters
                        
                        close all
                        
                        p.mode = 'cluster';
                        p.pThreshold = p.pvals(pv);
                        p.outputbase = fullfile(cur_outdir, ['con' num2str(conNo)], p.correction, num2str(p.pThreshold), p.mode);
                        
                        % Output text to scren
                        cprintf('blue','%s\n', p.outputbase);
                        
                        % If need to run
                        if ~exist(fullfile(p.outputbase, '.complete'), 'file') || reUpload2bieeglview || ~isempty(setxor(new_scanIDs, processed_scanIDs(~cellfun(@isempty, processed_scanIDs))))
                            
                            % Make directory
                            if isdir(p.outputbase)
                                system(['rm -rf ' p.outputbase]);
                            end
                            mkdir(p.outputbase);
                            
                            try
                                if(~strcmp(p.conStat,'T'))
                                    %                                     rmpath(genpath('/usr/local/spm8/'));
                                    %                                     addpath(genpath('/usr/local/spm5/'));
                                    makeHTML_spm5('cluster', p);
                                else
                                    makeHTML('cluster', p);
                                end
                                if(~strcmp(p.conStat,'T'))
                                    %                                     rmpath(genpath('/usr/local/spm5/'));
                                    %                                     addpath(genpath('/usr/local/spm8/'));
                                end
                                %                                makeHTML_FandT('cluster', p);
                                system(['echo >> ' fullfile(p.outputbase, '.complete')]);
                                cprintf('green','Mode: Cluster\nDONE\n');
                                close all
                            catch err
                                cprintf('red','Mode: Cluster\nERROR\n');
                            end
                            % If don't need to run
                        else
                            cprintf('green','NO NEW SUBJECTS FOUND\n');
                        end
                    end
                end
                
                %%% Make threshold plots
                if p.doThreshold
                    
                    p.outputbase = fullfile(cur_outdir, ['con' num2str(conNo)], p.correction, 'threshold');
                    p.mode = 'threshold';
                    threshfile = ['.threshold_' sprintf('%3.3f', p.pvals)];
                    
                    % Make directory
                    if ~isdir(p.outputbase)
                        mkdir(p.outputbase);
                    end
                    
                    % Output text to screen
                    cprintf('blue','%s\n', p.outputbase);
                    
                    % If need to run
                    if ~exist(fullfile(p.outputbase, threshfile), 'file') || length(dir(fullfile(p.outputbase, '*.jpg')))~= 3 || reUpload2bieeglview || ~isempty(setxor(new_scanIDs, processed_scanIDs(~cellfun(@isempty, processed_scanIDs))))
                        
                        [~, ~] = system(['rm ' fullfile(p.outputbase, '.thresh*')]);
                        
                        % For each view
                        for v = 1:length(p.views)
                            p.view = p.views{v};
                            try
                                if(~strcmp(p.conStat,'T'))
                                    %                                     rmpath(genpath('/usr/local/spm8/'));
                                    %                                     addpath(genpath('/usr/local/spm5/'));
                                    makeHTML_spm5('cluster', p);
                                else
                                    makeHTML('cluster', p);
                                end
                                system(['echo >> ' fullfile(p.outputbase, threshfile)]);
                                cprintf('green','View: %s\tDONE\n',p.view);
                            catch err
                                [~, ~] = system(['rm ' fullfile(p.outputbase, '*.jpg')]);
                                cprintf('red','ERROR\n');
                                continue;
                            end
                        end
                        % If don't need to run
                    else
                        cprintf('green','NO NEW SUBJECTS FOUND\n');
                    end
                end
            end
            
            
            %%% Save the workspace to a .mat file
            save((fullfile(cur_outdir, 'bieeglview.mat')));
            writeScanIDs(cur_outdir, new_scanIDs);
        end
        
        function out_imgs = coregister_images(ref_imgs, source_imgs, other_imgs, prefix)
            %
            % Author : Taihao Jin
            %
            % Creation Date : 3/15/2016
            
            % Purpose : Coregisters the source_imgs to the ref_imgs using SPM8. The co-registration
            %           parameters are applied to the other_imgs.
            %
            % Usage : CommonMethods.coregister_images(ref_imgs, source_imgs, other_imgs)
            %
            % Inputs :
            %   anatdir - Input folder containing anatomical scan.
            %   outdir - The output directory to write the co-registered images.
            %   ref_imgs: the reference image.
            %   source_imgs: the source images.
            %   other_imgs: The co-registration parameters are applied to the other_imgs. other_imgs
            %   should be in the same space as the source_imgs.
            
            % Outputs : out_imgs, the images generated in outdir. They are
            % in the space of ref_imgs.
            %
            
            if(~exist('prefix', 'var'))
                prefix='r';
            end
            addpathtj();
            [outdir, name, ext]=fileparts(source_imgs{1});
            outdir_tmp=fullfile(outdir, 'tmp_coreg');
            mkdir(outdir_tmp);
            %do the coregistration in a temp folder
            for i=1:length(source_imgs)
                [path, name, ext]=fileparts(source_imgs{i});
                system(['cp ' source_imgs{i} ' ' outdir_tmp]);
                source_imgs{i}=fullfile(outdir_tmp, [name ext]);
            end
            
            spm('defaults', 'fmri')
            spm_jobman('initcfg');
            
            load /home/harshad/Scripts/jobTemplates/Prodrome_coregister_SPM8.mat
            
            matlabbatch{1}.spm.spatial.coreg.estwrite.ref = ref_imgs;
            matlabbatch{1}.spm.spatial.coreg.estwrite.source = source_imgs;
            matlabbatch{1}.spm.spatial.coreg.estwrite.other = other_imgs;
            matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = prefix;
            
            spm_jobman('run_nogui', matlabbatch);
            
            %getPrefixedFiles(fnames, prefix)
            
            imgs=horzcat(source_imgs, other_imgs);
            out_imgs=cell(1,length(source_imgs));
            for i=1:length(source_imgs)
                [path, name, ext]=fileparts(imgs{i});
                img=[prefix name ext];
                out_img=fullfile(outdir, img);
                system(['cp ' fullfile(outdir_tmp, img) ' ' outdir]);
                out_imgs{i}=out_img;
            end
            out_imgs=horzcat(out_imgs, CommonMethods.getPrefixedFiles(other_imgs, prefix));
            system(['rm ' fullfile(outdir_tmp, '*.*')]);
            rmdir(outdir_tmp);
        end
        
        function [filepgroups, uniqueFiles] = getDuplicatedPresetationFiles(fnames)
            len=length(fnames);
            filepgroups={};
            uniqueFiles={};
            if(length(fnames)<2)
                fnames=fnames;
            end
            while(len >1)
                gps=1:len;
                fname1=fnames{1};
                uniqueFiles{end+1}=fname1;
                for i=2:len
                    fname2=fnames{i};
                    if(CommonMethods.AreDuplicatedPresentationFiles(fname1, fname2))
                        gps(i)=1;
                    end
                end
                idx=find(gps==1);
                if(length(idx)>1)
                    filepgroups{end+1}=fnames(idx);
                end
                fnames=fnames(gps~=1);
                len=length(fnames);
            end
            if(~isempty(fnames))
                uniqueFiles{end+1}=fnames{1};
            end
        end
        
        function dn = getPresentationFileDatenum(infile)
            [~, sceDate, sceTime] = read_presentation_log_hdr(infile);
            nums=[CommonMethods.str2nums_decimal(sceDate) CommonMethods.str2nums_decimal(sceTime)];
            dn=datenum(nums);
        end
        
        
        function IS = IdenticalCellarray(ca1, ca2)
            %assumption: 1D cellarray, and the elements are numbers or
            %strings.
            IS=true;
            len1=length(ca1);
            len2=length(ca2);
            if(len1~=len2)
                IS=false;
                return;
            end
            
            idx1=find(cellfun(@(x) isempty(x), ca1));
            idx2=find(cellfun(@(x) isempty(x), ca2));
            if(any(idx1~=idx2))
                IS=false;
                return;
            end
            ca1(idx1)={0};
            ca2(idx1)={0};
            
            idxn1=find(cellfun(@(x) isa(x,'double'), ca1));
            idxn2=find(cellfun(@(x) isa(x,'double'), ca2));
            idxs1=find(cellfun(@(x) ischar(x), ca1));
            idxs2=find(cellfun(@(x) ischar(x), ca2));
            if(any(idxn1~=idxn2)||any(idxs1~=idxs2))
                IS=false;
                return;
            end
            if((length(idxn1)+length(idxn2)+length(idxs1)+length(idxs2))~=len1+len2)
                error('ca1 or/and ca2 are not cellarray of string and number in CommmonMethods.IdenticalCellarray(ca1, ca2)');
            end
            
            if(any(arrayfun(@(x) ~strcmp(ca1{x},ca2{x}), idxs1)))
                IS=false;
                return;
            end
            
            if(isempty(idxn1))
                return;
            end
            
            ca1=ca1(idxn1);
            ca2=ca2(idxn1);
            idx1=find(cellfun(@(x) isnan(x), ca1));
            idx2=find(cellfun(@(x) isnan(x), ca2));
            if(length(idx1)~=length(idx2))
                IS=false;
                return;
            end
            
            if(any(idx1~=idx2))
                IS=false;
                return;
            end
            ca1(idx1)={0};
            ca2(idx1)={0};
            
            %now ca1 and ca2 are numeric cellarry
            if(any(arrayfun(@(x) ca1{x}~=ca2{x}, 1:length(ca1))))
                IS=false;
                return;
            end
        end
        
        function IS = AreDuplicatedPresentationFiles(infile1, infile2)
            IS=true;
            datmat1 = read_presentation_log(infile1);
            datmat2 = read_presentation_log(infile2);
            [~, sceDate1, sceTime1] = read_presentation_log_hdr(infile1);
            [~, sceDate2, sceTime2] = read_presentation_log_hdr(infile2);
            if(~strcmp(sceDate1, sceDate2)||~strcmp(sceTime1, sceTime2))
                IS=false;
                return;
            end
            dim1=size(datmat1);
            dim2=size(datmat2);
            if(dim1(1)~=dim2(1)||dim1(1)~=dim2(1))
                IS=false;
                return;
            end
            for c=1:dim1(2)
                ttime1=datmat1(2:dim1(1),c);
                ttime2=datmat2(2:dim2(1),c);
                if(~CommonMethods.IdenticalCellarray(ttime1,ttime2))
                    IS=false;
                    return
                end
            end
        end
        function IS = IsCellarrayOfString(ca)
            IS=true;
            if(~iscell(ca))
                IS=false;
                return;
            end
            if(~all(cellfun(@(x) ischar(x),ca)))
                IS=false;
                return;
            end
        end
        function IS = IsQuitedPresentationFile(infile)
            IS=false;
            datmat = read_presentation_log(infile);
            fields=datmat(:,3);
            fields(cellfun(@(x)isempty(x),fields))={'empty'};
            if(CommonMethods.IsCellarrayOfString(fields))
                if(any(ismember(fields, 'Quit')))
                    IS=true;
                end
            else
                return;
            end
        end
        
        function fnames = sortPresentationFiles(infiles)
            dns=cellfun(@(x)CommonMethods.getPresentationFileDatenum(x), infiles);
            [tmp, order]=sort(dns);
            fnames=infiles(order);
        end
        
        function m=rangeAddition(m1, m2)
            m=m1+m2;
            idx=find(m1<0.*m2>0);
            m(idx)=m1(idx)-m2(idx);
        end
        
        function nifti3to4d(folder)
            %to concatenate all 3D nii images in the directory foler into a 4D image
            lst=dir(folder);
            fnames={};
            newName=fullfile(folder,'f.nii');
            for n=1:length(lst)
                fname=lst(n).name;
                [path,name,ext]=fileparts(fname);
                if(isempty(ext))
                    continue;
                end
                if(~strcmpi(ext,'.nii'))
                    continue;
                end
                if(strcmpi(fname,'f.nii'))
                    continue;
                end
                fnames{end+1}=fullfile(folder,fname);
            end
            fnames=sort(fnames);
            len=length(fnames);
            if(len<1)
                return;
            elseif(len==1)
                copyfile(fnames{1}, newName);
                return;
            end
            dt=0;%to keep the same data type
            try
                V4=spm_file_merge(fnames,newName,dt);
            catch err
                cprintf('red', '%s\n', err.message);
            end
        end
        
        function testDisp()
            close(gcf);
            xTickLabels={'hc1' 'hc2' 'fs3' 'sfsfafs' ;'ct1' 'ct2' 'safs' 'fsasfsa'};
            starts={'*' '***' '**' '*'};
            means=5*rand(2,4)+5;
            sems=rand(2,4);
            title='Structures having significant group difference in volume change';
            
            CommonMethods.dispMeanSemComparison('xTickLabels', xTickLabels, 'means', means, 'sems', sems, 'title', title, 'stars', starts);
        end
        
        function dispMeanSemComparison(varargin)
            len=length(varargin);
            if(mod(len,2)~=0)%the input areguments has to be a series of name and value pairs.
                return;
            end
            l2=len/2;
            vnames=arrayfun(@(x) varargin{2*x+1}, 0:(l2-1), 'UniformOutput', false);
            values=arrayfun(@(x) varargin{2*x+2}, 0:(l2-1), 'UniformOutput', false);
            
            lvnames=cellfun(@(x)lower(x),vnames,'UniformOutput',false);
            
            inds=find(ismember(lvnames,'xticklabels'));
            if(length(inds)~=1)
                return;
            end
            ind=inds(1);
            xTickLabels0=values{ind};
            
            
            inds=find(ismember(lvnames,'means'));
            if(length(inds)~=1)
                return;
            end
            ind=inds(1);
            means=values{ind};
            
            inds=find(ismember(lvnames,'title'));
            if(length(inds)~=1)
                titleT={};
            else
                ind=inds(1);
                titleT=values{ind};
            end
            
            inds=find(ismember(lvnames,'stars'));
            if(length(inds)==1)
                ind=inds(1);
                stars=values{ind};
            else
                stars={};
            end
            
            inds=find(ismember(lvnames,'fname'));
            if(length(inds)==1)
                ind=inds(1);
                fname=values{ind};
            else
                fname={};
            end
            
            inds=find(ismember(lvnames,'ylabel'));
            if(length(inds)==1)
                ind=inds(1);
                yLabel=values{ind};
            else
                yLabel={};
            end
            
            inds=find(ismember(lvnames,'sems'));
            if(length(inds)~=1)
                return;
            end
            ind=inds(1);
            sems=values{ind};
            
            %             inds=find(ismember(lvnames,'pvalues'));
            %             if(length(inds)~=1)
            %                 return;
            %             end
            %             ind=inds(1);
            %             pVs=vnames{ind};
            
            inds=find(ismember(lvnames,'hAxis'));
            if(length(inds)~=1)
                hAxis=gca();
            else
                ind=inds(1);
                hAxis=values{ind};
            end
            
            if(isempty(hAxis))
                hAxis=gca;
            end
            
            axes(hAxis);
            
            rows=size(means,1);
            cols=size(means,2);
            
            x=0;
            
            xTickLabels={};
            xTick=[];
            fc=[0 0 1];
            
            m=CommonMethods.rangeAddition(means,sems);
            yn=min(min(m));
            yx=max(max(m));
            ex=max(max(sems));
            en=min(min(sems));
            if(yn>0)
                yn=0;
            end
            if(yx<0)
                yx=0;
            end
            margin=0.05*(yx-yn);
            
            e=0;
            for c=1:cols
                %               fc=rand(1,3);
                xc=0;
                yxt=yn;
                if(c>1)
                    x=x+1;
                end
                for r=1:rows
                    x=x+1;
                    xc=xc+x;
                    %                    xTick(end+1)=x-0.5;
                    xTick(end+1)=x;
                    y=means(r,c);
                    e=sems(r,c);
                    if(y+e>yxt)
                        yxt=y+e;
                        if(yxt<0)
                            yxt=0;
                        end
                    end
                    label=CommonMethods.makeFigureLabel(xTickLabels0{r,c});
                    hbar=bar(x,y);
                    set(hbar,'facecolor',fc);
                    %                    fc=get(hbar,'facecolor');
                    %                     if(~isempty(stars))
                    %                         text(x,y+1.5*e,stars{c});
                    %                     end
                    hold(hAxis,'on');
                    if(e~=0)
                        
                        nege=zeros(size(e));
                        pose=zeros(size(e));
                        nege(y<=0)=e(y<=0);
                        pose(y>0)=e(y>0);
                        herr=errorbar(x,y,nege,pose);
                        
                        %                         herr=errorbar(x,y,e);
                        set(herr,'LineStyle','none');
                        set(herr,'LineWidth',3);
                        set(herr,'color',fc);
                        hold(hAxis,'off');
                        set(hAxis,'fontsize',14);
                        hold(hAxis,'on');
                    end
                    xTickLabels{end+1}=label;
                end
                xc=xc/rows;
                yxt=yxt+1.3*margin;
                shiftx=[0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
                if(~isempty(stars))
                    th=text(xc,yxt,CommonMethods.makeFigureLabel(stars{c}));
                    set(th,'HorizontalAlignment', 'center');
                end
            end
            
            xlim([0 x+1]);
            
            yr=[yn-margin yx+2.5*margin];
            if(yr(1)>0)
                yr(1)=0;
            end
            ylim(yr);
            
            set(hAxis,'xTick',xTick);
            set(hAxis,'xTickLabels',xTickLabels);
            xticklabel_rotate([],90);
            if(~isempty(titleT))
                t=title(titleT);
                set(t, 'fontsize', 14);
            end
            if(~isempty(yLabel))
                yl=ylabel(hAxis,yLabel);
                set(yl,'fontsize', 14);
            end
            fig=gcf();
            %           set(hAxis,'box','off')
            %           set(hAxis,'ticklength',[0 0])
            if(~isempty(fname))
                fig.PaperPositionMode = 'auto';
                print('-bestfit',fname,'-dpdf');
            end
        end
        
        function dispOneSampleTTest(data, names, title, ylabel, fname)
            figure();
            means=cellfun(@(x)mean(x), data);
            sems=cellfun(@(x)std(x)/sqrt(length(x)),data);
            len=length(names);
            stars=cell(1,len);
            for i=1:len
                [~,p]=ttest(data{i});
                if(p > 0.05)
                    stars{i}='p > 0.05';
                elseif(p >= 0.01)
                    stars{i}='p < 0.05';
                elseif(p > 0.001)
                    stars{i}='p < 0.01';
                elseif(p > 0.0001)
                    stars{i}= 'p < 0.001';
                else
                    stars{i}= 'p < 0.0001';
                end
            end
            CommonMethods.dispMeanSem('hAxis', gca(),'xTickLabels', names, 'means', means, 'sems', sems, 'title', title, 'stars', stars, 'yLabel', ylabel,'fname',fname);
        end
        
        function dispMeanSem(varargin)
            len=length(varargin);
            if(mod(len,2)~=0)%the input areguments has to be a series of name and value pairs.
                return;
            end
            l2=len/2;
            vnames=arrayfun(@(x) varargin{2*x+1}, 0:(l2-1), 'UniformOutput', false);
            values=arrayfun(@(x) varargin{2*x+2}, 0:(l2-1), 'UniformOutput', false);
            
            lvnames=cellfun(@(x)lower(x),vnames,'UniformOutput',false);
            
            inds=find(ismember(lvnames,'xticklabels'));
            if(length(inds)~=1)
                return;
            end
            ind=inds(1);
            xTickLabels0=values{ind};
            
            
            inds=find(ismember(lvnames,'means'));
            if(length(inds)~=1)
                return;
            end
            ind=inds(1);
            means=values{ind};
            
            inds=find(ismember(lvnames,'title'));
            if(length(inds)~=1)
                titleT={};
            else
                ind=inds(1);
                titleT=values{ind};
            end
            
            inds=find(ismember(lvnames,'stars'));
            if(length(inds)==1)
                ind=inds(1);
                stars=values{ind};
            else
                stars={};
            end
            
            inds=find(ismember(lvnames,'fname'));
            if(length(inds)==1)
                ind=inds(1);
                fname=values{ind};
            else
                fname={};
            end
            
            inds=find(ismember(lvnames,'ylabel'));
            if(length(inds)==1)
                ind=inds(1);
                yLabel=values{ind};
            else
                yLabel={};
            end
            
            inds=find(ismember(lvnames,'sems'));
            if(length(inds)~=1)
                return;
            end
            ind=inds(1);
            sems=values{ind};
            
            %             inds=find(ismember(lvnames,'pvalues'));
            %             if(length(inds)~=1)
            %                 return;
            %             end
            %             ind=inds(1);
            %             pVs=vnames{ind};
            
            inds=find(ismember(lvnames,'hAxis'));
            if(length(inds)~=1)
                hAxis=gca();
            else
                ind=inds(1);
                hAxis=values{ind};
            end
            
            if(isempty(hAxis))
                hAxis=gca;
            end
            
            axes(hAxis);
            
            rows=size(means,1);
            cols=size(means,2);
            
            x=0;
            
            xTickLabels={};
            xTick=[];
            fc=[0 0 1];
            
            m=CommonMethods.rangeAddition(means,sems);
            yn=min(min(m));
            yx=max(max(m));
            ex=max(max(sems));
            en=min(min(sems));
            if(yn>0)
                yn=0;
            end
            if(yx<0)
                yx=0;
            end
            margin=0.05*(yx-yn);
            
            e=0;
            for c=1:cols
                %               fc=rand(1,3);
                xc=0;
                yxt=yn;
                if(c>1)
                    x=x+1;
                end
                for r=1:rows
                    x=x+1;
                    xc=xc+x;
                    %                    xTick(end+1)=x-0.5;
                    xTick(end+1)=x;
                    y=means(r,c);
                    e=sems(r,c);
                    if(y+e>yxt)
                        yxt=y+e;
                        if(yxt<0)
                            yxt=0;
                        end
                    end
                    label=CommonMethods.makeFigureLabel(xTickLabels0{r,c});
                    hbar=bar(x,y);
                    set(hbar,'facecolor',fc);
                    %                    fc=get(hbar,'facecolor');
                    %                     if(~isempty(stars))
                    %                         text(x,y+1.5*e,stars{c});
                    %                     end
                    hold(hAxis,'on');
                    if(e~=0)
                        nege=zeros(size(e));
                        pose=zeros(size(e));
                        nege(y<=0)=e(y<=0);
                        pose(y>0)=e(y>0);
                        herr=errorbar(x,y,nege,pose);
                        set(herr,'LineStyle','none');
                        set(herr,'LineWidth',3);
                        set(herr,'color',fc);
                        hold(hAxis,'off');
                        set(hAxis,'fontsize',14);
                        hold(hAxis,'on');
                    end
                    xTickLabels{end+1}=label;
                end
                xc=xc/rows;
                yxt=yxt+1.3*margin;
                shiftx=[0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
                if(~isempty(stars))
                    th=text(xc,yxt,CommonMethods.makeFigureLabel(stars{c}));
                    set(th,'HorizontalAlignment', 'center');
                end
            end
            
            xlim([0 x+1]);
            
            yr=[yn-margin yx+2.5*margin];
            if(yr(1)>0)
                yr(1)=0;
            end
            ylim(yr);
            
            set(hAxis,'xTick',xTick);
            set(hAxis,'xTickLabels',xTickLabels);
            xticklabel_rotate([],90);
            if(~isempty(titleT))
                t=title(titleT);
                set(t, 'fontsize', 14);
            end
            if(~isempty(yLabel))
                yl=ylabel(hAxis,yLabel);
                set(yl,'fontsize', 14);
            end
            fig=gcf();
            %           set(hAxis,'box','off')
            %           set(hAxis,'ticklength',[0 0])
            if(~isempty(fname))
                fig.PaperPositionMode = 'auto';
                print('-bestfit',fname,'-dpdf');
            end
        end
        
        function shiftUIControls(hObject, dx, dy)
            objs=findall(hObject);
            uics={};
            for i=1:length(objs)
                obj=objs(i);
                if(strcmp(get(obj,'Type'),'uicontrol'))
                    uics{end+1}=obj;
                    pos=get(obj,'Position');
                end
            end
            
            for i=1:length(uics)
                uic=uics{i};
                pos=get(uic,'Position');
                pos(2)=pos(2)+dy;
                pos(1)=pos(1)+dx;
                set(uic,'Position', pos);
            end
        end
        
        function selectUnitControl(hObject, dx, dy)
            objs=findall(hObject);
            uics={};
            for i=1:length(objs)
                obj=objs(i);
                if(strcmp(get(obj,'Type'),'uicontrol'))
                    uics{end+1}=obj;
                    pos=get(obj,'Position');
                end
            end
            
            for i=1:length(uics)
                uic=uics{i};
                pos=get(uic,'Position');
                pos(2)=pos(2)+dy;
                pos(1)=pos(1)+dx;
                set(uic,'Position', pos);
            end
        end
        
        function status=subtractImgages(imgb,imga,changeImg)
            %save the subtraction image (the first subtract the second) as changeImg
            status=0;
            V=spm_vol(imgb);
            Vc=V;
            Vc.dt(1)=16;
            Vc.fname=changeImg;
            try
                spm_write_vol(Vc, spm_read_vols(V)-spm_read_vols(spm_vol(imga)));
            catch err
                return;
            end
            status=1;
        end
        
        function status=divideImgages(imgb,imga,Imgq)
            %save the image division(the first divide the second) as changeImg
            status=0;
            V=spm_vol(imgb);
            Vq=V;
            Vq.dt(1)=16;
            Vq.fname=Imgq;
            try
                spm_write_vol(Vq, spm_read_vols(V)./spm_read_vols(spm_vol(imga)));
            catch err
                return;
            end
            status=1;
        end
        
        function updateListBoxValue(hObject)
            objs=findall(hObject);
            for i=1:length(objs)
                obj=objs(i);
                try
                    style=get(obj,'Style');
                catch err
                    continue;
                end
                if(strcmp(style,'listbox'))
                    len=length(get(obj,'String'));
                    val=get(obj,'Value');
                    if(val>len)
                        set(obj,'Value',len);
                    end
                end
            end
        end
        function IS=IsOverlappingRanges(range1, range2)
            IS=range1(1) > range2(2) || range1(2) < range2(1);
        end
        
        function name=getDefaultColumnNames(cols)
            str='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
            label='0123456789ABCDEFGHIJKLMNOP';
            
            numc=dec2base(cols-1,26);
            len=length(numc);
            name=char(len);
            for i=1:len
                name(i)=str(strfind(label,numc(i)));
            end
        end
        
        function header=getDefaultTableHeader(width)
            header=arrayfun(@(x)CommonMethods.getDefaultColumnNames(x), 1:width,'UniformOutput',false);
        end
        
        function [m,keep]=removeDuplicatedRows(m0)
            rows=size(m0,1);
            m1=m0-min(min(m0))+1;
            mx=max(max(m1))+1;
            m2=mx*mx*(m1(:,1))+mx*(m1(:,2))+m1(:,3);
            
            keep=ones(1,rows);
            for r=1:rows-1
                if(~keep(r))
                    continue;
                end
                tmp=m2(r+1:rows)-m2(r);
                idx=find(tmp==0);
                for t=1:length(idx)
                    keep(idx(t)+r)=0;
                end
            end
            m=m0(find(keep),:);
        end
        
        function X=refineDatabase(X)
            len=length(X.DB);
            DB_sorted=cell(1,len);
            DB_eroded=cell(1,len);
            DB_sorted_eroded=cell(1,len);
            se=ones(3,3,3);%the structurl element for the eroion
            mat=X.mat;
            for i=1:len
                db=X.DB{i};
                anatomy=db.anatomy;
                anatomy(cellfun(@(x)isempty(x),anatomy))={'unspecified'};
                %%creating a data base with sorted anatomy names
                [anatomy_sorted, inds]=sort(anatomy);
                dbs.anatomy=anatomy_sorted;
                dbs.type=db.type;
                len1=length(anatomy);
                dim=size(db.mnilist);
                
                sortedlist=zeros(dim);
                for id=1:len1
                    ids=inds(id);
                    tmp=zeros(dim);
                    tmp(db.mnilist==ids)=id;
                    sortedlist=sortedlist+tmp;
                end
                dbs.mnilist=sortedlist;
                DB_sorted{i}=dbs;
                
                %%creating a data base with erroded rois
                
                liste=zeros(dim);
                for id=1:len1
                    tmp=zeros(dim);
                    tmp(db.mnilist==id)=id;
                    tmp=imerode(tmp,se);
                    liste=liste+tmp;
                end
                dbe.anatomy=db.anatomy;
                dbe.type=db.type;
                dbe.mnilist=liste;
                DB_eroded{i}=dbe;
                
                %%creating a data base with sorted anatomy names and
                %%erroded rois
                listse=zeros(dim);
                for id=1:len1
                    ids=inds(id);
                    tmp=zeros(dim);
                    tmp(db.mnilist==ids)=id;
                    tmp=imerode(tmp,se);
                    listse=listse+tmp;
                    name=anatomy_sorted{id};
                    name=strrep(name,' ','_');
                    name=strrep(name,'-','_');
                    [x y z]=ind2sub(dim,find(tmp));
                    coor=[x y z];
                    mnis=CommonMethods.ind2coor_mat(coor,mat);
                    cmd=['wholeMaskMNIAll_eroded.' name ' =mnis;'];
                    eval(cmd);
                end
                
                dbse.anatomy=dbs.anatomy;
                dbse.type=db.type;
                dbse.mnilist=listse;
                DB_sorted_eroded{i}=dbse;
            end
            X.DB_sorted=DB_sorted;
            X.DB_eroded=DB_eroded;
            X.DB_sorted_eroded=DB_sorted_eroded;
            X.wholeMaskMNIAll_eroded=wholeMaskMNIAll_eroded;
        end
        
        function strca=num2str_cellarray(ca)
            dim=size(ca);
            strca=ca;
            idx=find(cellfun(@(x)isnumeric(x),ca));
            for i=1:length(idx)
                ind=idx(i);
                strca{ind}=num2str(ca{ind});
            end
        end
        function q=getVariableCorrectedQuantity(q,var0,include)
            %this function do the var control for q.table
            %var0 should be a collum vector
            if(~exist('include','var'))
                include=true(length(var0),1);
            end
            
            table0=q.table;
            names=q.colNames;
            dim=size(table0);
            table=nan(dim);
            figure;
            num=0;
            axes=cell(1,2);
            for c=1:dim(2)
                y0=table0(:,c);
                idx=find(~isnan(y0).*include);
                y=y0(idx);
                var=var0(idx);
                
                idx1=1:2:length(var);
                y1=y(idx1);
                var1=var(idx1);
                [r,m,b]=regression(var1',y1');
                [r1,p]=corr(var1,y1);
                if(p>0.05)
                    %                    continue;
                end
                
                num=num+1;
                if(num>16)
                    figure;
                    num=1;
                end
                subplot(4,4,num);
                mina=min(var);
                maxa=max(var);
                %                scatter(age1, y1);
                scatter(var0, y0);
                y=y0-m*(var0-mina);
                hold(gca);
                line([mina;maxa],[b+m*mina;b+m*maxa]);
                num=num+1;
                axes{1}=gca;
                xlim([55 85]);
                xlabel('age');
                ylabel('score');
                title(CommonMethods.makeFigureLabel(names{c}));
                
                subplot(4,4,num);
                scatter(var0, y);
                xlim([55 85]);
                xlabel('age');
                ylabel('score');
                title(CommonMethods.makeFigureLabel([names{c} '_v. c.']));
                axes{2}=gca;
                CommonMethods.setEqualYScales(axes);
                table(:,c)=y;
            end
            q.table_varc=table;
        end
        function v=calVoxelVolume(fname)
            V=spm_vol(fname);
            mat=V.mat;
            v=abs(det(mat(1:3,1:3)));
            %              c1=[0 0 0 1]';
            %              c2=[1 0 0 1]';
            %              c3=[0 1 0 1]';
            %              c4=[0 0 1 1]';
            %              coor1=mat*c1;
            %              coor1=coor1(1:3);
            %              coor2=mat*c2;
            %              coor2=coor2(1:3);
            %              coor3=mat*c3;
            %              coor3=coor3(1:3);
            %              coor4=mat*c4;
            %              coor4=coor4(1:3);
            %              l1=dist(coor1',coor2);
            %              l2=dist(coor1',coor3);
            %              l3=dist(coor1',coor4);
            %              v=l1*l2*l3;
        end
        
        function autoAdjustColumnWidths(fname)
            e=actxserver('excel.application');
            eW=e.Workbooks;
            eF=eW.Open(fname); % your filename here
            eS=eF.ActiveSheet;
            cols=eS.UsedRange.Columns.Count;
            rows=eS.UsedRange.Rows.Count;
            r=eS.UsedRange.Address;
            range=eS.Range(r);
            range.EntireColumn.AutoFit;
            eF.Save;
            eF.Close; % close the file
            e.Quit; % close Excel entirely
        end
        
        function idim=getImageDimension(fname);
            %this function returns a structure containing the image
            %dimension related information.
            %dim: dimension
            %mat: the matrix for volxel index to coordinate transition
            %lwd: the voxel size in mm along the firt, second and the
            %third direction.
            %voxelV: the voxel size in mm^3.
            %vol: the size of the entire image in mm^3
            %size: the size of the image (length, width and the height) in
            %mm
            V=spm_vol(fname);
            mat=V.mat;
            dim=V.dim;
            dim(4)=length(V);
            v=abs(det(mat(1:3,1:3)));
            c1=[0 0 0 1]';
            c2=[1 0 0 1]';
            c3=[0 1 0 1]';
            c4=[0 0 1 1]';
            coor1=mat*c1;
            coor1=coor1(1:3);
            coor2=mat*c2;
            coor2=coor2(1:3);
            coor3=mat*c3;
            coor3=coor3(1:3);
            coor4=mat*c4;
            coor4=coor4(1:3);
            l1=CommonMethods.dist(coor1,coor2);
            l2=CommonMethods.dist(coor1,coor3);
            l3=CommonMethods.dist(coor1,coor4);
            idim.dim=dim;
            idim.lwd=[l1 l2 l3];
            idim.voxelV=v;
            idim.size=[l1*dim(1) l2*dim(2) l3*dim(3)];
            idim.Vol=dim(1)*dim(2)*dim(3)*v;
            idim.mat=V.mat;
        end
        
        function shape=getShape(shape)
            %this function returns the voxel indexes included in the shape
            %specified by the structure shape.
            %the shape should have two input parameters (fields): fname and
            %size
            shape.err=[];
            if(~isfield(shape,'idim'))
                if(~isfield(shape,'fname'))
                    shape.err='The structure shape should contain at least one of the two fields "idim" and "fname"';
                    return;
                else
                    shape.idim=CommonMethods.getImageDimension(shape.fname);
                end
            end
            if(~isfield(shape,'name'))
                shape.err='The name of the shape has not been specified';
                return;
            end
            
            lwd=shape.idim.lwd;
            radius=shape.size;%woud be half of the edge length of the squre of the radius of the sphere
            %shape.size is in mm
            sizeInVoxels=[ceil(radius/lwd(1)) ceil(radius/lwd(2)) ceil(radius/lwd(3))];
            tmp=numel(zeros(2*sizeInVoxels+1));
            vIndexes=zeros(tmp,3);
            dist2=zeros(tmp,1);%the distance from the center to each voxel
            radius2=radius*radius;
            num=0;
            for l1=-sizeInVoxels(1):sizeInVoxels(1)
                for l2=-sizeInVoxels(2):sizeInVoxels(2)
                    for l3=-sizeInVoxels(3):sizeInVoxels(3)
                        num=num+1;
                        vIndexes(num,1)=l1;
                        vIndexes(num,2)=l2;
                        vIndexes(num,3)=l3;
                        d=[l1,l2,l3].*lwd;
                        d=sum(d.*d);
                        dist2(num)=d;
                    end
                end
            end
            switch shape.name
                case 'sphere'
                    shape.vIndexes=vIndexes(dist2<=radius2,:);
                case 'square'
                    shape.vIndexes=vIndexes;
                otherwise
                    shape.err='The name of the shape must be "spere" or "square"';
                    return;
            end
        end
        
        function shape=updateShapePosition(shape,center)
            %cetner: in voxel index
            %curVIndexes: the voxel indexes of the shape centered at the
            %"center"
            %validIndexe: indicating whether any of curVIndexes are
            %outsize of the dimension of the image (shape.idim.dim).
            %see CommonMethods.getShape(shape) for other fields of shape.
            
            %          if(~isfield(shape, 'curVIndexes'))
            dim=size(shape.vIndexes);
            shape.curVIndexes=zeros(dim);
            %          end
            validIndexes=1;
            
            for i=1:3
                inds=shape.vIndexes(:,i)+center(i);
                shape.curVIndexes(:,i)=inds;
                validIndexes=validIndexes.*(inds>=1).*(inds<=shape.idim.dim(i));
            end
            shape.curVIndexes=shape.curVIndexes(validIndexes>0,:);%a mx3 matrix
            shape.curVIndexes1=sub2ind(shape.idim.dim(1:3),shape.curVIndexes(:,1),shape.curVIndexes(:,2),shape.curVIndexes(:,3));
        end
        
        function fm=getRandomFrequencyMap(shape,num,mask)
            %this function generates a frequency map after placing "num"
            %"shapes" with uniform probability within "maskInds";
            dim=shape.idim.dim;
            if ~exist('mask','var')
                mask=ones(dim);
            end
            maskInds=find(mask);
            len=length(maskInds);
            inds=randi([1,len],1,num);
            centerInds=maskInds(inds);
            [x,y,z]=ind2sub(dim,centerInds);
            centers=[x y z];
            fm=zeros(dim);
            volsum=zeros(1,num);
            indexRatio=zeros(1,num);
            cs=zeros(1,num);
            for i=1:num
                center=centers(i,:);
                shape=CommonMethods.updateShapePosition(shape,center);
                a=fm(shape.curVIndexes1);
                fm(shape.curVIndexes1)=fm(shape.curVIndexes1)+1;
                b=fm(shape.curVIndexes1);
                c=sum(b-a)-length(shape.curVIndexes1);
                if i>1100
                    a=1;
                end
                volsum(i)=sum(sum(sum(fm)))/length(shape.vIndexes);
                indexRatio(i)=length(shape.curVIndexes1)/length(shape.vIndexes);
                cs(i)=c;
            end
            fm=fm/num;
        end
        function fm=getFrequencyMap(shape,centers)
            %this function generates a frequency map after placing "num"
            %"shapes" with uniform probability within "maskInds";
            %coors: a N by 3 matrix, each row is a 3D coordinates
            dim=shape.idim.dim;
            fm=zeros(dim);
            num=size(centers,1);
            for i=1:num
                center=centers(i,:);
                shape=CommonMethods.updateShapePosition(shape,center);
                fm(shape.curVIndexes1)=fm(shape.curVIndexes1)+1;
            end
            fm=fm/num;
        end
        
        function makeFileList(folder,fname,exts)
            %this function write the list of files in folder into the file
            %fname.
            %folder: the directory
            %fname: into which the file list will be written to
            
            lst=dir(folder);
            if(~exist('exts','var'))
                exts=[];
            end
            len=length(lst);
            fid=fopen(fullfile(folder,fname),'wt');
            for i=1:len
                tmp=lst(i).name;
                [~,name,ext]=fileparts(tmp);
                if(length(name)<5)
                    continue;
                end
                if(~isempty(exts))
                    if(~any(ismember(exts,ext)))
                        continue;
                    end
                end
                fprintf(fid,'%s\n',fullfile(folder,tmp));
            end
            fclose(fid);
        end
        
        function fileList=getFilesInFolder(folder,exts)
            %this function returns all files in folder, filtering by the
            %extension in "exts";
            if(~exist('exts','var'))
                exts=[];
            end
            if(ischar(exts))
                exts={exts};
            end
            lst=dir(folder);
            len=length(lst);
            fileList={};
            for i=1:len
                tmp=lst(i).name;
                [~,name,ext]=fileparts(tmp);
                if(length(name)<5)
                    continue;
                end
                if(~isempty(exts))
                    if(~any(ismember(exts,ext)))
                        continue;
                    end
                end
                fileList{end+1}=fullfile(folder,tmp);
            end
        end
        
        function fileList=getSubfolderList(folder)
            %returns all directories under folder
            lst=dir(folder);
            len=length(lst);
            fileList={};
            for i=1:len
                tmp=lst(i);
                if(~tmp.isdir)
                    continue;
                end
                name=tmp.name;
                if(strcmp(name(1),'.')||strcmp(name(1),'..'))
                    continue;
                end
                fileList{end+1}=fullfile(folder,name);
            end
        end
        
        function fileList=getSubfolderList_recursive(root)
            subfolders=CommonMethods.getSubfolderList(root);
            fileList=subfolders;
            for s=1:length(subfolders)
                sfs=CommonMethods.getSubfolderList_recursive(subfolders{s});
                fileList=horzcat(fileList, sfs);
            end
        end
        function st=getMeanSemStr(a,precision1, precision2)
            len=length(a);
            st=CommonMethods.meanSem2str(mean(a),std(a)/sqrt(len-1),precision1, precision2);
        end
        function str=meanSem2str(mean,sem,precision1, precision2)
            if(~exist('precision1','var'))
                precision1=2;
            end
            if(~exist('precision2','var'))
                precision2=3;
            end
            
            str=[CommonMethods.num2str_decimal(mean,precision1) ' (' CommonMethods.num2str_decimal(sem,precision2) ')'];
        end
        function [n,m]=getNumberOfDigits(a)
            n=0;
            tol=0.00001;
            a=abs(a);
            a0=fix(a);
            b=a-a0;
            while(a0>0)
                n=n+1;
                a0=fix(a0/10);
            end
            m=0;
            while(b>0.00001)
                m=m+1;
                b=b*10;
                b=b-fix(b+tol);
            end
        end
        function [n,m]=getDigits(a)
            n='';
            a=abs(a);
            tol=0.00001;
            a0=fix(a+tol);
            b=a-a0;
            while(a0>0)
                n=[n num2str(mod(a0,10))];
                a0=fix(a0/10);
            end
            nt=n;
            len=length(n);
            for i=1:len
                n(i)=nt(len+1-i);
            end
            m='';
            while(b>0.00001)
                b=b*10;
                tmp=fix(b+tol);
                m=[m num2str(tmp)];
                b=b-tmp;
            end
        end
        function str=num2str_decimal(num,digits)
            %turning num into a string with digits after the decimal point
            if(~exist('digits','var'))
                digits=2;
            end
            [n,m]=CommonMethods.getDigits(num);
            len=length(m);
            if(len>digits)
                m=m(1:digits);
            else
                for i=1:digits-len
                    m=[m '0'];
                end
            end
            if(isempty(n))
                n='0';
            end
            str=n;
            if(~isempty(m))
                str=[str '.' m];
            end
            if(num < 0)
                str=['-' str];
            end
        end
        
        function fnames=getFileList_recursive(root,exts,pattern)
            folders=horzcat(root, CommonMethods.getSubfolderList_recursive(root));
            fnames={};
            if ~exist('pattern','var')
                pattern='';
            end
            if ~exist('exts','var')
                exts={};
            end
            for f=1:length(folders)
                folder=folders{f};
                fnames=horzcat(fnames, CommonMethods.getFileList(folder,exts,pattern));
            end
        end
        
        function fnames=getFileList(folder,exts,pattern)
            %folder, extension (example '.nii'), regexp string
            if(~iscell(exts))
                exts={exts};
            end
            if ~exist('pattern','var')
                pattern='';
            end
            lst=dir(folder);
            if(~exist('exts','var'))
                exts={};
            end
            len=length(lst);
            fnames={};
            for i=1:len
                tmp=lst(i).name;
                [~,name,ext]=fileparts(tmp);
                if(strcmp(name,'.')||strcmp(name,'..'))
                    continue;
                end
                if(~isempty(exts))
                    if(~any(ismember(exts,ext)))
                        continue;
                    end
                end
                fnames{end+1}=fullfile(folder,[name ext]);
            end
            if ~isempty(pattern) && ~isempty(fnames)               
                idx=cellfun(@(x)isempty(regexp(x,pattern,'ONCE')),fnames);
                fnames=fnames(idx);
            end
            idx=cellfun(@(x)exist(x,'file')==2,fnames);
            fnames=fnames(idx);
        end
        
        function p=F_Table(X,df_num,df_den)
            %look up the F table
            p=fcdf(X,df_num,df_den,'upper');
        end
        
        
        function write_vols_4d_fnames(fnames,fname4D)
            %getVoxelValues(roi_imgs, data_imgs,hold)
            %fname: a cellarray of 3D imgs
            %fname4D: the file name of the 4D output image
            img1=fnames{1};
            V1=spm_vol(img1);
            len=length(fnames);
            vols=zeros([V1.dim len]);
            for i=1:len
                img=fnames{i};
                V(i)=V1;
                V(i).fname=fname4D;
                volst=CommonMethods.getVoxelValues(img1,img,1);
                vols(:,:,:,i)=volst{1};
            end
            CommonMethods.write_vols_4d(V,vols);
        end
        
        function write_vols_4d(V, vols)
            %this function write a 4D vols to a file specified by V, a
            %structure array (length is the number of volumes). vols is a
            %4D array
            %The 4D file name will be the same as V(1).fname
            V1=V(1);
            fname=V(1).fname;
            dim=V1.dim;
            [folder, name, ext]=fileparts(fname);
            tempfolder=fullfile(folder,'temp_4D_spm_1grsfsafs');%not duplicated
            if(exist(tempfolder,'file'))
                rmdir(tempfolder, 's');
            end
            mkdir(tempfolder);
            len=length(V);
            imgs_3d=cell(1,len);
            for i=1:len
                Vt=V1;
                fn=fullfile(tempfolder, ['tmp_' num2str(i) ext]);
                Vt.fname=fn;
                Vt.mat=V(i).mat;
                imgs_3d{i}=fn;
                volst=zeros(dim(1:3));
                volst(:,:,:)=vols(:,:,:,i:i);
                spm_write_vol(Vt,volst);
            end
            load('3D_to_4D_Conversion_spm12.mat');
            matlabbatch{1}.spm.util.cat.vols=imgs_3d';
            matlabbatch{1}.spm.util.cat.name=fname;
            matlabbatch{1}.spm.util.cat.dtype=V1.dt(1);
            spm_jobman('run', matlabbatch);
            rmdir(tempfolder, 's');
        end
        
        function pops=getObjProperties(obj)
            %This function obtains all properties of obj
            names=properties(obj);
            for i=1:length(names)
                name=names{i};
                eval(['pops.' name '=obj.' name ';']);
            end
        end
        
        function [tbl, desMat, desMatColNames]=getSPMDesMatTable(spmFile)
            %this function export the design matrix in the SPM.mat as a
            %cell array of lines. The first line is the header. It's tab
            %seperated version.
            load(spmFile);
            tab=sprintf('%s\t','');
            cm=tab;
            desMat=SPM.xX.X;
            desMatColNames=SPM.xX.name;
            cols=length(desMatColNames);
            rows=size(desMat,1);
            tbl={};
            line=desMatColNames{1};
            for c=2:length(desMatColNames)
                line=[line cm desMatColNames{c}];
            end
            tbl{end+1}=line;
            
            for r=1:rows
                line=num2str(desMat(r,1));
                for c=2:cols
                    line=[line cm num2str(desMat(r,c))];
                end
                tbl{end+1}=line;
            end
        end
        
        function [tbl, contrastNames, statMapNames]=getSPMContrastTable(spmFile)
            %this function export the contrast vecters as a cell array of
            %lines. The first line is the header. It's tab
            %seperated version.
            tab=sprintf('%s\t','');
            load(spmFile);
            cm=tab;
            con=SPM.xCon;
            num=length(con);
            tbl={};
            line='name';
            vlen=length(con(1).c);
            desMatColNames=SPM.xX.name;
            for i=1:vlen
                line=[line cm desMatColNames{i}];
            end
            line=[line cm 'fname, contrast'];
            tbl{end+1}=line;
            
            statMapNames={};
            contrastNames={};
            for i=1:num
                construct=con(i);
                rows=size(construct.c,2);
                cols=size(construct.c,1);
                for r=1:rows
                    if(r==1)
                        line=construct.name;
                    else
                        line=' ';
                    end
                    for c=1:cols;
                        line=[line cm num2str(construct.c(c))];
                    end
                    contrastNames{end+1}=construct.name;
                    statMapNames{end+1}=construct.Vspm.fname;
                    line=[line cm construct.Vspm.fname ': ' construct.name];
                    tbl{end+1}=line;
                end
            end
        end
        
        function outSPM_StatSummary(spmFile, coors)
            if(~exist('voxels', 'var'))
                coors=[-7.5 34.5 27; 19.5 6 69];;
            end
            tab=sprintf('%s\t','');
            cm=tab;
            if(~exist(spmFile,'file'))
                return;
            end
            load(spmFile);
            if(~isempty(coors))
                sCoors=arrayfun(@(x)['v_' sprintf('%6.1f', coors(x,1)) ' '  sprintf('%6.1f', coors(x,2)) ' '  sprintf('%6.1f', coors(x,3))], 1:size(coors,1), 'UniformOutput', false);
            else
                sCoors=[];
            end
            out_dir=SPM.swd;
            fid=fopen(fullfile(out_dir,'SPM_StatsSummary.tsv'), 'wt');
            line=[out_dir];
            nV=size(coors,1);
            
            fprintf(fid,'%s\n',line);
            
            [tbl, conNames, conStatsFiles]=CommonMethods.getSPMContrastTable(spmFile);
            imgs=cellfun(@(x)fullfile(out_dir,x),conStatsFiles, 'UniformOutput', false);
            
            rows=length(tbl);
            
            for r=1:rows
                line=[' ' cm tbl{r}];
                if(r==1)
                    for i=1:nV
                        line=[line cm sCoors{i}];
                    end
                else
                    voxels=CommonMethods.getVoxels_1Img(imgs{r-1}, coors');
                    for i=1:nV
                        line=[line cm num2str(voxels(i))];
                    end
                end
                fprintf(fid,'%s\n',line);
            end
            
            tbl=CommonMethods.getSPMDesMatTable(spmFile);
            rows=length(tbl);
            imgs=SPM.xY.P;
            voxelsT=[];
            for r=1:rows
                line=[' ' cm ' ' cm tbl{r}];
                if(exist('scanIDs','var'))
                    if(r==1)
                        line=[line cm 'scanID'];
                    else
                        line=[line cm scanIDs{r-1}];
                    end
                end
                if(exist('imgs','var'))
                    if(r==1)
                        line=[line cm 'Img'];
                    else
                        line=[line cm imgs{r-1}];
                    end
                end
                if(r==1)
                    for i=1:nV
                        line=[line cm sCoors{i}];
                    end
                else
                    voxels=CommonMethods.getVoxels_1Img(imgs{r-1}, coors');
                    voxelsT=[voxelsT;voxels];
                    
                    for i=1:nV
                        line=[line cm num2str(voxels(i))];
                    end
                end
                fprintf(fid,'%s\n',line);
            end
            fclose(fid);
        end
        function closeAll(exceptions)
            %this function close all figures, excluding the exceptions
            %exceptions is a cellarray of figure handles
            %            hs=findall(0);
            hs=allchild(0);
            if(~exist('exceptions','var'))
                exceptions={};
            end
            if(~iscell(exceptions))
                exceptions={exceptions};
            end
            
            ne=cellfun(@(x)~isempty(x),exceptions);
            exceptions=exceptions(ne);
            len=length(exceptions);
            for i=1:length(hs)
                h=hs(i);
                excepted=false;
                for j=1:len
                    ept=exceptions{j};
                    if(ischar(ept))
                        excepted=~isempty(strfind(get(h,'Name'),ept));
                    else
                        excepted=h==ept;
                    end
                    if(excepted)
                        break;
                    end
                end
                if(excepted)
                    continue;
                end
                try
                    delete(hs(i));
                catch err
                end
            end
            system('pkill jxbrowser-chrom');
        end
        function line=applyingLineIndentation(line, level, dilimiter)
            for i=1:level
                line=[' ' dilimiter line];
            end
        end
        function tbl=applyingLineIndentationToTable(tbl, level, dilimiter)
            line=' ';
            line=CommonMethods.applyingLineIndentation(line, level-1, dilimiter);
            rows=size(tbl,1);
            col=cell(rows,1);
            col(:,1)={line};
            tbl=CommonMethods.insertColumnToTable(tbl,col,1,dilimiter);
        end
        function table=insertColumnToTable(table, col,c,dilimiter)
            %table cell array of lines, each line contain equal number of
            %dilimiters. col is a cell array of strings. The length of c and
            %lenth of the (the number of lines) table are the same
            %c is the column number of col. col will become the first
            %column if c==1, and will become the last column if c is negative.
            if(isempty(table))
                table=col;
                return;
            end
            rows=length(table);
            for r=1:rows
                line=table{r};
                pos=strfind(line,dilimiter);
                cols=length(pos)+1;
                len=length(line);
                if(c==1)%inserting as the first column
                    line1=[col{r} dilimiter line];
                elseif(c<0)%appending as the last colmn
                    line1=[line dilimiter col{r}];
                else
                    line1=[line(1:pos(c-1)) dilimiter line(pos(c-1)+1:len)];
                end
                table{r}=line1;
            end
        end
        function table=insertRowToTable(table, line0,r,dilimiter)
            %table cell array of lines, each line contain equal number of
            %dilimiters. line is a string may contain dilimiters.
            %line will be inserted at the first line if r==1, and will be
            %the last line if r<0.
            %each cell in teh table will be having the same number of
            %dilimiter
            if(isempty(table))
                table={line0};
                return;
            end
            rows=length(table);
            line=table{1};
            strs=strsplit(line,dilimiter);
            cols=length(strs);
            
            strst=strsplit(line0,dilimiter);
            cols0=length(strst);
            if(cols>cols0)
                for i=1:(cols-cols0)
                    line0=[line0 dilimiter ' '];
                end
            elseif(cols0>cols)
                st=' ';
                for i=1:(cols0-cols-1)
                    st=[st dilimiter ' '];
                end
                col=cell(rows,1);
                col(:,1)={st};
                table=CommonMethods.insertColumnToTable(table, col,-1,dilimiter);
            end
            if(r==1)
                table=vertcat(line0,table);
            elseif(r<0)
                table=vertcat(table,line0);
            else
                table=vertcat(table(1:r-1), line0, table(r: rows));
            end
        end
        function tbl=mergeTables(tbl1,tbl2,dilimiter)
            len=length(tbl1);
            tbl=tbl2;
            for i=1:len
                tbl=CommonMethods.insertRowToTable(tbl,tbl1{len-i+1}, 1, dilimiter);
            end
        end
        function tbl=horzcatTables(tbl1,tbl2)
            len1=length(tbl1);
            len2=length(tbl2);
            if(len1==0||len2==0)
                len1=len1;
            end
            strs1=cellfun(@(x)strsplit(x,CommonMethods.tab),tbl1, 'UniformOutput', false);
            strs2=cellfun(@(x)strsplit(x,CommonMethods.tab),tbl2, 'UniformOutput', false);
            if(len2>len1)
                for l=len1+1:len2
                    tbl1=CommonMethods.insertRowToTable(tbl1, ' ',-1,CommonMethods.tab);
                end
            else
                for l=len2+1:len1
                    tbl2=CommonMethods.insertRowToTable(tbl2, ' ',-1,CommonMethods.tab);
                end
            end
            len=length(tbl1);
            tbl=cell(len,1);
            for l=1:len
                tbl{l}=[tbl1{l} CommonMethods.tab tbl2{l}];
            end
            strs=cellfun(@(x)strsplit(x,CommonMethods.tab),tbl, 'UniformOutput', false);
        end
        function testInsertRow()
            %dil=sprintf('%s\n','');
            dil=',';
            tbl={['sfds1' dil 'fsafsf' dil 'fsfaf']; ['sfds2' dil 'fsafsf' dil 'fsfaf']; ['sfds3' dil 'fsafsf' dil 'fsfaf']; ['4sfds' dil 'fsafsf' dil 'fsfaf']};
            line='line';
            tbl=CommonMethods.insertRowToTable(tbl,line,1,dil);
            tbl=CommonMethods.insertRowToTable(tbl,line,-1,dil);
            tbl=CommonMethods.insertRowToTable(tbl,line,4,dil);
            line=['line2' dil ' ' dil ' ' dil ' ' dil ' '];
            tbl=CommonMethods.insertRowToTable(tbl,line,4,dil);
        end
        function newStr=insertAfter(str,pos,str1)
            newStr=[str(1:pos) str1 str(pos+1:length(str))];
        end
        function str=cell2str(strs,connector)
            len=length(strs);
            if(len<1)
                str='';
                return;
            end
            str=strs{1};
            for i=2:len
                str=[str connector strs{i}];
            end
        end
        function cls=cell2tablecells(strs,pos,tw,strEmpty)
            %this function return a cellarray of strings.
            %tw: table width, pos: position of the first element in strs,
            %strEmpy: filling the space with s
            cls=cell(1,tw);
            if(ischar(strs))
                strs={strs};
            end
            for i=1:pos-1
                cls{i}=strEmpty;
            end
            for i=1:length(strs);
                ii=pos+i-1;
                if(ii>tw)
                    break;
                end
                cls{ii}=strs{i};
            end
            for i=pos+length(strs):tw
                cls{i}=strEmpty;
            end
        end
        function strs=replaceEmptyCells(strs0,emp)
            %this function replaces every empty cell in the string cell
            %array strs wit emp
            strs0(cellfun(@(x)isempty(x),strs0))={emp};
            strs=strs0;
        end
        function line= cell2tableline(strs,pos,tw,del,strSpace)
            %this function return a delimitted line using a cellarray of
            %strings.
            %tw: table width, pos: position of the first element in strs,
            %strSpace: filling the space with strSpace
            if(ischar(strs))
                strs={strs};
            end
            if(pos==1)
                pos0=2;
                line=strs{1};
            else
                line=strSpace;
                pos0=1;
                for i=2:(pos-1)
                    line=[line del strSpace];
                end
            end
            for i=pos0:length(strs);
                line=[line del strs{i}];
            end
            for i=pos+length(strs):tw
                line=[line del strSpace];
            end
        end
        
        function h=getSPM_extension()
            h=CommonMethods.getFigureByName('spm_extension');
        end
        
        function h=getFigureByName(name0)
            len=length(name0);
            h=[];
            hs=allchild(0);
            for i=1:length(hs)
                name=hs(i).Name;
                if(length(name)<len)
                    continue;
                end
                if(strcmpi(name(1:len), name0))
                    h=hs(i);
                    break
                end
            end
        end
        
        function str=getSelectedString(hObject)
            strs=get(hObject,'String');
            str=strs{get(hObject,'Value')};
        end
        
        function spmExt=updateSPMExt(coorStr,fnameSPM)
            
            xjv=xjViewHandler.get_xjView();
            if(~isempty(xjv))
                set(0,'CurrentFigure',xjv);
                spm_orthviews('reposition');
            else
                return;
            end
            
            spmExt=CommonMethods.getFigureByName('spm_extension');
            
            spmExtGuiData=guidata(spmExt);
            if(exist('fnameSPM','var'))
                fnameSPM=strtrim(fnameSPM);
                if(exist(fnameSPM,'file'))
                    exception={fnameSPM};
                    spmExtGuiData.clearAssociatedObjects(spmExt,exception);
                    spmExtGuiData.importAssociatedObject(spmExt,fnameSPM);
                end
            end
            coor=CommonMethods.str2nums(coorStr, ' ')';
            if(length(coor)~=3)
                return;
            end
            spmExtGuiData.updateCoordinates_mni(spmExt,coor(1:3));
        end
        function interval=getInterval(date1, date2)
            interval=datenum(date2)-datenum(date1);
        end
        function y=nan2mean(x)
            %this function replace NaN in X with the mean of other elements
            %in x.
            selection=~isnan(x);
            m=mean(x(selection));
            y=x;
            y(~selection)=m;
        end
        function changeSPMPath(fnameSPM,foldero,foldern)
            load(fnameSPM);
            SPM.xY.P=cellfun(@(x)strrep(x,foldero,foldern),SPM.xY.P,'UniformOutput',false);
            save(fnameSPM,'SPM');
        end
        function z=getFisherZ(r)
            z = 0.5*log((1+r)/(1-r));
        end
        
        function showVolComparison(vols, inds, par)
            %this function display orthogonal views of vols (cellarray of
            %3D matrices) at position inds (array of 1D or 3D indexes).
            %par: a structure specifies display options
            %inds should be a column vector
            len=length(inds);
            if(~exist('par','var'))
                par=[];
            end
            if(~isfield(par,'zoom'))
                par.zoom=1;
            end
            if(~isfield(par,'margin'))
                par.margin=5;
            end
            if(~isfield(par,'ws'))
                par.ws=20;
            end
            if(~isfield(par,'w'))
                par.w=floor(sqrt(len));
            end
            if(~isfield(par,'mag'))
                par.mag=5;
            end
            if(~isfield(par,'fh'))
                par.fh=figure;
            end
            if(~isfield(par,'h'))
                par.h=ceil(sqrt(len));
                if(par.h*par.w<len)
                    par.h=par.h+1;
                end
            end
            
            lent=2*par.ws+1;
            no=0;
            for i=1:par.h
                for j=1:par.w
                    no=no+1;
                    if(no>len)
                        continue;
                    end
                    m=CommonMethods.getOrthviews(vols{no},inds{no},par.ws,par.margin);
                    %                   m=CommonMethods.magnify(m,par.mag);
                    subplot(par.h,par.w,no);
                    imshow(mat2gray(m));
                end
            end
        end
        function revreseNormalize(img_inMNI,DeformFieldImgs)
            selection=cellfun(@(x)exist(x,'file')>0,DeformFieldImgs);
            DeformFieldImgs=DeformFieldImgs(selection);
            %fname='/home/tjin/MyProjects/MatlabProjects_2017b/MRI/Scripts/Jobtemplates/sAC_PC_Marks.nii';
            load('/home/tjin/MyProjects/MatlabProjects_2017b/MRI/Scripts/Jobtemplates/Deformation_pullback.mat');
            
            len=length(DeformFieldImgs);
            matlabbatch{1}.spm.util.defs.out{1}.pull.fnames={img_inMNI};
            jobs=cell(1,len);
            jobs(:)=matlabbatch;
            parfor i=1:len
                %for i=1:1
                matlabbatch=jobs(i);
                cprintf('blue','%s\n',['The ' num2str(i) '-th (of ' num2str(len) ')']);
                defImg=DeformFieldImgs{i};
                folder=fileparts(defImg);
                cd(folder);
                %    continue;
                matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr={folder};
                matlabbatch{1}.spm.util.defs.comp{1}.def={defImg};
                try
                    spm_jobman('run', matlabbatch);
                catch err
                    cprintf('red','%s\n',folder);
                end
            end
        end
        function makeAC_PC_imgs(imgs)
            spm('defaults', 'FMRI');
            if(true)
                len=length(imgs);
                DeformFieldImgs=cell(1,len);
                for i=1:len
                    img=imgs{i};
                    [folder,name,ext]=fileparts(img);
                    DeformFieldImgs{i}=fullfile(folder,['iy_' name ext]);
                end
                selection=cellfun(@(x)exist(x,'file')>0,DeformFieldImgs);
                %DeformFieldImgs=DeformFieldImgs(selection);
                fname='/home/tjin/MyProjects/MatlabProjects_2017b/MRI/Scripts/Jobtemplates/sAC_PC_Marks.nii';
                load('/home/tjin/MyProjects/MatlabProjects_2017b/MRI/Scripts/Jobtemplates/Deformation_pullback.mat');
                
                len=length(DeformFieldImgs);
                matlabbatch{1}.spm.util.defs.out{1}.pull.fnames={fname};
                jobs=cell(1,len);
                jobs(:)=matlabbatch;
                parfor i=1:len
                    %for i=1:1
                    matlabbatch=jobs(i);
                    cprintf('blue','%s\n',['The ' num2str(i) '-th (of ' num2str(len) ')']);
                    defImg=DeformFieldImgs{i};
                    folder=fileparts(defImg);
                    cd(folder);
                    %    continue;
                    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr={folder};
                    matlabbatch{1}.spm.util.defs.comp{1}.def={defImg};
                    try
                        spm_jobman('run', matlabbatch);
                    catch err
                        cprintf('red','%s\n',folder);
                    end
                end
            end
            
            masks=cellfun(@(x)fullfile(fileparts(x), 'wsAC_PC_Marks.nii'), imgs, 'UniformOutput', false);
            len=length(masks);
            for i0=1:length(masks)
                cprintf('blue','%s\n',['The ' num2str(i0) '-th (of ' num2str(len) ')']);
                mask=masks{i0};
                fname_maskVolInds=fullfile(fileparts(mask),'maskVolInds.txt');
                vols=spm_read_vols(spm_vol(mask));
                vols(isnan(vols))=0;
                lx=imregionalmax(vols);
                inds=find(lx);
                [x,y,z]=ind2sub(size(vols),inds);
                vals=vols(lx);
                [~,inds]=sort(vals);
                volInds=[x(inds) y(inds) z(inds)];%vol inds for AC, PC and the midplane point
                fid=fopen(fname_maskVolInds,'wt');
                for i=1:3
                    line=[num2str(volInds(i,1)) ' ' num2str(volInds(i,2)) ' ' num2str(volInds(i,3))];
                    fprintf(fid,'%s\n',line);
                end
                fclose(fid);
            end
            
        end
        
        function img=applyRigidTransform(img0, R4, affix,overwrite)
            %this function applys a rigid tranform (R4) to the image
            if(~exist('affix','var'))
                affix='';
            end
            if ~exist('overwrite','var')
                overwrite=true;
            end
            
            [folder,name,ext]=fileparts(img0);
            img=fullfile(folder,['R4_' name affix ext]);
            
            if exist(img,'file') && ~overwrite
                return;
            end
            
            V=spm_vol(img0);           
            vols=spm_read_vols(V);
            
            V.fname=img;
            V.mat=R4*V.mat;
            spm_write_vol(V,vols);
        end
        
        function separateClusters(img)
            [folder,name,ext]=fileparts(img);
            V=spm_vol(img);
            vols=spm_read_vols(V);
            idx=find(vols>0);
            dim=size(vols);
            [ix,iy,iz]=ind2sub(dim,idx);
            coors=[ix iy iz];
            cluster_coors=CommonMethods.findClusters(coors');
            num=length(cluster_coors);
            for i=1:num
                cur_vols=zeros(dim);
                ixyz=cluster_coors{i}';
                inds=sub2ind(dim,ixyz(:,1),ixyz(:,2),ixyz(:,3));
                cur_vols(inds)=1;
                V.fname=fullfile(folder,[name '_cluster' num2str(i) '_' num2str(length(inds)) 'voxels' ext]);
                spm_write_vol(V,cur_vols);
            end
        end
        
        function changeYData(hFig,Ydata0, Ydata)
        end
        function data=getFigData(hFig,Ydata0, Ydata)
            axesObjs = get(hFig, 'Children');
            dataObjs = get(axesObjs, 'Children');
            objTypes = get(dataObjs, 'Type');
            xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
            ydata = get(dataObjs, 'YData');
            zdata = get(dataObjs, 'ZData');
            data=[xdata ydata zdata];
        end
        function tpm=get_spm_templates()
            tpm.TPM=fullfile(CommonEnvMethods.getMatlabProjectPath(),'External Software','spm12','tpm','TPM.nii');
            tpm.avg152T1=fullfile(CommonEnvMethods.getMatlabProjectPath(),'External Software','spm12','canonical','avg152T1.nii');
            tpm.single_subj_T1=fullfile(CommonEnvMethods.getMatlabProjectPath(),'External Software','spm12','canonical','single_subj_T1.nii');
        end
        function bw=boundaryMask3(L)
            %bw 1 for boundary
            persistent strel;
            dim=size(L);
            bw=zeros(dim);
            if(isempty(strel))
                strel=cell(1,26);
                num=0;
                for x=-1:1
                    for y=-1:1
                        for z=-1:1
                            if(x==0&&y==0&&z==0)
                                continue;
                            end;
                            num=num+1;
                            strel{num}=[x y z];
                        end
                    end
                end
                %put [0 0 0] to the first position;
                temp=strel{14};
                strel{14}=strel{1};
                strel{1}=temp;
            end
            
            for x=1:dim(1)
                for y=1:dim(2)
                    for z=1:dim(3)
                        pos0=[x y z];
                        val0=L(x,y,z);
                        for i=1:26
                            pos=pos0+strel{i};
                            if(any(pos<=0)||any(pos>dim))
                                continue;
                            end
                            val=L(pos(1),pos(2),pos(3));
                            if(val0~=val)
                                bw(x,y,z)=1;
                                break;
                            end
                        end
                    end
                end
            end
        end
        function [volInds, coors]=getMaskCenters(img)
            %return the centers (topological) of masks in img, assumed to
            %be banary image.
            if ischar(img)
                V=spm_vol(img);
                mat=V.mat;
                img=spm_read_vols(V);
            end
            dm=bwdist(~img>0);
            dim=size(img);
            [~,MaxPos0,~,MinPos0]=CommonMethods.localExtrema(dm);
            volInds0=MaxPos0;
            %removing duplicated 
            num=size(volInds0,1);
            regions=zeros(dim);
            volInds={};
            coors={};%world coordinates
            for i=1:num
                pos=volInds0(i,:);
                if(regions(pos(1),pos(2),pos(3)))
                    continue;
                end
                [~,m]=regionGrowing(dm,pos,0.);
                poss=find(m>0);%finding the center of the positions have equal dist matrix values
                [x,y,z]=ind2sub(dim,poss);
                
                cnt=CommonMethods.findCenter([x y z]');
                volInds{end+1}=cnt;
                %marking the mask so that not finding other local maxima of
                %dist matrix as the center.
                [~,m]=regionGrowing(img,pos,0.);
                regions(m>0)=i;
                c=mat*[cnt;1];
                coors{end+1}=c(1:3);
            end
        end
        function [uvars,tbl]=buildDummyVarTable(vars)
            %This function convert a list of categorical variables into a
            %dummy variable table. Each column of the bable is a category
            %label, and each row indicate indicates the presence of the
            %label using 0 (Not present) and 1 (present)
            
            uvars=unique(vars);
            cols=length(uvars);
            rows=length(vars);
            tbl=nan(rows,cols);
            for c=1:cols
                tbl(:,c)=ismember(vars,uvars{c});
            end
        end
        function gunzipImages(folder)
            if ~exist('folder','var')
                folder=pwd;
            end
            lst=CommonMethods.getFileList(folder,'.gz');
            for i=1:length(lst)
                fname=lst{i};
                strs=strsplit(fname,'.');
                if ~strcmp(strs{end-1},'nii')
                    continue;
                end
                gunzip(fname);
            end
        end
        function testGLM()
            load fisheriris
            X = meas(51:end,:); % versicolor and virginica
            y = strcmp('versicolor',species(51:end));
            mdl = fitglm(X,y,'linear',...
                'distr','binomial');
            a=1;
        end
        function [inds, indsr]=getReflectionIndexes(dim)
            %returning 1-D index of array with dimenssionarity dim and the
            %reflection along the first axis
            tmp=ones(dim);
            inds=find(tmp);
            [x,y,z]=ind2sub(dim,inds);
            len=dim(1);
            xr=arrayfun(@(x)(len-x+1),x);
            indsr=sub2ind(dim,xr,y,z);
            %[xr,yr,zr]=ind2sub(dim,indsr);
        end
        function Vr=getReflection_3D(V)
            dim=size(V);
            Vr=zeros(dim);
            [inds,indsr]=CommonMethods.getReflectionIndexes(dim);
            Vr(inds)=V(indsr);
        end
        function testFreeSurfer()
            fnamev='C:\Users\tjin\MyProjects\FreeSurfer\FryerCDA\SM1001\mri\brainmasko.nii';
            fnames_p='C:\Users\tjin\MyProjects\FreeSurfer\FryerCDA\SM1001\surf\lh.pial';
            fnames_w='C:\Users\tjin\MyProjects\FreeSurfer\FryerCDA\SM1001\surf\lh.white';
            folder='C:\Users\tjin\MyProjects\FreeSurfer\FryerCDA\SM1001\surf';
            
            fnames=fnames_p;
            [v, f] = freesurfer_read_surf(fnames);
            V=spm_vol(fnamev);
            inds=CommonMethods.mni2ind_cols(v', V.mat)';
            inds1=sub2ind(V.dim,inds(:,1),inds(:,2),inds(:,3));
            vols=zeros(V.dim);
            vols(inds1)=1;
            V.fname=fullfile(folder,'lh_pial.nii');
            spm_write_vol(V,vols);
            
            fnames=fnames_w;
            [v, f] = freesurfer_read_surf(fnames);
            V=spm_vol(fnamev);
            inds=CommonMethods.mni2ind_cols(v', V.mat)';
            inds1=sub2ind(V.dim,inds(:,1),inds(:,2),inds(:,3));
            vols=zeros(V.dim);
            vols(inds1)=1;
            V.fname=fullfile(folder,'lh_white.nii');
            spm_write_vol(V,vols);
        end
     end
end