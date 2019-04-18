function [ output_args ] = exportSignificantTables( input_args )
%EXPORTSIGNIFICANTTABLES Summary of this function goes here
%   Detailed explanation goes here
imgNames={'R_ALFF' 'R_fALFF' 'ReHo27' 'ReHo19' 'ReHo7'};
norms={'' 'std' 'zScore'};
for imgn=1:length(imgNames)
    for nm=1:length(norms)
        nr=norms{nm};
        if(isempty(nr))
            type=imgNames{imgn};
        else
            type=[imgNames{imgn} '_' norms{nm}];
        end
        if(nm>1||imgn>1)
%            continue;
        end
        for fi=1:8
            %    fname=fullfile('/home1/data/IRB13278/data/Resting/group/analysisLogs','SPM_Summaries_Rest_ALFF_fALFF_ReHo_GroupAnalysis_12-Sep-2017-19-04-11_log1.tsv');
            name=['SPM_Summaries_Rest_ALFF_fALFF_ReHo_GroupAnalysis_12-Sep-2017-19-04-11_log' num2str(fi) '.tsv'];
            folder='/home1/data/IRB13278/data/Resting/group/analysisLogs';
            fname=fullfile(folder,name);
            fnameo=fullfile(folder,['significantTables_' type '.tsv']);%output file
%            fnameo=fullfile(folder,['significantTables_' 'all' '.tsv']);
            fo=fopen(fnameo,'at');
            fname=fullfile(folder,name);
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
                fnames{i}=strs{1};
            end
            %           fnames=fnames(cellfun(@(x)~isempty(strfind(x,filesep)), fnames));
            idx=cellfun(@(x)~isempty(strfind(x,filesep)), fnames);
            inds=find(idx);
            
            len=length(inds);
            lent=length(lines);
            for i0=1:len
                lI=inds(i0);
                if(i0<len)
                    lF=inds(i0+1)-1;
                else
                    lF=lent;
                end
                linest=lines(lI:lF);%the lines for a single SPM.mat
                head=linest{1};
                strs=strsplit(head,filesep);
                imgType=strs{8};
                if(~strcmp(imgType, type))
                    continue;
                end
                fprintf(fo,'%s\n',head);
                idxt=cellfun(@(x)~isempty(strfind(x,'spmT_')),linest);
                indst=find(idxt);
                
                numcon=length(indst);
                lentt=length(linest);
                nums=0;
                for c=1:numcon
                    lI=indst(c);
                    if(c<numcon)
                        lF=indst(c+1)-1;
                    else
                        lF=lentt;
                    end
                    if((lF-lI)<=2)
                        continue;%no significant clusters;
                    end
                    cdisp=linest{lI};
                    if(~isempty(strfind(cdisp,'age'))||~isempty(strfind(cdisp,'ticv'))||~isempty(strfind(cdisp,'all')))
                        continue;
                    end
                    nums=nums+1;
%                     if(nums==1)
%                         fprintf(fo,'%s\n',head);
%                     end
                    
                    strst=strsplit(linest{lI},tab);
                    strst{3}=imgType;
                    line='';
                    for s=1:length(strst)
                        line=[line tab strst{s}];
                    end
                    fprintf(fo,'%s\n',line);
                    fprintf(fo,'%s\n',[' ' tab ' ' tab linest{lI+1}]);
                    for lt=lI+2:lF
                        fprintf(fo,'%s\n',[' ' tab linest{lt}]);
                    end
                    fprintf(fo,'%s\n',' ');
                end
                %                fprintf(fo,'%s\n','');
            end
            fclose(fo);
        end
    end
end
end

