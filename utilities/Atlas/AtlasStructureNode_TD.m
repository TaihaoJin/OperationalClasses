classdef AtlasStructureNode_TD
    %ATLASSTRUCTURENODE_TD Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        level;%1 for hemisphere, 2 lobe, 3 gyrus ..
        type;%name of the level: hemishere, lobe ...
        dim;
        mat;
        anatomy;%the name of the structure
        index;%the index to the anatomical region in the TD Atlas
        voxelIndexes;%3D indexes, cell array, 1: left, 2: inner, 3: rihgt
        voxelIndexes1;%the 1D indexes
        coor;%mni coors
        overlappingStructures;
        %cell array. The length of the cellarray is the same as the level
        %of the atlas. Each element contain the structures (at each level) overlapping with
        %"this".         
        size;%number of vexels, left, center and right counted separately. 
        size_selected;
        overlappingStructures_selected;
        voxelIndexes1_selected;
    end
    
    methods
        function obj=AtlasStructureNode_TD(tdaHandler, level, index)
            obj.level=level;
            obj.index=index;
            obj.mat=tdaHandler.mat;
            obj.dim=tdaHandler.dim;
            TD=tdaHandler.Atlas_TD;
            obj.type=TD{level}.type;
            mnilist=TD{level}.mnilist;
            obj.anatomy=TD{level}.anatomy{index};
            
            inds=find(mnilist==index);
            [idx1,idx2,idx3]=ind2sub(obj.dim,inds);
            vinds=[idx1 idx2 idx3];
            mni=CommonMethods.ind2coor_mat(vinds, obj.mat);
            xs=mni(:,1);
            ps=find(xs>0);
            ns=find(xs<0);
            zs=find(xs==0);
            obj.voxelIndexes={vinds(ns,:) vinds(zs,:) vinds(ps,:)};
            obj.coor={mni(ns,:) mni(zs,:) mni(ps,:)};
            obj.voxelIndexes1={inds(ns,:) inds(zs,:) inds(ps,:)};
            obj.size=[nnz(ns) nnz(zs) nnz(ps)];
            
            len=length(TD);
            obj.overlappingStructures=cell(1,len);
            for l=1:len
                mnilist=TD{l}.mnilist;
                oss=cell(1,3);%overlapping structures
                for i=1:3%left, middle, and right
                    sinds=arrayfun(@(x)mnilist(x),obj.voxelIndexes1{i});
                    u=unique(sinds);
                    lent=length(u);
                    osst=cell(1,lent);
                    for s=1:lent
                        sind=u(s);
                        if(sind > 0)
                            osst{s}=struct('level', l, 'sindex', sind, 'anatomy', TD{l}.anatomy{sind}, 'noVoxels', nnz(sinds==sind));
                        else
                            osst{s}=struct('level', l, 'sindex', sind, 'anatomy', 'undefined', 'noVoxels', nnz(sinds==sind));
                        end
                    end
                    oss{i}=osst;
                end
                obj.overlappingStructures{l}=oss;
            end
        end
        function updateSelection(obj,tdaHandler)            
            %this method update the selection status (whether the voxels are selected or not) of the voxels within
            %the anatomical region based on the 3D selection matrix of
            %tdaHandler. 
           obj.size_selected=zeros(1,3);
            obj.voxelIndexes1_selected=cell(1,3);
            for i=1:3
                inds=obj.voxelIndexes1{i};
                selected=find(tdaHandler.selectedMniList(inds));
                obj.voxelIndexes1_selected{i}=inds(selected);
                obj.size_selected(i)=length(selected);
            end
            
            len=length(tdaHandler.Atlas_TD);
            obj.overlappingStructures=cell(1,len);
            for l=1:len
                mnilist=tdaHandler.Atlas_TD{l}.mnilist;
                oss=cell(1,3);%overlapping structures
                for i=1:3
                    inds=obj.voxelIndexes1_selected{i};
                    if(isempty(inds))
                        continue;
                    end
                    sinds=arrayfun(@(x)mnilist(x),inds);
                    u=unique(sinds);
                    lent=length(u);
                    osst=cell(1,lent);
                    for s=1:lent
                        sind=u(s);
                        if(sind > 0)
                            osst{s}=struct('level', l, 'sindex', sind, 'anatomy', tdaHandler.Atlas_TD{l}.anatomy{sind}, 'noVoxels', nnz(sinds==sind));
                        else
                            osst{s}=struct('level', l, 'sindex', sind, 'anatomy', 'undefined', 'noVoxels', nnz(sinds==sind));
                        end
                    end
                    oss{i}=osst;
                end
                obj.overlappingStructures{l}=oss;
            end
        end
        function IS=IsSelected(obj)
            IS=sum(obj.size_selected)>0;
        end
        function st=getParentStructure(obj,level)
            st=cell(1,3);%the main overlapping structures with the left, middle and right part of the structure
            if level==obj.level
                st={obj,obj,obj};
            elseif level < obj.level
                for ind_sideness=1:3 %left, center, and right
                    ps=[];
                    osst=obj.overlappingStructures{level}{ind_sideness};%overlapping structure with the left
                    if ~isempty(osst)
                        vols=cellfun(@(x)x.noVoxels, osst);
                        [~,I]=sort(vols,'descend');
                        ps=osst{I(1)};
                    end
                    st{ind_sideness}=ps
                end
            end                    
        end
        function IS=IsSubStructureOf(obj,AnatomyName,level)
            IS=false;
            st=obj.getParentStructure(level);
            for i=1:3
                stp=st{i};
                if isempty(stp)
                    continue;
                end
                if strcmp(stp.anatomy,AnatomyName)
                    IS=true;
                    break
                end
            end
        end
    end
    
end

