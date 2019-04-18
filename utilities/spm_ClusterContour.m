classdef spm_ClusterContour <handle
    %SPM_CLUSTERCONTOUR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        st;%an spm structure;
        coor;%the current point;
        voxelIndex;%the voxel index of the current point;
        coors;%a collection of voxel coordinates, m by 3 matrix;
        mat;%voxel indexes to coordinates transform matrix;
        %This class will draw contours of objects filling coors. The
        %contours will be drawn on the three orth planes crossing at coor.
        voxelIndexes%The voxel indexes of coors
        ContourLines;%contour lines, a cellarry (1,3) for the line objects on the three axes showing the three orth views
        %the three lines corresponding to the three axis in st.vols{1}.ax
        ContourLines_vInd;%contour lines in voxel indexes;
        ContourLines_coor;%contour lines in coordinates
        TransformPars;%the parameters needed to transform contour lines in
        %world coordinates (such as mni) to the axis units of the three
        %axes
    end
    
    methods
        function obj=spm_ClusterContour()
            obj.ContourLines={};
        end
        function update(obj, st, coor, coors,mat)
            obj.st=st;
            obj.coor=coor;
            obj.coors=coors;
            obj.mat=mat;
            obj.buildContourLines();
            obj.drawContourLines();
        end
        
        function buildContourLines(obj)
            if(isempty(obj.coors))
                return;
            end
            obj.voxelIndexes=CommonMethods.mni2ind_cols(obj.coors', obj.mat)';
            obj.voxelIndex=CommonMethods.mni2ind_cols(obj.coor', obj.mat)';
            obj.ContourLines_vInd=cell(1,3);
            obj.ContourLines_coor=cell(1,3);
            os=[1 2 3];
            for a=1:3%looping through the three axes
                inds=obj.voxelIndexes(obj.voxelIndexes(:,a)==obj.voxelIndex(a),:);
                if(isempty(inds))
                    continue;
                end
                os2=os(os~=a);
                inds2=[inds(:,os2(1)) inds(:,os2(2))]; %the 2D voxel indexes on the plane defined by the current point and the axis a
                
                if(size(inds2,1)==1)
                    maxs=max([inds2;inds2]);
                else
                    maxs=max(inds2);
                end
                
                tmp=zeros(maxs(1)+2,maxs(2)+2);
                for t=1:size(inds2,1)
                    tmp(inds2(t,1),inds2(t,2))=1;
                end
                
                level=[1 1];
                cline=contourc(tmp,level);
                %contourc treat row number as x and column number as x, and
                %will be the first and the second rows in cline,
                %respectively. This will be reflected in how to insert
                %these voxel indexes into the 3D contour line of voxel
                %indexes.
                idx=find(cline(1,:)==level(1));
                len=length(idx);%number of contour lines
                clines_vind=cell(1,len);
                clines_coor=cell(1,len);
                for i=1:len
                    ind=idx(i);
                    len1=cline(2,ind);
                    vinds=zeros(len1,3);
                    vinds(:,a)=obj.voxelIndex(a);
                    vinds(:,os2(1))=cline(2,ind+1:ind+len1);
                    vinds(:,os2(2))=cline(1,ind+1:ind+len1);
                    %
                    clines_vind{i}=vinds;
                    coorline=CommonMethods.ind2coor_mat(vinds, obj.mat);
                    clines_coor{i}=coorline;
                end
                obj.ContourLines_vInd{a}=clines_vind;
                obj.ContourLines_coor{a}=clines_coor;
            end
        end
        function clearContourLines(obj)
            if(isempty(obj.ContourLines))
                return;
            end
            for a=1:3
                clines=obj.ContourLines{a};
                for l=1:length(clines)
                    delete(clines{l});
                end
            end
            obj.ContourLines={};
        end
        function drawContourLines(obj)
            obj.clearContourLines();
            obj.ContourLines=cell(1,3);
            bb   = obj.st.bb;
            Dims = round(diff(bb)'+1);
            is   = inv(obj.st.Space);
            cent = is(1:3,1:3)*obj.st.centre(:) + is(1:3,4);
            
            clines_coor=obj.ContourLines_coor{3};
            len=length(clines_coor);
            clines=cell(1,len);            
            axis=obj.st.vols{1}.ax{1}.ax;
            for l=1:len
                cline_coor=clines_coor{l};
                len1=size(cline_coor,1);
                xData=zeros(1,len1);
                yData=zeros(1,len1);
                for i=1:len1
%                    cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
                    cent = is(1:3,1:3)*cline_coor(i,:)' + is(1:3,4);
                    xData(i)=cent(1)-bb(1,1)+1;
                    yData(i)=cent(2)-bb(1,2)+1;
                end
                clines{l}=line(axis,xData,yData, 'LineWidth', 0.5, 'Color', [0 0 1]);
%                     set(st.vols{i}.ax{1}.lx,'HitTest','off',...
%                         'Xdata',[0 TD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
%                     set(st.vols{i}.ax{1}.ly,'HitTest','off',...
%                         'Ydata',[0 TD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));                                        
            end
            obj.ContourLines{3}=clines;
%            drawnow;
%            return;

            clines_coor=obj.ContourLines_coor{1};
            len=length(clines_coor);
            clines=cell(1,len);            
            axis=obj.st.vols{1}.ax{2}.ax;
            for l=1:len
                cline_coor=clines_coor{l};
                len1=size(cline_coor,1);
                xData=zeros(1,len1);
                yData=zeros(1,len1);
                for i=1:len1
%                    cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
                    cent = is(1:3,1:3)*cline_coor(i,:)' + is(1:3,4);
                    xData(i)=cent(1)-bb(1,1)+1;
                    yData(i)=cent(3)-bb(1,3)+1;
                end
%                     set(st.vols{i}.ax{2}.lx,'HitTest','off',...
%                         'Xdata',[0 CD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
%                     set(st.vols{i}.ax{2}.ly,'HitTest','off',...
%                         'Ydata',[0 CD(2)]+0.5,'Xdata',[1 1]*(cent(1)-bb(1,1)+1));
                clines{l}=line(axis,xData,yData, 'LineWidth', 0.5, 'Color', [0 0 1]);
            end
            obj.ContourLines{2}=clines;
                        
            clines_coor=obj.ContourLines_coor{2};
            len=length(clines_coor);
            clines=cell(1,len);            
            axis=obj.st.vols{1}.ax{3}.ax;
            for l=1:len
                cline_coor=clines_coor{l};
                len1=size(cline_coor,1);
                xData=zeros(1,len1);
                yData=zeros(1,len1);
                for i=1:len1
                    %                    cent = is(1:3,1:3)*st.centre(:) + is(1:3,4);
                    cent = is(1:3,1:3)*cline_coor(i,:)' + is(1:3,4);
                    if obj.st.mode ==0
                        xData(i)=cent(3)-bb(1,3)+1;
                        yData(i)=cent(2)-bb(1,2)+1;
                        %                         set(st.vols{i}.ax{3}.lx,'HitTest','off',...
                        %                             'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(2)-bb(1,2)+1));
                        %                         set(st.vols{i}.ax{3}.ly,'HitTest','off',...
                        %                             'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(cent(3)-bb(1,3)+1));
                    else
                        xData(i)=bb(2,2)+1-cent(2);
                        yData(i)=cent(3)-bb(1,3)+1;
                        %                         set(st.vols{i}.ax{3}.lx,'HitTest','off',...
                        %                             'Xdata',[0 SD(1)]+0.5,'Ydata',[1 1]*(cent(3)-bb(1,3)+1));
                        %                         set(st.vols{i}.ax{3}.ly,'HitTest','off',...
                        %                             'Ydata',[0 SD(2)]+0.5,'Xdata',[1 1]*(bb(2,2)+1-cent(2)));
                    end
                end
                clines{l}=line(axis,xData,yData, 'LineWidth', 0.5, 'Color', [0 0 1]);
            end 
            obj.ContourLines{1}=clines;
            drawnow;
        end
    end
end

