classdef ClusterHandler < handle
    %CLUSTERHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mni; %mni coordinates. It is a 3xcols matrix. each column is one 3d coordinate.
        mat; %volIndex to mni conversion matrix
        %dim depends on the range of mni, its a minimum volumes that can
        %hold all coordinates in mni.;
        vol;%the 3D data matrix;
        vClusterIndex;%a 3d array (dim) holding cluster indexes.
        %the offset of vol indexes;
        indexOffset;
        shapeHandler;
        stats;% the cluster statistics
        vMorphDist;% label 
    end
    
    methods
        function obj=ClusterHandler(mni, mat, vol)
            obj.mni=mni;
            obj.mat=mat;
            obj.vol=vol;
            obj.update();
        end
        
        function update(obj,mni, mat, vol)
            if(~exist('mni','var'))
                mni=obj.mni;
                mat=obj.mat;
                vol=obj.vol;
            end
            
            if(size(mni, 1)~=3)
                error('mni for ''ClusterHandler'' should be a 3xcols matrix.');
            end
            
            obj.buildStats(mni, mat, vol);
        end
        
        function [mnis, center, morphDist] = getClusterMNIs(obj, mni)
            %this function returns the mni coordinates of the cluster
            %including mni.
            mnis=[];
            center=[];
            ind=CommonMethods.mni2ind(mni, obj.mat);
            cind=obj.vClusterIndex(ind(1),ind(2),ind(3));
            if(cind~=-1)
                mnis=obj.stats.mnis{cind};
                center=obj.stats.centers(:,cind);
            end
            coor=CommonMethods.mni2ind(mni,obj.mat);
            
            morphDist=obj.vMorphDist(coor(1), coor(2), coor(3));
        end
        
        function cind = getClusterIndx(obj,mni)
            ind=CommonMethods.mni2ind(mni, obj.mat);
            cind=obj.vClusterIndex(ind(1),ind(2),ind(3));
        end
        
        function num = getSize(obj)
            num=size(obj.mni,2);
        end
        
        function buildStats(obj,mni, mat, vol)
            len=size(mni,2);%the number of the voxels
            dim=size(vol);
            volIndexes=CommonMethods.mni2ind_cols(mni, mat);%3xlen 3D indexes
            volInd1=sub2ind(dim,volIndexes(1,:)',volIndexes(2,:)',volIndexes(3,:)');
            %the one D index
            clusterIndexes=spm_clusters(volIndexes);
            %the index indicating the cluster of the voxels belong to
            
            obj.vClusterIndex=-1*ones(dim);
            
            obj.vMorphDist=ones(dim);
            
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
            
            [~, I]=sort(sizes0);
            nonEmpty=find(~sizes0==0);
            
            lenc=length(nonEmpty);
            sizes=zeros(1,lenc);
            stdc=zeros(1,lenc);
            center=zeros(3,lenc);
            meanc=zeros(1,lenc);
            maxc=zeros(1,lenc);
            minc=zeros(1,lenc);
            medianc=zeros(1,lenc);
            description=cell(1,lenc);
            indexesc=cell(1,lenc);
            mnic=cell(1,lenc);
            
            for i=1:lenc
                i0=I(lenc0-i+1);
                ind=nonEmpty(i0);
                indexes=indexes0{ind};
                indexesc{i}=indexes;
                
                %voxels in the cluster
                volIndt=volInd1(indexes);
                intensities=vol(volIndt);
                sizes(i)=sizes0(ind);
                obj.vClusterIndex(volIndt)=i;
                [c, md]=CommonMethods.getMorphCenter(volIndt,dim);
                %c is the morphological center.
                %md is the morphologicla distance of each voxel pointed by
                %volIndt.
                obj.vMorphDist(volIndt)=md;
                
                cmni=mat*[c;1];
                center(:,i)=cmni(1:3);
                meanc(i)=mean(intensities);
                stdc(i)=std(intensities);
                maxc(i)=max(intensities);
                minc(i)=min(intensities);
                medianc(i)=median(intensities);
                mnic{i}=arrayfun(@(x) [mni(1,x); mni(2,x); mni(3,x)], indexes, 'UniformOutput', false);
                mnic{i}=cell2mat(mnic{i});
                description{i}=[sprintf('%4.0f',i) ': ' CommonMethods.mni2str(center(:,i)) ' size: ' sprintf('%4.0f', sizes(i)) ' mean: ' sprintf('%10.4f', meanc(i))];
            end
            
            obj.stats.mnis=mnic;
            obj.stats.num=len;
            obj.stats.sizes=sizes;
            obj.stats.centers=center;
            obj.stats.means=meanc;
            obj.stats.median=medianc;
            obj.stats.numVoxels=sum(sizes);
            obj.stats.indexes=indexesc;
            obj.stats.descriptions=description;
            obj.stats.std=stdc;
        end
    end
    
end

