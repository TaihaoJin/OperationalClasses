classdef Shape_Handler < handle
    %SHAPE_HANDLER Summary of this class goes here
    %   Detailed explanation goes here
    %%Author: Taihao Jin
    %%Date: 2/27/2015
    %%Purpose: Shape analysis:
    properties
        vol;% the original binary image
        dim;
        se;%structural element
        vol_MorphDist;% the images in which the morpholocial distance are marked.
        computeInner;%compute the morphdist of inner points when computeInner is true
        computeOuter;%compute the morphdist of outer points when computeOuter is true
        Gap=11;
        morphCenter;%the inner point has the smallest morphDist;
    end
    
    methods
        function obj = Shape_Handler(vol0,computeInner, computeOuter, se)
            if(exist('computeInner', 'var'))
                if(numel(computeOuter)~=1)
                    error('woring input for ''computeInner''');
                end
            else
                computeInner=true;
            end
            
            if(exist('computeOuter', 'var'))
                if(numel(computeOuter)~=1)
                    error('woring input for ''computeOuter''');
                end
            else
                computeOuter=true;
            end
            
            if(~exist('vol0','var'))
                vol0=zeros(15,16,14);
                vol0(4:10,4:11,4:10)=1;
            end
            
            if(~exist('stel','var'))
                se=strel(ones(3,3,3));
            end
            
            obj.vol=vol0;
            obj.dim=size(vol0);
            obj.computeInner=computeInner;
            obj.computeOuter=computeOuter;
            obj.se=se;
            obj.calMorphDist();
        end
        
        function calMorphDist(obj,computeInner, computeOuter,se)
            if(exist('computeInner', 'var'))
                if(numel(computeOuter)~=1)
                    error('woring input for ''computeInner''');
                end
            else
                computeInner=obj.computeInner;
            end
            
            if(exist('computeOuter', 'var'))
                if(numel(computeOuter)~=1)
                    error('woring input for ''computeOuter''');
                end
            else
                computeOuter=obj.computeOuter;
            end
            
            if(~exist('se','var'))
                se=obj.se;
            end
            
            obj.computeInner=computeInner;
            obj.computeOuter=computeOuter;
            obj.se=se;
            vol0=obj.vol;
            volm=ones(size(vol0));
            volm(vol0==1)=0;
            vol_dial=imdilate(vol0,se);
            indexes=find(vol_dial-vol0);
            
            if(computeOuter)
                MorphDist=1;
                while(~isempty(indexes))
                    volm(indexes)=MorphDist;
                    vol0=vol_dial;
                    vol_dial=imdilate(vol0,se);
                    indexes=find(vol_dial-vol0);
                    MorphDist=MorphDist+1;
                end
            end
            
            if(computeInner)
                vol0=obj.vol;
                vol_erode=imerode(vol0,se);
                indexes=find(vol0-vol_erode);
                MorphDist=-1;
                while(~isempty(indexes))
                    volm(indexes)=MorphDist;
                    vol0=vol_erode;
                    vol_erode=imerode(vol0,se);
                    indexes=find(vol0-vol_erode);
                    MorphDist=MorphDist-1;
                end
            end
            
            minMD=min(min(min(volm)));
            inds=find(volm==minMD);
            [s1,s2,s3]=ind2sub(obj.dim,inds);
            coors=zeros(3,length(inds));
            coors(1,:)=s1;
            coors(2,:)=s2;
            coors(3,:)=s3;
            obj.morphCenter=CommonMethods.findCenter(coors);
            obj.vol_MorphDist=volm;
        end
        
        function volDiff = calVolDiff(obj, vol)
            volDiff=zeros(size(vol));
            diffIndexes=find(vol-obj.vol);
            volDiff(diffIndexes)=obj.vol_MorphDist(diffIndexes);
        end
        
        function diffStats = calDiffStats(obj, vol)
            volDiff = obj.calVolDiff(vol);
            pIDs= volDiff>0;
            nIDs= volDiff<0;
            
            distsP=volDiff(pIDs);
            distsN=volDiff(nIDs);
            
            diffStats.numP=numel(distsP);
            diffStats.maxP=max(distsP);
            diffStats.meanP=mean(distsP);
            diffStats.medianP=median(distsP);
            
            diffStats.numN=numel(distsN);
            diffStats.minN=min(distsN);
            diffStats.meanN=mean(distsN);
            diffStats.median=median(distsN);
        end
        
        function writeMorphDist(obj,V,fname,vol)
            if(~exist('vol','var'))
                vol=obj.vol_MorphDist;
            end
            V.dt(1)=8;
            if(isfield(V,'pinfo')) V=rmfield(V,'pinfo'); end;
            V.fname=fname;
            vn=vol<0;
            vol=vol-obj.Gap*vn;
            spm_write_vol(V,vol);
        end
        
    end
    
    methods (Static = true)
        function [vol1, vol2] = getDefaultVols()
            vol1=zeros(18,16,17);
            vol1(4:10,4:10,4:10)=1;            
            vol2=zeros(size(vol1));
            vol2(8:14,4:10,4:10)=1;
        end
        
        function [volDiff, vol1, vol2] = volDiffTest()
            [vol1, vol2]=Shape_Handler.getDefaultVols();
            handler=Shape_Handler(vol1);
            volDiff=handler.calVolDiff(vol2);
        end  
        
        function diffStats = diffStatsTest()
            [vol1, vol2]=Shape_Handler.getDefaultVols();
            handler=Shape_Handler(vol1);
            diffStats=handler.calDiffStats(vol2);
        end
        
        function [vec, namesFields] = diffStats2Vec(stats)
            namesFields=filednames(stats);
            vect=cellfun(@(x) stats.(x),namesFields);
        end
    end        
end

