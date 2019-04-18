classdef regionNodeHandler
    %REGIONNODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rnodes;% a cell array that holds region nodes
        rinds;% an array that holds the region indexes of all maxima (peaks)
        seh;% the structure element handler
    end
    
    methods        
        function obj=regionNodeHandler(lsNode)
            %lsNode an object of VovlumeLandscapeNode;
            
            peakinds=unique(lsNode.st.vol_lxind);
            peakinds=peakinds(peakinds>0);
            len=length(peakinds);
            obj.rnodes=cell(1,len);
            pindx=max(peakinds);
            obj.rinds=zeros(1,pindx);
            obj.rinds(peakinds)=1:len;
            obj.seh=StructureElementHandler(size(lsNode.st.vol_lxind));
            %building region nodes
            for r=1:len
                rnode.rind=r;%the region index used to identify it. 
                rnode.peakind=peakinds(r);% the peak location of the region
                rnode.neibors=obj.getNeighborRegions(lsNode,peakinds(r));
                rnode=obj.calRegionOrthProfile(lsNode,rnode);
                obj.rnodes{r}=rnode;
            end
        end
        
        function rinds=getNeighborRegions(obj,lsNode,peakInd)
            inds=find(lsNode.st.vol_lxind==peakInd);
            len=length(inds);
            pinds=[];
            for i=1:len
                nbrs=obj.seh.getNeighbors(inds(i));
                pinds=[pinds; nbrs];
            end
            pinds=lsNode.st.vol_lxind(pinds);
            pinds=pinds(pinds>0);
            pinds=unique(pinds);
            rinds=obj.rinds(pinds);
        end
        
        
        function rnode=calRegionOrthProfile(obj,lsNode,rnode)            
            orthProfileNode=cell(1,3);%the three orthoganal profile values
            ws=20;
            uvec=eye(3);
            peak=rnode.peakind;
            trail=[peak];
            for i=1:3
                d=-uvec(i,:);
                pos=obj.seh.getInd(peak,d);
                num=1
                while(num<=ws)
                    pos=obj.seh.getInd(pos,d);
                    trail=[pos trail];
                    num=num+1;
                end
                d=uvec(i,:);
                pos=obj.seh.getInd(peak,d);
                num=1;
                while(num<=ws)
                    trail=[trail pos];
                    pos=obj.seh.getInd(pos,d);
                    num=num+1;
                end
                orthProfileNode{i}=lsNode.st.vols(trail);
                rnode.orthProfile=orthProfile;
            end
        end
        function rnode=calRegionOrthProfile_o(obj,lsNode,rnode)            
            orthProfileNode.trail=cell(1,3);%the three orthoganal profile positions
            orthProfileNode.value=cell(1,3);%the three orthoganal profile values
            orthProfileNode.center=zeros(1,3);%the indexes of the peak in the three profiles
            
            orthProfileNode.peak=[];%the index of the peak in the profile
            uvec=eye(3);
            peak=rnode.peakind;;
            for i=1:3
                d=-uvec(i,:);
                pos=obj.seh.getInd(peak,d);
                trail=[peak];
                while(lsNode.st.vol_lxind(pos)==peak)
                    trail=[pos trail];
                    pos=obj.seh.getInd(pos,d);
                end
                orthProfileNode.center(i)=length(trail);
                d=uvec(i,:);
                pos=obj.seh.getInd(peak,d);
                while(lsNode.st.vol_lxind(pos)==peak)
                    trail=[trail pos];
                    pos=obj.seh.getInd(pos,d);
                end
                orthProfileNode.trail{i}=trail;
                orthProfileNode.value{i}=lsNode.st.vols(trail);
                rnode.orthProfileNode=orthProfileNode;
            end
        end
    end   
    
    methods (Static=true)
        function test()
            info=COH1_DTI_AnalysisInfo();
            imgs=info.Images.FA;
            mask=info.mask_spmICV;
            fname=imgs{1};
            lsNode=VolumeLandscapeNode(fname,mask);
            obj=regionNodeHandler(lsNode);
            a=1;
        end
    end
end

