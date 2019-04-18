classdef VolumeLandscapeNodeHandler < handle
    %VOLUMELANDSCAPENODEHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lsnodes;%a cell array of volume landscape nodes
        refNode;
    end
    
    methods
        function obj=VolumeLandscapeNodeHandler(imgs, mask, refImg)
            %fnames: the image file names
            len=length(imgs);
            lsnodes=cell(1,len);
            for i=1:len
                cprintf('blue','%s\n',['loading lsv for img: ' imgs{i}]);
                lsnodes{i}=VolumeLandscapeNode(imgs{i},mask);
            end
            obj.refNode=VolumeLandscapeNode(refImg,mask);
            obj.lsnodes=lsnodes;
        end
        function s=calRegionSimilarity(obj,vol1,rnode1,vol2, rnode2)
            %computing the similarity between the two region node
            
            ws=50;
            inds=[1:ws (ws+2):(2*ws+1)];
            profile=VolumeLandscapeNodeHandler.getOrthprofile(vol1,rnode1.peakind,ws);
            v1t=[profile(1,:) profile(2,inds) profile(3,inds)]/profile(1,ws+1);
            profile=VolumeLandscapeNodeHandler.getOrthprofile(vol2,rnode2.peakind,ws);
            v2t=[profile(1,:) profile(2,inds) profile(3,inds)]/profile(1,ws+1);
            
            selection=(~isnan(v1t).*~isnan(v2t))>0;
            v1t=v1t(selection);
            v2t=v2t(selection);
            if(length(v1t) < 3)
                s=0;
            else
                [~,p]=corr(v1t', v2t');
                s=-log(p);
            end
        end
        
        
        function s=calRegionSimilarity_o(obj,rnode1, rnode2)
            %computing the similarity between the two region node
            %It is the average z of the correlation between the profiles of
            %the two regions
            v1t=[];
            v2t=[];
            for i=1:3
                v1=rnode1.orthProfileNode.value{i};%the profile along i-th dimension
                ind1=rnode1.orthProfileNode.center(i);%the cener(peak)of the profile
                v2=rnode2.orthProfileNode.value{i};
                ind2=rnode2.orthProfileNode.center(i);
                v1=v1/v1(ind1);%normalized profile
                v2=v2/v2(ind2);
                ln=min(ind1,ind2);
                rn=min(length(v1)-ind1,length(v2)-ind2);
%                 v1t=[v1t v1((1+(ind1-ln)):(ind1+rn))];%the two profiles are now are the same length with the peaks are at the same location.
%                 v2t=[v2t v2((1+(ind2-ln)):(ind2+rn))];
                v1t=[v1t v1((1+(ind1-ln)):(ind1-1))];%the two profiles are now are the same length with the peaks are at the same location.
                v2t=[v2t v2((1+(ind2-ln)):(ind2-1))];
                v1t=[v1t v1((1+(ind1+1)):(ind1+rn))];%the two profiles are now are the same length with the peaks are at the same location.
                v2t=[v2t v2((1+(ind2+1)):(ind2+rn))];
            end
            if(length(v1t) < 3)
                s=0;
            else
                [~,p]=corr(v1t', v2t');
                s=-log(p);
            end
        end
        
        function [rnode, rnodes, so]=chooseRegionNode(obj,vol, rnodes,volref, refnode)            
            ss=cellfun(@(x)obj.calRegionSimilarity(volref, refnode,vol, x),rnodes);
            [so,inds]=sort(ss,'descend');
            rnodes=rnodes(inds);
            rnode=rnodes{1};
        end
        
        function v=getRegionPeakValue(obj,rnode)
            vec=rnode.orthProfileNode.value{1};
            ind=rnode.orthProfileNode.center(1);
            v=vec(ind);
        end
        
        function exportCorrespodingLXVolumes(obj)
            %this function export aligned peak voxel images for all images
            len=length(obj.lsnodes);            
            refNode=obj.refNode;
            volref=refNode.st.vols;
            refname=refNode.st.fname;
            [~,refname,ext]=fileparts(refname);
            refVol_lxind=refNode.st.vol_lxind;
            refRegions=refNode.st.regionHandler.rnodes;%region nodes in the reference image
            rlen=length(refRegions);
            refRegionVolIndexes=cellfun(@(x)find(refVol_lxind==x.peakind), refRegions, 'UniformOutput', false); %voxel indexes of all voxels belong to each region
            for i=1:len                
                vlsNode=obj.lsnodes{i};
                lsh=vlsNode.st.regionHandler;
                V=spm_vol(vlsNode.st.fname);
                fname=V.fname;
                msg=['exporting the lx volume the file: ' fname];
                cprintf('blue', '%s\n', msg);
                [folder,name,ext]=fileparts(fname);
                vol=vlsNode.st.vols;
                vols=zeros(V.dim);
                vols_lxind=zeros(V.dim);
                for r=1:rlen
                    refnode=refRegions{r};
                    pind=vlsNode.st.vol_lxind(refnode.peakind);                    
                    if(pind==0)
                        continue;
                    end
                    rind=lsh.rinds(pind);
                    rnode=lsh.rnodes{rind};
                    rnodes=arrayfun(@(x)lsh.rnodes{x},rnode.neibors,'UniformOutput',false);
                    [rnode,rnodes,so]=obj.chooseRegionNode(vol,rnodes,volref,refnode);
%                    testing=mod(r,20)==1;
                    testing=false;
                    if(testing)
                        inds=cellfun(@(x)x.peakind,rnodes,'UniformOutput',false);
                        inds=horzcat({refnode.peakind refnode.peakind}, inds);
                        volst=cell(1,length(inds));
                        volst(2:length(volst))={vlsNode.st.vols};
                        volst{1}=refNode.st.vols;
                        vv=VolumeVisualizer(volst(1:9),inds(1:9));
                    end
                    v=obj.getRegionPeakValue(rnode);
                    vols(refRegionVolIndexes{r})=v;
                    vols_lxind(refRegionVolIndexes{r})=rnode.peakind;
                end
                V.fname=fullfile(folder,[name '_' refname '_lx' ext]);
                spm_write_vol(V,vols);
                V.fname=fullfile(folder,[name '_' refname '_lxind' ext]);
                spm_write_vol(V,vols_lxind);
            end
        end
    end  
    
    methods (Static=true)
        function test()
            clc;
            info=COH1_DTI_AnalysisInfo();
            imgs=info.Images.FA;
            refImg=fullfile(info.groupdir,'compImgs','meanFA4D.nii');
            mask=info.mask_spmICV;
            obj=VolumeLandscapeNodeHandler(imgs, mask, refImg);
            obj.exportCorrespodingLXVolumes;
        end
        
        function profile=getOrthprofile(vol,pos,ws)
            %this function compute orthorgonal profile centered by pos
            dim=size(vol);
            if(length(pos)==1)
                [x, y, z]=ind2sub(dim,pos);
                pos=[x y z];
            end
            lent=2*ws+1;
            profile=nan(3,lent);
            poss=zeros(3,lent);
            
            poss(1,:)=(pos(1)-ws):(pos(1)+ws);
            poss(2,:)=pos(2);
            poss(3,:)=pos(3);
            selection=min(((poss>0).*(poss<=repmat(dim',[1,size(poss,2)]))))>0;
            for i=1:length(selection)
                if(selection(i))
                    profile(1,i)=vol(poss(1,i),poss(2,i),poss(3,i));
                end
            end
            
            poss(1,:)=pos(1);
            poss(2,:)=(pos(2)-ws):(pos(2)+ws);
            poss(3,:)=pos(3);
            selection=min(((poss>0).*(poss<=repmat(dim',[1,size(poss,2)]))))>0;
            for i=1:length(selection)
                if(selection(i))
                    profile(2,i)=vol(poss(1,i),poss(2,i),poss(3,i));
                end
            end
            
            poss(1,:)=pos(1);
            poss(2,:)=pos(2);
            poss(3,:)=(pos(3)-ws):(pos(3)+ws);
            selection=min(((poss>0).*(poss<=repmat(dim',[1,size(poss,2)]))))>0;
            for i=1:length(selection)
                if(selection(i))
                    try
                        profile(3,i)=vol(poss(1,i),poss(2,i),poss(3,i));
                    catch err
                        a=1;
                    end
                end
            end
        end
        function displayCorrespondingVoxels(coor,imgInds)
            persistent st;
            
            if(isempty(st))
                info=COH1_DTI_AnalysisInfo();
                imgs=info.Images.FA;
%                imgs=imgs(1:8);
                refImg=fullfile(info.groupdir,'compImgs','meanFA4D.nii');
                mask=info.mask_spmICV;

                obj=VolumeLandscapeNodeHandler(imgs, mask, refImg);
                st.obj=obj;
                [refFolder,refName,ext]=fileparts(refImg);
                st.refImg=obj.refNode.st.fname;
                st.refImg_lxind=obj.refNode.st.fname_lxind;
                st.refV=obj.refNode.st.V;
                st.refVol=obj.refNode.st.vols;
                st.refVol_lx=obj.refNode.st.vol_lx;
                st.refVol_lxind=obj.refNode.st.vol_lxind;
                len=length(imgs);
                st.imgs=imgs;
                st.imgs_lxind=cell(1,len);
                st.imgs_reflxind=cell(1,len);
                st.vols=cell(1,len);
                st.vols_lxind=cell(1,len);
                st.vols_reflxind=cell(1,len);
                for i=1:len
                    img=imgs{i};
                    lsnode=obj.lsnodes{i};
                    cprintf('blue','%s\n',['reading image: ' img]);
                    [folder,name,ext]=fileparts(img);
                    st.imgs_lxind{i}=fullfile(folder,[name '_lxind' ext]);
                    st.imgs_reflx{i}=fullfile(folder,[name '_' refName '_lx' ext]);
                    st.imgs_reflxind{i}=fullfile(folder,[name '_' refName '_lxind' ext]);
                    st.vols{i}=lsnode.st.vols;
                    st.vols_lxind{i}=lsnode.st.vol_lxind;
                    st.vols_reflx{i}=spm_read_vols(spm_vol(st.imgs_reflx{i}));
                    st.vols_reflxind{i}=spm_read_vols(spm_vol(st.imgs_reflxind{i}));
                end
                st.vols=horzcat({st.refVol}, st.vols);
                st.vols_lxind=horzcat({st.refVol_lxind}, st.vols_lxind);
                st.vols_reflx=horzcat({st.refVol_lx}, st.vols_reflx);
                st.vols_reflxind=horzcat({st.refVol_lxind}, st.vols_reflxind);
            end
            if(~exist('coor','var'))
                coor=[-31 30 -37];
            end
            if(~exist('imgInds','var'))
                imgInds=1:9;
            end
            if(size(coor,1)==3)
                coor=coor';
            end
            inds=round(CommonMethods.coor2ind_mat(coor, st.refV.mat));
            ind=sub2ind(st.refV.dim,inds(1),inds(2),inds(3));
            inds_lx=cellfun(@(x)round(x(ind)),st.vols_lxind,'UniformOutput',false);
            inds_reflx=cellfun(@(x)round(x(ind)),st.vols_reflxind,'UniformOutput',false);
            [x,y,z]=ind2sub(st.refV.dim,cell2mat(inds_lx));
            xyz_lx=[x',y',z'];
            [x,y,z]=ind2sub(st.refV.dim,cell2mat(inds_reflx));
            xyz_reflx=[x',y',z'];
            len=length(inds_lx);
            inds_o=cell(1,len);
            inds_o(:)=inds_lx(1);
            v1=arrayfun(@(x)st.vols{x}(inds_lx{x}),1:len);
            v2=arrayfun(@(x)st.vols{x}(inds_reflx{x}),1:len);
            v3=arrayfun(@(x)st.vols{x}(inds_o{x}),1:len);
            v4=arrayfun(@(x)st.vols_reflx{x}(inds_o{x}),1:len);
            vv=VolumeVisualizer(st.vols(imgInds), inds_lx(imgInds));
            vv1=VolumeVisualizer(st.vols(imgInds), inds_reflx(imgInds));
            vv2=VolumeVisualizer(st.vols(imgInds), inds_o(imgInds));
            tmp=st.obj.refNode.st.vol_lxind(ind)
        end
    end
end

