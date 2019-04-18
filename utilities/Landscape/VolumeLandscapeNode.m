classdef VolumeLandscapeNode < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        st
    end
    
    methods
        function obj=VolumeLandscapeNode(fname,mask)
            [folder,name,ext]=fileparts(fname);
            obj.st.fname=fname;%the main image
            obj.st.mask=mask;%the mask
            [~, mname, ~]=fileparts(mask);
            obj.st.fname_lxind=fullfile(folder,[name '_lxind' ext]);%index of the local maximum associated with voxel;
            obj.st.fname_lx=fullfile(folder,[name '_lx' ext]);%index of the local maximum associated with voxel;
            obj.st.fname_region=fullfile(folder,[name '_' mname '_region.mat']); %region index of each voxel;
            obj.st.version=1.2;%The version number will be matched with that of saved object. It will be rebuilt if they dont's match
%            obj.st.version=1.3;%The version number will be matched with that of saved object. It will be rebuilt if they dont's match
            ok=obj.loadNode();
            if(~ok)
                obj.buildNode();   
                save(obj.st.fname_region, 'obj');
            end
            if(~isfield(obj.st,'regionHandler'))
                obj.st.regionHandler=regionNodeHandler(obj);
                save(obj.st.fname_region, 'obj');
            end
        end
        function ok=loadNode(obj)
            try
                X=load(obj.st.fname_region);
                ok=obj.st.version==X.obj.st.version;
                if(~ok)
                    return;
                end
                obj.st=X.obj.st;
            catch err
                ok=false;
            end
        end
        function buildNode(obj)
            obj.st.V=spm_vol(obj.st.fname);
            obj.st.vols=spm_read_vols(obj.st.V);
            m=CommonMethods.getVoxelValues(obj.st.fname,obj.st.mask,0);
            obj.st.indxm=find(m{1}>0); %voxel indexes under the mask
            tmp=obj.st.vols(obj.st.indxm);%the voxel value under the mask
            x=max(tmp);
            n=min(tmp);
            obj.st.volst=(n-1)*ones(obj.st.V.dim);
            obj.st.volst(obj.st.indxm)=tmp;
            
            %[peak, poss] = findPeak(vol, pos, sign);
            [x, y, z]=ind2sub(obj.st.V.dim,obj.st.indxm);
            xyz=[x y z];
            obj.st.vol_lxind=zeros(obj.st.V.dim);     
            obj.st.vol_lx=zeros(obj.st.V.dim);
            
            sign=1;
            tmp=[97 157 74];
            for i=1:length(obj.st.indxm)
                if(obj.st.vol_lxind(obj.st.indxm(i))>0)
                    %the position has been labeled already
                    continue;
                end
                if(all(~(xyz(i,:)-tmp)))
                    i=i;
                end
                [peak, poss]=CommonMethods.findPeak(obj.st.volst,xyz(i,:),sign);
                poss=[peak;poss];
                inds=sub2ind(obj.st.V.dim,poss(:,1),poss(:,2),poss(:,3));
                ind=inds(1);
                lx=obj.st.vols(ind);
                for j=1:size(poss,1)
                    obj.st.vol_lxind(inds(j))=ind;
                    obj.st.vol_lx(inds(j))=lx;
                end
            end
        end   
        function outputImgs(obj)
            V=obj.st.V;
            V.fname=obj.st.fname_lx;
            spm_write_vol(V,obj.st.vol_lx);
            V.fname=obj.st.fname_lxind;
            V.dt(1)=8;
            spm_write_vol(V,obj.st.vol_lxind);
        end
        function saveNode(obj)
            save(obj.st.fname_region,'obj');
        end
        function dcoor=getCoorShift(obj,coor)
            %this function computes the differences in coordinates between
            %coor and the it's associated peak location
            %coor and dcoors are 1x3 row vectors
            inds=CommonMethods.mni2ind(coor,obj.st.V.mat);
            peak=obj.st.vol_lxind(inds(1),inds(2),inds(3));
            [x,y,z]=ind2sub(obj.st.V.dim,peak);
            xyz=[x y z];
            coor1=CommonMethods.ind2coor_mat(xyz, obj.st.V.mat);
            dcoor=coor1-coor;
        end
    end
    
    methods (Static=true)
        function test()
            info=COH1_DTI_AnalysisInfo();
            imgs=info.Images.FA;
            meanImg=fullfile(info.groupdir,'compImgs','meanFA4D.nii');
%            CommonMethods.spm_check_registration({meanImg});
            vol=spm_read_vols(spm_vol(meanImg));
            dim=size(vol);
            ind=fix(dim/2);
            ws=max(ind);
            CommonMethods.showVolComparison({vol vol vol vol},{ind ind ind ind});
            m(isnan(m))=max(max(m));
            fh=figure();
            imshow(mat2gray(m));
            mask=info.mask_spmICV;
            obj=VolumeLandscapeNode(meanImg,mask);
            obj.outputImgs();
            if(5>3)
                return;
            end
            CommonMethods.spm_check_registration(horzcat({meanImg},imgs(1:3)));
            parfor i=1:length(imgs)
                fname=imgs{i};
                obj=VolumeLandscapeNode(fname,mask);
                obj.outputImgs();
                fprintf(1,'%s\n',['processed: ' fname]);
            end
            %            CommonMethods.spm_check_registration({obj.st.fname obj.st.fname_lx});
        end
    end
end

