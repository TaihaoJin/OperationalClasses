classdef FractalDimension < handle
    %FRACTALDIMENSION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vols;
        info;
        results;
    end
    
    methods
        function obj=FractalDimension(vols,boxSizeOption)
            if(~exist('boxSizeOption','var'))
                boxSizeOption='exp';
            end
            
            if(~exist('vols','var'))
                rc=true;
                if(rc)
                    vols=randcantor(0.7, 2^7, 3);%using random cantor
                else
                    %using a spherical object (r=20)
                    shape.size=20;
                    shape.name='sphere';
                    shape.fname='/home/tjin/MyProjects/MatlabProjects_2017b/External Software/spm12/canonical/avg152T1.nii';
                    shape=CommonMethods.getShape(shape);
                    V=spm_vol(shape.fname);
                    vols=zeros(V.dim);
                    center=[round(V.dim(1)/2), round(V.dim(2)/2), round(V.dim(3)/2)];
                    shape=CommonMethods.updateShapePosition(shape,center);
                    inds=sub2ind(V.dim,shape.curVIndexes(:,1), shape.curVIndexes(:,2), shape.curVIndexes(:,3));
                    vols(inds)=1;
                end
            end
            obj.vols=vols;
            obj.info.boxSizeOption=boxSizeOption;
            obj.results.err='';
            obj.realignVoxels();
            if(~isempty(obj.results.err))
                return;
            end
            obj.calBoxSizes();
            if(~isempty(obj.results.err))
                return;
            end
            obj.adjustVols();
            obj.calFD();
            if(~isempty(obj.results.err))
                return;
            end
            obj.refineFD();
%            obj.displayFD();
        end
        
        function realignVoxels(obj)
            %this function realign (rotate) voxels such as that the sides
            %of minimal volume bounding box will be aligned with the three
            %axis
            %[rotmat,cornerpoints,volume,surface,edgelength] = minboundbox(x,y,z,metric,level);
            inds=find(obj.vols>0);
            [x,y,z]=ind2sub(size(obj.vols),inds);
            obj.info.center=CommonMethods.findCenter([x';y';z'])';
            obj.info.size=nnz(length(inds));
            try
                [mat,cp,v,s,els]=minboundbox(x,y,z,'v',3);
            catch err
                obj.results.err=err;
                return;
            end
            xyz1=round((mat*[x y z]')');%coordinates of the non-zero voxels after roation. Now the bounding boxe should be parallel 
            l=min(xyz1);%low corner
            t=max(xyz1);%high corner
            dim1=t-l+1;
            x2=xyz1(:,1)-l(1)+1;
            y2=xyz1(:,2)-l(2)+1;
            z2=xyz1(:,3)-l(3)+1;
            inds1=sub2ind(dim1,x2,y2,z2);
            vols1=zeros(dim1);
            vols1(inds1)=1;
            vols1=imfill(vols1,'holes');%posible holes created by roundoff erros of the rotation.
            obj.vols=vols1;
        end
        
        function calBoxSizes(obj)
            %option: exp or lin
            option=obj.info.boxSizeOption;
            inds=find(obj.vols>0);
            [x,y,z]=ind2sub(size(obj.vols),inds);
            xyz=[x y z];
            if(size(xyz,1)==1)
                xyz=[xyz;xyz];
            end
            obj.info.range=[min(xyz); max(xyz)];
            obj.info.size=nnz(obj.vols);
            
            d=min(obj.info.range(2,:)-obj.info.range(1,:))+1;
            if(strcmp(option,'exp'))
                rx=round(0.4*d);
                rx=min(rx, 64);
                nx=floor(log2(rx));
                obj.info.r=arrayfun(@(x)2^x,0:nx);
            elseif(strcmp(option,'lin'))
                rx=round(0.4*d);
                rx=2*floor(rx/2);
                rx=min(rx, 16);
                nx=rx/2;
                obj.info.r=[1 arrayfun(@(x)2*x,1:nx)];
            else
                err('Box size option has to be "exp" or "lin"');
            end
            nx=length(obj.info.r);
            if(nx < 3)
                obj.results.err='the object is too small';
                return;
            end
        end
        
        function adjustVols(obj)
            %readjusting vols to leaving room for offsets and to remove
            %empty porsion of vols to save the computation time. 
            rx=obj.info.r(end);
            dim0=size(obj.vols);
            dim=dim0+rx;
            
            %making rooms for offset
            volst=zeros(dim);
            volst(rx+1:dim(1),rx+1:dim(2),rx+1:dim(3))=obj.vols;
            obj.vols=volst;
            
            nOffsets=5;
            obj.info.offsets=randi([1 rx-1],[3,nOffsets]);
            obj.info.range=obj.info.range+rx;
        end
        
        function calFD(obj)   
            nx=length(obj.info.r);
            if(nx < 3)
                obj.results.err='the object is too small';
                return;
            end
            
            nOffsets=size(obj.info.offsets,2);
            len=length(obj.info.r);
            obj.results.r=zeros(1,len);
            obj.results.nr=zeros(nOffsets,len);
            obj.results.counts=zeros(nOffsets,len);
            obj.results.Vpbs=cell(nOffsets,len);
            obj.results.lacs=zeros(nOffsets,len);
            obj.results.D0=zeros(1,nOffsets);
            obj.results.D1=zeros(1,nOffsets);
            for i=1:len;
                r=obj.info.r(i);
                obj.results.r(i)=r;
                for o=1:nOffsets
                    offset=obj.info.offsets(:,o);
                    [I,counts,vpb]=obj.boxcounts(r,offset);
                    obj.results.nr(o,i)=I;
                    obj.results.counts(o,i)=counts;
                    obj.results.Vpbs{o,i}=vpb;
                    lac=std(vpb)/mean(vpb);
                    obj.results.lacs(o,i)=lac*lac;
                end
            end
            obj.results.Lambda=mean(obj.results.lacs(:,3));
%            obj.results.Lambda=mean(mean(obj.results.lacs));
            for o=1:nOffsets
                x=-arrayfun(@(r)log(r),obj.results.r)';
                y=arrayfun(@(I)log(I),obj.results.nr(o,:))';
                y2=arrayfun(@(I)log(I),obj.results.counts(o,:))';
                mdl=fitlm(x,y);
                mdl2=fitlm(x,y2);
                obj.results.D0(o)=mdl2.Coefficients(2,1).Estimate;
                obj.results.D1(o)=mdl.Coefficients(2,1).Estimate;
            end
        end
        function refineFD(obj)
            if(~isempty(obj.results.err))
                return;
            end
            nOffsets=size(obj.info.offsets,2);
            len=length(obj.info.r);  
            idx=1:len;
            if(len>3)
                idx=idx(1:len-1);
            end
            
            len=length(idx);  
            idx=1:len;
            if(len>3)
                idx=idx(2:len);
            end
            
            len=length(idx);
            obj.info.idx=idx;
            obj.results.r_r=obj.results.r(idx);
            obj.results.nr_r=zeros(nOffsets,len);
            obj.results.counts_r=zeros(nOffsets,len);
            obj.results.D0_r=zeros(1,nOffsets+3);
            obj.results.D1_r=zeros(1,nOffsets+3);
            
            
            for o=1:nOffsets
                x=-arrayfun(@(r)log(r),obj.results.r(idx))';
                y=arrayfun(@(I)log(I),obj.results.nr(o,(idx)))';
                y2=arrayfun(@(I)log(I),obj.results.counts(o,(idx)))';
                mdl=fitlm(x,y);
                mdl2=fitlm(x,y2);
                obj.results.nr_r(o,:)=obj.results.nr(o,idx);
                obj.results.counts_r(o,:)=obj.results.counts(o,idx);
                obj.results.D0_r(o)=mdl2.Coefficients(2,1).Estimate;
                obj.results.D1_r(o)=mdl.Coefficients(2,1).Estimate;
            end
            x=-arrayfun(@(r)log(r),obj.results.r(idx))';
            y=arrayfun(@(x)log(min(obj.results.nr(:,x))),idx)';
            y2=arrayfun(@(x)log(min(obj.results.counts(:,x))),idx)';
            mdl=fitlm(x,y);
            mdl2=fitlm(x,y2);
            obj.results.nr_r(nOffsets+1,:)=min(obj.results.nr(:,idx));
            obj.results.counts_r(nOffsets+1,:)=min(obj.results.counts(:,idx));
            obj.results.D0_r(nOffsets+1)=mdl2.Coefficients(2,1).Estimate;
            obj.results.D1_r(nOffsets+1)=mdl.Coefficients(2,1).Estimate;
            
            x=-arrayfun(@(r)log(r),obj.results.r(idx))';
            y=arrayfun(@(x)log(mean(obj.results.nr(:,x))),idx)';
            y2=arrayfun(@(x)log(mean(obj.results.counts(:,x))),idx)';
            mdl=fitlm(x,y);
            mdl2=fitlm(x,y2);
            obj.results.nr_r(nOffsets+2,:)=mean(obj.results.nr(:,idx));
            obj.results.counts_r(nOffsets+2,:)=mean(obj.results.counts(:,idx));
            obj.results.D0_r(nOffsets+2)=mdl2.Coefficients(2,1).Estimate;
            obj.results.D1_r(nOffsets+2)=mdl.Coefficients(2,1).Estimate;
            
            x=-arrayfun(@(r)log(r),obj.results.r(idx))';
            y=arrayfun(@(x)log(max(obj.results.nr(:,x))),idx)';
            y2=arrayfun(@(x)log(max(obj.results.counts(:,x))),idx)';
            mdl=fitlm(x,y);
            mdl2=fitlm(x,y2);
            obj.results.nr_r(nOffsets+3,:)=max(obj.results.nr(:,idx));
            obj.results.counts_r(nOffsets+3,:)=max(obj.results.counts(:,idx));
            obj.results.D0_r(nOffsets+3)=mdl2.Coefficients(2,1).Estimate;
            obj.results.D1_r(nOffsets+3)=mdl.Coefficients(2,1).Estimate;            
        end
        
        function [I,counts,vpb]=boxcounts(obj,r,offset)
            %scanning the object (non-overlapping) and computing entropy
            %(I), box counts (counts), and voxel per box (vpb).
            dim=size(obj.vols);
            num=obj.info.size;%number of non-zero voxels
            I=0;
            counts=0;
            vpb=[];
            for i1=offset(1):r:dim(1)
                i1F=min(i1+r-1,dim(1));
                for i2=offset(2):r:dim(2)
                    i2F=min(i2+r-1,dim(2));
                    for i3=offset(3):r:dim(3)
                        i3F=min(i3+r-1,dim(3));
                        n=nnz(obj.vols(i1:i1F,i2:i2F,i3:i3F));                        
                        if n > 0
                            counts=counts+1;
                            p=n/num;
                            I=I-p*log(p);
                            vpb(end+1)=n;
                        end
                    end
                end
            end           
        end
        function displayFD(obj)
            figure()
            subplot(2,2,1);
            dim=size(obj.results.counts);
            rows=dim(1);
%             x=-arrayfun(@(r)log(r),obj.results.r_r)';
%             y=arrayfun(@(I)log(I),obj.results.counts_r(rows,:))';           
            x=-arrayfun(@(r)log(r),obj.results.r)';
            y=arrayfun(@(I)log(I),obj.results.counts(rows,:))';           
            scatter(x,y);
            lsline()
            title(['D0: ' num2str(obj.results.D0_r(rows))])
            xlabel('ln(-r)');
            ylabel('ln(N(r))');
            
            subplot(2,2,2);
            x=-arrayfun(@(r)log(r),obj.results.r)';
            y=arrayfun(@(c)mean(obj.results.lacs(:,c)),1:size(obj.results.lacs,2))';           
            scatter(x,y);
            lsline()
            title(['Lambda: ' num2str(obj.results.Lambda)])
            xlabel('ln(-r)');
            ylabel('mean Lambda');
            
            subplot(2,2,3);
            len=size(obj.info.offsets,2)+3;
            x=1:len';
            y=obj.results.D0_r';           
            scatter(x,y);
            title(['D0: ' CommonMethods.getMeanSemStr(y,3,3)])
            xlabel('offset index');
            ylabel('D0');
            
            subplot(2,2,4);
            len=size(obj.info.offsets,2);
            y=obj.results.lacs;           
            x=(1:len)';
            hold on;
            for i=1:size(y,2)
                scatter(x,y(:,i));
            end
            title(['Lambda: ' CommonMethods.getMeanSemStr(reshape(y,[1 numel(y)]),3,3)])
            xlabel('offset index');
            ylabel('lambda');
        end
        function [center, size, D0, Lambda]=getFD(obj)
            center=obj.info.center;
            size=obj.info.size;
            if isempty(obj.results.err)
                D0=obj.results.D0_r(end);
                Lambda=obj.results.Lambda;
            else
                D0=nan;
                Lambda=nan;
            end
        end
    end
    methods (Static=true)        
        function main()
%            FractalDimension.buildTestObjs();
        end
        function [centers, sizes, FDs, Lambdas]=computeFD(vols,boxSizeOption)
            %this function compute FD for each cluster in vols
            if(~exist('boxSizeOption','var'))
                boxSizeOption='exp';
            end
            idx=find(vols>0);
            dim=size(vols);
            [ix,iy,iz]=ind2sub(dim,idx);
            coors=[ix iy iz];
            cluster_coors=CommonMethods.findClusters(coors');
            num=length(cluster_coors);
            sizes=zeros(1,num);
            FDs=zeros(1,num);
            centers=cell(1,num);
            Lambdas=zeros(1,num);
            for i=1:num
                cur_vols=zeros(dim);
                ixyz=cluster_coors{i}';
                inds=sub2ind(dim,ixyz(:,1),ixyz(:,2),ixyz(:,3));
                cur_vols(inds)=1;
                obj=FractalDimension(cur_vols,boxSizeOption); 
                [center,v,FD,Lambda]=obj.getFD();
                centers{i}=center;     
                sizes(i)=v;
                FDs(i)=FD;
                Lambdas(i)=Lambda;
            end
        end
        function buildTestObjs()
            num=50;
            objs=cell(1,num);
            for r=1:num
                objs{r}=FractalDimension();
            end
            fname='/home3/data/Salgia_Fractal/group/objs_3D_RandomCantorSet.mat';
            save(fname,'objs');
        end
        function examTestObjs()
            fname='/home3/data/Salgia_Fractal/group/objs_3D_RandomCantorSet.mat';
            load(fname);
            runs=length(objs);
            for i=1:length(objs)
                obj=objs{i};
                obj.refineFD();
                objs{i}=obj;
            end
            
            obj=objs{1};
            figure()
            subplot(2,2,1);
            dim=size(obj.results.counts);
            rows=dim(1);
            x=-arrayfun(@(r)log(r),obj.results.r_r)';
            y=arrayfun(@(I)log(I),obj.results.counts_r(rows,:))';           
            scatter(x,y);
            lsline()
            title(['FD: ' num2str(obj.results.D0_r(rows))])
            xlabel('ln(1/r)');
            ylabel('ln(N(r))');
            
            subplot(2,2,2);
            x=-arrayfun(@(r)log(r),obj.results.r_r)';
            y=arrayfun(@(I)log(I),obj.results.nr_r(rows,:))';           
            scatter(x,y);
            lsline()
            title(['D1: ' num2str(obj.results.D1_r(rows))])
            xlabel('ln(1/r)');
            ylabel('ln(I(r))');
            
            subplot(2,2,3);
            y=cellfun(@(x)x.results.D0_r(end),objs);
            len=length(y);
            x=(1:len)';
            scatter(x,y);
            title(['FD: ' CommonMethods.getMeanSemStr(y,3,3)])
            xlabel('Run');
            ylabel('FD');
            
            subplot(2,2,4);
            y=cellfun(@(x)x.results.D1_r(end),objs);
            len=length(y);
            x=(1:len)';
            scatter(x,y);
            title(['D1: ' CommonMethods.getMeanSemStr(y,3,3)])
            xlabel('Run');
            ylabel('D1');
        end
    end
end

