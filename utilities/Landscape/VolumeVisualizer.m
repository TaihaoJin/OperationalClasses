classdef VolumeVisualizer < handle
    %VOLUMEVISUALIZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vols;%3D arraies to be visualized, a cell array
        inds;%the volume indeces, a cell array
        par;%visualization parameters
        fh;%the figure handle
        M;%the 2D master matrix composed of the orthview oxel values of the volumes
        ind2s;%the 2D inds corresponding to inds in the master matrix M
    end
    
    methods
        function obj=VolumeVisualizer(vols, inds, par)
            obj.vols=vols;
            if(length(inds{1})==1)%converting the indeces 1D to 3D
                for i=1:length(inds)
                    [x,y,z]=ind2sub(size(vols{i}),inds{i});
                    inds{i}=[x y z];
                end
            end
            obj.inds=inds;
            len=length(inds);
            if(~exist('par','var'))
                par=[];
            end
            if(~isfield(par,'zoom'))
                par.zoom=3;
            end
            if(~isfield(par,'margin_orth'))
                par.margin_orth=5;%the spacing between the three orthoganal views of the same volue
            end
            if(~isfield(par,'margin'))
                par.margin=15;
            end
            if(~isfield(par,'ws'))%number of voxels to display on each side of the voxel pointed by the voxel index
                par.ws=50;
            end
            if(~isfield(par,'w'))
                par.w=floor(sqrt(len));%number of volumes shown in one row
            end
            if(~isfield(par,'mag'))%magnification (how many matrix elements in the master matrix for a single voxel value)
                par.mag=5;
            end
            if(~isfield(par,'h'))%number of volumes shown in one row
                par.h=ceil(sqrt(len));
                if(par.h*par.w<len)
                    par.h=par.h+1;
                end
            end
            obj.par=par;
            obj.composeMasterMatrix();
            obj.showVolComparison();
        end
        
        function composeMasterMatrix(obj)
            %this function cosntructs the master matrix
            %
            lent=2*(2*obj.par.ws+1)+obj.par.margin_orth;%the with/height of a orthview matrix
            W=obj.par.w*(lent+obj.par.margin)-obj.par.margin;
            H=obj.par.h*(lent+obj.par.margin)-obj.par.margin;
            len=length(obj.vols);
            obj.ind2s=cell(1,len);
            M=nan(H,W);
            no=0;
            for i=1:obj.par.h     
                si=(i-1)*(lent+obj.par.margin);
                for j=1:obj.par.w      
                    sj=(j-1)*(lent+obj.par.margin);
                    no=no+1;   
                    if(no>len)
                        continue;
                    end
                    [m,centers]=VolumeVisualizer.getOrthviews(obj.vols{no},obj.inds{no},obj.par.ws,obj.par.margin_orth);
                    dim=size(m);
                    for t1=1:dim(1)
                        for t2=1:dim(2)
                            M(si+t1,sj+t2)=m(t1,t2);
                        end
                    end
                    for c=1:3
                        centers{c}(1)=centers{c}(1)+si;
                        centers{c}(2)=centers{c}(2)+sj;
                    end                    
                    obj.ind2s{no}=cellfun(@(x)obj.par.zoom*x,centers,'UniformOutput',false);
                end
            end
            mx=max(max(M));
            mn=min(min(M));
            d=0.01*(mx-mn);
            M(isnan(M))=mx+d;
            obj.M=VolumeVisualizer.magnify(M,obj.par.zoom);
        end
        function showVolComparison(obj)
            dim=size(obj.M);
            obj.fh=figure();
            set(obj.fh, 'position', [15,15,dim(2),dim(1)]);
            imshow(mat2gray(obj.M));
            hold on;
            num=length(obj.vols);
            
            ll=obj.par.ws;%the length of the crosshair
            for i=1:num
                for j=1:3
                    center=obj.ind2s{i}{j};
                    x=[center(1)-ll center(1)+ll];
                    y=[center(2) center(2)];
                    line(x,y,'Color','blue');
                    x=[center(1) center(1)];
                    y=[center(2)-ll center(2)+ll];
                    line(x,y,'Color','blue');
                end
            end
        end
    end
    
    methods (Static=true)
        function [m,centers]=getOrthviews(vol,ind,ws,margin)
            %this function get an orthogonal views of a 3D volume at ind
            %ws: the window size
            len=2*ws+1;
            if(~exist('margin','var'))
                margin=5;
            end
            L=2*len+margin;
            
            m=nan(L,L);
            dim=size(vol);
            if(length(ind)==1)
                [x,y,z]=ind2sub(dim,ind);
                ind=[x y z];
            end
            shift=[[0 0 len+margin];[0 len+margin 0]];
            
            moves={[0 1 0; 0 0 1] [1 0 0; 0 0 1] [1 0 0; 0 1 0]};
            centers=cell(1,3);%the positions of the poxel pointed by ind in m
            
            for a=1:3
                mv=moves{a};
                m1=mv(1,:);
                m2=mv(2,:);
                lt=ind-ws*m1-ws*m2;
                si=shift(1,a);
                sj=shift(2,a);
                for i=1:len                    
                    for j=1:len
                        indt=lt+m1*(i-1)+m2*(j-1);
                        if(all((indt>0).*(indt<=dim)))
                            try
                                m(si+i,sj+j)=vol(indt(1),indt(2),indt(3));
                            catch err
                                continue;
                            end
                        end                        
                    end
                end
                n=min(min(m));
                x=max(max(m));
                d=0.001*(x-n);
                m(isnan(m))=n-d;
                m((len+1):(len+margin),1:L)=nan;
                m(1:L,(len+1):(len+margin))=nan;
                
                centers{a}=[si+ws+1,sj+ws+1];
            end            
         end
        function M=magnify(m,nmag)
             dim=size(m);
             Dim=nmag*dim;
             lt=[0 0];
             M=zeros(Dim);
             for i=1:dim(1)
                 lt(1)=nmag*(i-1);
                 for j=1:dim(2)
                     
                     lt(2)=nmag*(j-1);
                     val=m(i,j);
                     M((lt(1)+1):(lt(1)+nmag),(lt(2)+1):(lt(2)+nmag))=val;
                 end
             end
         end
    end
    
end

