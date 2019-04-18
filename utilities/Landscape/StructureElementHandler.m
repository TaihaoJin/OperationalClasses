classdef StructureElementHandler
    %STRUCTUREELEMENTHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim; %the dimention of volume in which the structure element would be used
        strel;%n-dimensional index relative to the current point
        strel1;% one dimensional index
        indShifts; %shift of 1d index
        numel;%number of elements in strel
    end
    
    methods
        function obj=StructureElementHandler(dim)
            obj.dim=dim;
            
            m=length(dim);
            v0=[-1 0 1]';
            for i=2:m
                v1=[-1 0 1];
                len0=size(v0,1);
                len1=length(v1);
                v=zeros(len0*len1,i);
                num=0;
                for i0=1:len0                    
                    for i1=1:len1
                        num=num+1;
                        v(num,:)=[v0(i0,:) v1(i1)];
                    end
                end    
                v0=v;
            end
            obj.numel=num;
            obj.strel=v;
            
            len=size(v,1);            
            v1=zeros(1,len);
            pos0=2*ones(1,m);
            
            
            line='dimt=[3';
            for i=2:m
                line=[line ' 3'];
            end
            line=[line '];'];
            eval(line);
            
            indshifts=zeros(dimt);
            for i=1:len
                pos1=pos0+v(i,:);
                line0=['sub2ind(dim, ' num2str(pos0(1))];
                line1=['sub2ind(dim, ' num2str(pos1(1))];
                line2=['indshifts(' num2str(pos1(1))];
                for j=2:m
                    line0=[line0 ', ' num2str(pos0(j))];
                    line1=[line1 ', ' num2str(pos1(j))];
                    line2=[line2 ', ' num2str(pos1(j))];
                end
                line0=[line0 ')'];
                line1=[line1 ')'];
                line2=[line2 ')=v1(i);'];
                v1(i)=eval(line1)-eval(line0);
                eval(line2);               
            end           
            obj.strel1=v1';
            obj.indShifts=indshifts;
        end
        
        function neighbors=getNeighbors(obj,pos)
            %pos is a scalar or a 1 by n vector, n==length(obj.dim)
            cols=length(pos);
            if(cols==1)
                el=obj.strel1;
            else
                el=obj.strel;
            end
            neighbors=repmat(pos,[obj.numel,1])+el;
        end
        
        function ind=getInd(obj, ind,shift0)
            %shift is a 1 by m vector for the shift
            m=length(obj.dim);
            shift=shift0+2*ones(1,m);
            if(m==3)
                ind=ind+obj.indShifts(shift(1),shift(2),shift(3));
            else
                line=['obj.indShifts(' num2str(shift(1))];
                line1='dim=[3';
                for i=2:m
                    line=[line ', ' num2str(shift(i))];
                    line1=[line1 ', 3'];
                end
                line=[line ')'];
                line1=[line1 '];'];
                eval(line1);
                ind=ind+eval(line);
            end           
        end
    end
    
    methods (Static=true)
        function test()
            %testing three dimension
            dim=[255 176 382];
            pos=[27 96 124];
            ind=sub2ind(dim,pos(1),pos(2),pos(3));
            obj=StructureElementHandler(dim);
            nbrs=obj.getNeighbors(pos);
            nbrs1=obj.getNeighbors(ind);
            nbrs3=zeros(27,3);
            nbrs1t=zeros(27,1);
            for i=1:27
                [x,y,z]=ind2sub(dim,nbrs1(i));
                nbrs3(i,:)=[x y z];
                nbrs1t(i)=obj.getInd(ind,obj.strel(i,:));
            end
            tt3=nbrs-nbrs3;
            tt3t=nbrs1-nbrs1t;
            
            %testing four dimension
            dim=[255 176 382 678]; 
            pos=[27 96 124 556];
            ind=sub2ind(dim,pos(1),pos(2),pos(3),pos(4));
            obj=StructureElementHandler(dim);
            nbrs=obj.getNeighbors(pos);
            nbrs1=obj.getNeighbors(ind);
            nbrs4=zeros(81,4);
            nbrs1t=zeros(81,1);
            for i=1:81
                [x,y,z,z1]=ind2sub(dim,nbrs1(i));
                nbrs4(i,:)=[x y z z1];
                nbrs1t(i)=obj.getInd(ind,obj.strel(i,:));
            end
            tt4=nbrs-nbrs4;
            tt4t=nbrs1-nbrs1t;
            %both cases passed. 1/17/2017
        end
    end
end

