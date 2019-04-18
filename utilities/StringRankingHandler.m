classdef StringRankingHandler <handle
    %STRINGRANKINGHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name;
        useAlias;
        theList;
        theAlias;
        %a cell array store all string members.
        thePositions;
        %the i-th element of thePositions is the index pointing to the cell
        %that holds the string whose ranking is i; The newest stored string
        %element has the highest ranking. 
        maxNEL;%maximum number of element to store
        delimiter;
    end
    
    methods
        function obj=StringRankingHandler(name,maxE,useAlias)
            if(~exist('maxE','var'))
                maxE=Inf;
            end
            if(~exist('useAlias','var'))
                useAlias=false;
            end
            obj.name=name;
            obj.maxNEL=maxE;
            obj.theList={};
            if(useAlias)
                obj.theAlias={};
                obj.useAlias=true;
            end
            obj.thePositions=[];
            obj.name=obj.name;
        end
        
        function store(obj,str,alias)
            idx=find(ismember(obj.theList, str));
            if(isempty(idx))
                obj.add(str,alias);
            else
                if(obj.useAlias)
                    obj.theAlias{idx(1)}=alias;
                end
                obj.promote(idx(1));
            end
        end
        
        function add(obj, str,alias)
            nel=length(obj.theList);
            if(length(obj.theList)<obj.maxNEL)
                pos=nel+1;
                obj.theList{pos}=str;
                if(exist('alias','var'))
                    obj.theAlias{pos}=alias;
                end
                obj.thePositions(pos)=nel+1;
            else
                pos=obj.thePositions(1);
                obj.theList{pos}=str;
                if(exist('alias','var'))
                    obj.theAlias{pos}=alias;
                end
                obj.thePositions(1:nel-1)=obj.thePositions(2:nel);
                obj.thePositions(nel)=pos;
            end
        end
        
        function promote(obj,pos)
            idx=find(obj.thePositions==pos);
            id=idx(1);
            nel=length(obj.thePositions);
            pos=obj.thePositions(id);
            obj.thePositions(id:nel-1)=obj.thePositions(id+1:nel);
            obj.thePositions(nel)=pos;
        end
        
        function [els, alias]=getElements(obj)
            %the the newerly stored elements first (lower index)
            nel=length(obj.thePositions);
            positions=obj.thePositions(nel:-1:1);
            els=obj.theList(positions);
            alias=[];
            if(obj.useAlias)
                alias=obj.theAlias(positions);
            end
        end
        
        function setMaxEls(obj,me)
            obj.maxNEL=me;
        end
        
        function readElement(obj,line)
            fields=strsplit(StringRankingHandler.getDelimiter(), line);
            len=length(fields);
            if(len<=1)
                return;
            end            
            app=fields{1};
            if(~strcmp(app, obj.name))
                return;
            end
            str=fields{2};
            ua=len>3;
            if(~ua)
                obj.store(str);
            else
                obj.store(str,fields{4});
            end
        end
        
        function el = getElement_alias(obj,als)
            el=[];
            if(~obj.useAlias)
                return
            end
            idx=find(ismember(obj.theAlias,als));
            if(length(idx)~=1)
                return;
            end
            el=obj.theList{idx(1)};
        end
        
        function name = getAlias_name(obj,name)
            name=[];
            if(~obj.useAlias)
                return
            end
            idx=find(ismember(obj.theList, name));
            if(length(idx)~=1)
                return;
            end
            name=obj.theAlias{idx(1)};
        end
        
        function strs = toString(obj)
            [els, alias]=obj.getElements();
            ua=obj.useAlias;
            
            len=length(els);
            strs=cell(1,len);
            app=obj.name;
            mn=num2str(obj.maxNEL);
            del=StringRankingHandler.getDelimiter();
            
            for i=1:len
                if(~ua) 
                    line=[app del els{i} del mn];
                else
                    line=[app del els{i} del mn del alias{i}];
                end
                strs{i}=line;
            end
        end        
    end
    
    methods (Static=true)
        function del=getDelimiter()
            del='!@#$Srh1.0%^&*()';
        end
        function name=getAppName(line)
            name=[];
            fields=strsplit(StringRankingHandler.getDelimiter(), line);
            if(length(fields)>1)
                name=fields{1};
            end
        end
        function obj = initiateSRH(line)
            obj=[];
            fields=strsplit(StringRankingHandler.getDelimiter(), line);
            len=length(fields);
            if(len<=1)
                return;
            end            
            app=fields{1};
            mn=str2num(fields{3});
            ua=len>3;
            obj=StringRankingHandler(app,mn,ua);
        end
    end    
end

