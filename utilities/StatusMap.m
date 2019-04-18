classdef StatusMap <handle
    %STATUSMAP Summary of this class goes here
    %   Detailed explanation goes here
    % This class is for saving the status of an app
    properties
        appID;%name of the app
        statusFile;
        theMap;
    end
    
    methods
        function obj=StatusMap(app)
            homeD=CommonMethods.getHome();
            appDir=fullfile(homeD,'appData');
            if(~exist(appDir,'file'))
                mkdir(appDir);
            end
            obj.appID=app;
            obj.statusFile=fullfile(appDir,[app '.ini']);
            obj.theMap=containers.Map();
            if(exist(obj.statusFile,'file'))
                obj.load();
            end
        end
        function save(obj)
            fname=obj.statusFile;
            fid=fopen(fname,'wt');
            if(fid<=0)
                return;
            end
            keySet=keys(obj.theMap);
            for k=1:length(keySet)
                key=keySet{k};
                RSH=obj.theMap(key);
                strs=RSH.toString();
                %newer element first (lower index)
                len=length(strs);
                for l=1:len
                    id=len-l+1;
                    %storing the newest last;
                    fprintf(fid,'%s\n',strs{id});
                end
            end
            fclose(fid);
        end
        function tf=hasKey(obj,key)
            tf=isKey(obj.theMap, key);
        end
        function srh=getSRHandler(obj,key)
            srh=[];
            if(isKey(obj.theMap,key))
                srh=obj.theMap(key);
            end
        end
        function putSRHandler(obj,key,value)
            if(isempty(key))
                key=value.name;
            end
            obj.theMap(key)=value;
        end
        function putStatustString(obj,key, str)
            if(~isKey(obj.hasKey(obj,key)))
                theHandler=StringRankingHandler(key);
            end
            theHandler.store(str);
            obj.theMap(key)=theHandler;
        end
        function clearTheMap(obj)
%            keySet=keys(obj);
%            remove(obj,keySet);
            clear(obj.theMap);
        end
        function load(obj)
            lines=CommonMethods.getTextLines(obj.statusFile);
            len=length(lines);
            for l=1:len
                line=lines{l};
                key=StringRankingHandler.getAppName(line);
                if(isempty(key))
                    continue;
                end
                if(isKey(obj.theMap,key))
                    srh=obj.theMap(key);
                else
                    srh=StringRankingHandler.initiateSRH(line);
                    obj.theMap(key)=srh;
                end
                srh.readElement(line);
            end
        end
    end
    
end

