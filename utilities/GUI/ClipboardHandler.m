classdef ClipboardHandler
    %CLIPBOARDHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
    end
    methods (Static=true)
        function copyCellarray(data)
            tab=sprintf('\t');
            nl=sprintf('\n');
            len=size(data,1);
            line=CommonMethods.cell2str(data(1,:),tab);
            for r=2:len
                line=[line nl CommonMethods.cell2str(data(r,:),tab)];
            end
            clipboard('copy',line);
        end
    end
end

