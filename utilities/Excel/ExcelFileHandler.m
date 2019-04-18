classdef ExcelFileHandler <handle
    %EXCELFILEHANDLER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fname;
        eF;%the file
        eS;%the active sheet
        e;%the active server
    end
    
    methods
        function obj=ExcelFileHandler(fName)
            obj.fname=fName;
             obj.e=actxserver('excel.application');
             eW=obj.e.Workbooks;
             obj.eF=eW.Open(obj.fname); % your filename here
             obj.eS=obj.eF.ActiveSheet;
%              cols=eS.UsedRange.Columns.Count;
%              rows=eS.UsedRange.Rows.Count;
%              r=eS.UsedRange.Address;
%              range=eS.Range(r);
%              range.EntireColumn.AutoFit;
%              eF.Save;
%              eF.Close; % close the file
%              e.Quit; % close Excel entirely
        end
        function quit(obj)
            obj.e.Quit;
        end
        function save(obj)
            obj.eF.Save;
        end
        function close(obj)
            obj.eF.Close;
        end
        
        function autofitColumWidths(obj)
            r=obj.eS.UsedRange.Address;
            obj.eS.Range(r).EntireColumn.AutoFit;
        end
        function setHorizontalAlignment(obj, key)
            r=obj.eS.UsedRange.Address;
            obj.eS.Range(r).EntireColumn.HorizontalAlignment=key;
        end
        function cell=getCell(obj,r,c)
            %this function returns the cell at row r and col c of the
            %active sheet
            cell=get(obj.eS.cells, 'item', r,c);
        end
        function setCellInteriorColor(obj,rcs,color)
            color=obj.m2vbColor(color);
            len=size(rcs,1);
            for i=1:len
                clabel=obj.getCellLabel(rcs(i,1),rcs(i,2));
                obj.eS.Range(clabel).Interior.Color =color;% Set the color of cell "A1" of Sheet 1 to RED
                %               WB.Worksheets.Item(1).Range('A1').Interior.ColorIndex = 3;
            end
        end
        function label=getCellLabel(obj,r,c)
            label=[obj.getExcelColumnNames(c) num2str(r)];
        end
        function setBorders(obj,rows,cols,lineStyle)
            for r=1:length(rows)
                for c=1:length(cols)
                    cell=get(obj.eS.cells,'item',rows(r),cols(c));
                    hB=cell.borders;
                    hBL = get(hB,'item',1); %left border
                    hBR = get(hB,'item',2); %right
                    hBT = get(hB,'item',3); %top border
                    hBB = get(hB,'item',4); %top border
                    
                    set(hBL,'linestyle',lineStyle);
                    set(hBR,'linestyle',lineStyle);
                    set(hBT,'linestyle',lineStyle);
                    set(hBB,'linestyle',lineStyle);
                end
            end
        end
        function name=getExcelColumnNames(obj,c)
            str='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
            label='0123456789ABCDEFGHIJKLMNOP';
            
            numc=dec2base(c-1,26);
            len=length(numc);
            name=char(len);
            for i=1:len
                name(i)=str(strfind(label,numc(i)));
            end
        end
        function color = m2vbColor(obj,color)
            try
                % Convert color names to RBG triple (0-1) if not already in that format
                if ischar(color)
                    switch lower(color)
                        case {'y','yellow'}, color = [1,1,0];
                        case {'m','magenta'}, color = [1,0,1];
                        case {'c','cyan'}, color = [0,1,1];
                        case {'r','red'}, color = [1,0,0];
                        case {'g','green'}, color = [0,1,0];
                        case {'b','blue'}, color = [0,0,1];
                        case {'w','white',''}, color = [1,1,1]; % empty '' also sets white color
                        case {'k','black'}, color = [0,0,0];
                        otherwise, error(['Invalid color specified: ' color]);
                    end
                elseif ~isnumeric(color) | length(color)~=3 %#ok ML6
                    error(['Invalid color specified: ' color]);
                end
                
                % Convert to Microsoft decimal RGB format
                color = sum(floor(color*255) .* (256.^[0,1,2]));
            catch
                error(['Invalid color specified: ' lasterr]);
            end
        end
        
    end
    
end

