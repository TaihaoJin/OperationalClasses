function varargout = quickTabular_spmSummary(varargin)
% QUICKTABULAR_SPMSUMMARY MATLAB code for quickTabular_spmSummary.fig
%      QUICKTABULAR_SPMSUMMARY, by itself, creates a new QUICKTABULAR_SPMSUMMARY or raises the existing
%      singleton*.
%
%      H = QUICKTABULAR_SPMSUMMARY returns the handle to a new QUICKTABULAR_SPMSUMMARY or the handle to
%      the existing singleton*.
%
%      QUICKTABULAR_SPMSUMMARY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in QUICKTABULAR_SPMSUMMARY.M with the given input arguments.
%
%      QUICKTABULAR_SPMSUMMARY('Property','Value',...) creates a new QUICKTABULAR_SPMSUMMARY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before quickTabular_spmSummary_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to quickTabular_spmSummary_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help quickTabular_spmSummary

% Last Modified by GUIDE v2.5 02-Mar-2018 13:52:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @quickTabular_spmSummary_OpeningFcn, ...
                   'gui_OutputFcn',  @quickTabular_spmSummary_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

h=CommonMethods.getFigureByName('quickTabular_spmSummary');
if(~isempty(varargin))
    a=1;
end
% End initialization code - DO NOT EDIT


% --- Executes just before quickTabular_spmSummary is made visible.
function quickTabular_spmSummary_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to quickTabular_spmSummary (see VARARGIN)

% Choose default command line output for quickTabular_spmSummary
% set(hObject,'WindowButtonDownFcn',@clickcallback);
% set(hObject,'KeyPressFcn',@keypresscallback);
% set(hObject,'KeyReleaseFcn',@keyreleasecallback);
FM.addFigure(hObject,'StatSummary');
handles.output = hObject;
guiData=handles;
guiData.lines={};
guiData.table.header={'A' 'B' 'C' 'D' 'E' 'F'};
guiData.table.data=cell(6,6);
guiData.fname='';
guiData.dim=[1 1];
guiData.CoordinateColumn=1;
guiData.mainGui=hObject;
guiData.rebuildStatSummary=true;
% Update handles structure
guidata(hObject, guiData);
try 
%    displaySPMStatSummary(hObject);
catch err
    cprintf('red','%s\n',err.message);
end
%updateGUI(hObject);
% UIWAIT makes quickTabular_spmSummary wait for user response (see UIRESUME)
% uiwait(handles.figure1);



% --- Executes on button press in pasteFromClipboardTB.
function pasteFromClipboardTB_Callback(hObject, eventdata, handles)
% hObject    handle to pasteFromClipboardTB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%table=ocrTable();
content=clipboard('paste');
dt=double(content);
content=content(dt<=127);
len=length(content);
p=1;
p0=1;
row=1;
lines={};
while(p<=len)
    c=content(p);
    if(double(c)==10)%endline
        if(p>p0)
            lines{row}=content(p0:p-1);
            p0=p+1;
            row=row+1;
        end
    end
    p=p+1;
end

guiData=guidata(hObject);
guiData.lines=lines;
guidata(hObject,guiData);
updateGUI(hObject);

% --- Outputs from this function are returned to the command line.
function varargout = quickTabular_spmSummary_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in saveAsBT.
function saveAsBT_Callback(hObject, eventdata, handles)
% hObject    handle to saveAsBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadTableBT.
function loadTableBT_Callback(hObject, eventdata, handles)
% hObject    handle to loadTableBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(IsSPMStatSummaryTable(hObject))
    displaySPMStatSummary(hObject);
end

guiData=guidata(hObject);
tabledir=fullfile(CommonEnvMethods.getMatlabProjectPath(), 'Literature Tables');
if(~exist(tabledir,'file'))
    tabledir=pwd;
end
[name, path]=uigetfile(fullfile(tabledir,'*.txt;*.xls;*.xlsx;*.csv;*.tsv'),'select file');
fname=fullfile(path,name);
[~, ~, ext]=fileparts(fname);
tab=CommonMethods.tab;
switch ext
    case '.xlsx'
        [num,txt,raw]=xlsread(fname);
        guiData=guidata(hObject);
        guiData.table.data=raw;
        guiData.table.header=arrayfun(@(x) CommonMethods.getDefaultColumnNames(x), 1:size(raw,2), 'UniformOutput', false);
        guidata(hObject,guiData);        
    case '.csv'
    case '.tsv'
        fid=fopen(fname);
        t=textscan(fid,'%s','Delimiter','\n');
        t=t{1};
        rows=length(t);
        tbl=cellfun(@(x)strsplit(x,tab),t,'UniformOutput',false);
        lens=cellfun(@(x)length(x),tbl);
        cols=max(lens);
        data=cell(rows,cols);
        
        for r=1:rows
            tmp=tbl{r};
            len=length(tmp);
            for c=1:len
                data{r,c}=tmp{c};
            end
            for c=len+1:cols
                data{r,c}=' ';
            end
        end
%        obj=spmStatSummaryHandler();
        guiData.table.data=data;
        guiData.table.header=arrayfun(@(x) CommonMethods.getDefaultColumnNames(x), 1:cols, 'UniformOutput', false);
        guidata(hObject,guiData);    
end
updateGUI(hObject);



    function displaySPMStatSummary(hObject)
        guiData=guidata(hObject);
        
        if(get(guiData.MakingSelectionRB,'Value'))
            %Selection mode
            tbl=guiData.table;
            rows=length(guiData.rowSelection);
            cols=length(guiData.colSelection);
            tbl.rowNames=guiData.rowNames;
            tbl.colNames=guiData.colNames;
            tbl.data=cell(rows,cols);
            if(~isfield(guiData, 'rowSelection'))
                guiData.rowSelection=true(1,rows);
            end
            if(~isfield(guiData, 'colSelection'))
                guiData.colSelection=true(1,cols);
            end
            for r=1:rows
                for c=1:cols
                    if(guiData.rowSelection(r)&&guiData.colSelection(c))
                        tbl.data{r,c}='O';
                    else
                        tbl.data{r,c}='X';
                    end
                end
            end
        else
            %Table displaying mode
            p_cluster=str2num(get(guiData.pClusterTF,'String'));
            p_FWE=str2num(get(guiData.pFWETF,'String'));
            p_FDR=str2num(get(guiData.pFDRTF,'String'));
                        
            summaryFile=getSummaryFile();
            figName=['quickTabular_spmSummary: ' summaryFile];
            set(guiData.figure1,'Name', figName);
            if(isempty(summaryFile))
                setSummaryFile();
                summaryFile=getSummaryFile();
            end
            if(~isfield(guiData,'rebuildStatSummary'))
                guiData.rebuildStatSummary=true;
            end
            
            multipleComparisonCorrOption=getMultipleComparisonCorrOption(hObject);
            if(rebuildSummaryOutput(hObject, summaryFile))
                so=spmStatSummaryHandler(summaryFile,p_cluster,p_FWE,p_FDR,multipleComparisonCorrOption);
                guiData.so=so;
            else
                so=guiData.SSHandler;
            end
            %             if(guiData.rebuildStatSummary)
            %                 so=spmStatSummaryHandler(summaryFile,p_cluster,p_FWE,p_FDR,multipleComparisonCorOption);
            %             else
            %                 so=guiData.SSHandler;
            %             end
            
            tbl=so.resultsOutput.table_overview;
            guiData.rowNames=tbl.rowNames;
            guiData.colNames=tbl.colNames;
            rows=size(tbl.data,1);
            cols=size(tbl.data,2);
            
            newSelection=false;
            if(~isfield(guiData,'rowSelection'))
                newSelection=true;
            elseif(length(guiData.colSelection)~=cols||length(guiData.rowSelection)~=rows)
                newSelection=true;
            end
            
            if(newSelection)
                rowSelection=true(1,rows);
                colSelection=true(1,cols);
                guiData.rowSelection=rowSelection;
                guiData.colSelection=colSelection;
            else
                rowSelection=guiData.rowSelection;
                colSelection=guiData.colSelection;
            end
            
            if(get(guiData.SignificantRowsOnlyRB,'value'))
                sig=tbl.significant;
                rows=size(sig,1);
                selection=arrayfun(@(r)any(sig(r,colSelection)),1:rows);
                rowSelection=(rowSelection.*selection)>0;
            end
            
            if(get(guiData.SignificantColumnsOnlyRB,'value'))
                sig=tbl.significant;
                cols=size(sig,2);
                selection=arrayfun(@(c)any(sig(rowSelection,c)),1:cols);
                colSelection=(colSelection.*selection)>0;
            end
            
            tbl.data=tbl.data(rowSelection,colSelection);
            tbl.summaryNodeLookupTable=tbl.summaryNodeLookupTable(rowSelection,colSelection);
            tbl.colNames=tbl.colNames(colSelection);
            tbl.rowNames=tbl.rowNames(rowSelection);
            tbl.significant=tbl.significant(rowSelection,colSelection);
            guiData.SSHandler=so;%spmStatHandler
        end
        guiData.table=tbl;
        guidata(hObject,guiData);
        updateGUI(hObject);
        
        
        function multipleComparisonCorrOption=getMultipleComparisonCorrOption(hObject)
            guiData=guidata(hObject);
            mcro=get(guiData.MultCompCorrectionCB,'String');
            mcr=mcro(get(guiData.MultCompCorrectionCB,'Value')); %multiple comparison correction option
            guiData.MultCompCorrectionOption=mcr;
            if(strcmp(mcr,'None'))
                multipleComparisonCorrOption={};
            else
                multipleComparisonCorrOption.correctionOption=mcr;
                if(isfield(guiData,'rowSelection'))
                    multipleComparisonCorrOption.rowSelection=guiData.rowSelection;
                    multipleComparisonCorrOption.colSelection=guiData.colSelection;
                end
            end
            
        function IS=rebuildSummaryOutput(hObject, summaryFile)
            multipleComparisonCorrOption=getMultipleComparisonCorrOption(hObject);
            guiData=guidata(hObject);
            IS=true;
            if(~isfield(guiData,'so'))
                return;%first time
            end
            if(~strcmp(summaryFile, guiData.so.summaryFile))
                return;
            end
            if(~equivalentMultCompCorr(multipleComparisonCorrOption, guiData.so.multipleComparisonCorrOption))
                return;
            end
            p_cluster=str2num(get(guiData.pClusterTF,'String'));
            if(p_cluster~=guiData.so.p_cluster)
                return;
            end
            p_FWE=str2num(get(guiData.pFWETF,'String'));
            if(p_FWE~=guiData.so.p_v_FWE)
                return;
            end
            p_FDR=str2num(get(guiData.pFDRTF,'String'));
            if(p_FDR~=guiData.so.p_v_FDR)
                return;
            end
            IS=false
            
            function IS=equivalentMultCompCorr(option1,option2)
                IS=true;
                if(isempty(option1))
                    if(isempty(option2))
                        return;
                    else
                        IS=false;
                        return;
                    end
                elseif(isempty(option2))
                    if(isempty(option1))
                        return;
                    else
                        IS=false;
                        return;
                    end
                end
                if(~strcmp(option1.correctionOption,option2.correctionOption))
                    IS=false;
                    return;
                end
                if(CommonMethods.equalMatrices(option1.rowSelection, option2.rowSelection)&&CommonMethods.equalMatrices(option1.rowSelection, option2.colSelection))
                    return;
                else
                    IS=false;
                    return;
                end
                
            
function updated=updatefname(hObject)
guiData=guiDate(hObject);
fname=guiData.fname;
if(~isempty(fname))
    [initdir,~,~]=fileparts(fname);
else
    initdir=CommonEnvMethods.getMatlabProjectPath();
end

if(~isempty(initdir))
    cd(initdir);
end

dialogtitle='Select the file for the table';
%folder=uigetdir(initdir,dialogtitle);
%[name,path,filter]=uigetfile('*.csv',dialogtitle);
[name,path,filter]=uiputfile('*.csv',dialogtitle);
if(length(path)==1)
    if(path==0)
        specified=true;
    else
        specified=false;
    end
else
    specified=true;
end
if(specified)    
    guiData.fname=fullfile(path,name);
    updated=true;   
    set(guiData.fnameLB,'String',guiData.fname);
else
    updated=false;
end


    function updateGUI(hObject)
 %       hObject=getMainHObject(hObject);
 %  there seemstwo copies of handles for hObject. Turn this off to find it
 %  out on 3/1/2018
        updateTable(hObject);
        
    function hObject=getMainHObject(hObject)
        while(~MainHObject(hObject))
            hObject=get(hObject,'parent');
            if(isempty(hObject))
                break;
            end
        end
            
    function IS=MainHObject(hObject)
        IS=false;
        tmp=get(hObject);
        if(isfield(tmp,'Name'))
            IS=~isempty(strfind(tmp.Name, 'quickTabular_spmSummary'));
        end
        
        function updateTable(hObject);
            guiData=guidata(hObject);
            if(~isfield(guiData,'table'))
                return;
            end
            colNames=guiData.table.colNames;
            rowNames=guiData.table.rowNames;
            data=guiData.table.data;
            rows=size(data,1);
            cols=size(data,2);
            lens=zeros(rows,cols);%%needed for auto adjusting column widths
            fs=10;
            rows=size(data,1)+1;
            for r=1:rows
                for c=1:cols
                    if(r==1)
                        lens(r,c)=length(colNames{c});
                    else
                        lens(r,c)=length(data{r-1,c});
                    end
                end
            end
            lens=lens*fs;
            cellMaxLen=num2cell(max(lens));
            parent=get(guiData.tablePanel,'parent');
            position=getpixelposition(guiData.tablePanel,true);
            if(isfield(guiData, 'hTable'))
                delete(guiData.hTable);
            end
            hTable=uitable(parent,'Data',data,'Position', position,'Units','normalized');
            guiData.hTable=hTable;
            %            set(hTable,'ButtonDownFcn',@clickcallback);
            set(hTable,'KeyPressFcn',@keypresscallback);
            set(hTable,'KeyReleaseFcn',@keyreleasecallback);
            %            set(hTable,'Enable','Inactive');%Set 'Enable' off make it firing on left click on the table. Otherwise only firing on wright click.
            %            set(hTable,'Enable','On');%Set 'Enable' off make it firing on left click on the table. Otherwise only firing on wright click.
            jscrollpane = findjobj(hTable);
            jviewport = jscrollpane.getViewport;
            jtable = jviewport.getView;
            set(jtable, 'MouseClickedCallback', @clickcallback);
            set(hTable,'ColumnWidth',cellMaxLen,'columnName',colNames,'rowName', rowNames);
            
            if(isempty(get(hTable,'CellSelectionCallback')))
                set(hTable,'CellSelectionCallback',@uitable_CellSelectionCallback);
            end
            
            set(guiData.ClusterPstr,'String',num2str(guiData.table.p_cluster_m));
            set(guiData.VFWEPstr,'String',num2str(guiData.table.p_v_FWE_m));
            set(guiData.VFDRPstr,'String',num2str(guiData.table.p_v_FDR_m));
            guiData.table.hTable=hTable;
            guidata(hObject,guiData);
            function sns=getSelectedSummaryNodes(hObject)
                sn={};
                guiData=guidata(hObject);
                cellIndices=guiData.selectedCellIndices;
                if(~isfield(guiData,'so'))
                    return;
                end
                so=guiData.so;
                sns=arrayfun(@(x)so.statSummaryNodes{guiData.table.summaryNodeLookupTable(cellIndices(x,1),cellIndices(x,2))},1:size(cellIndices,1),'UniformOutput',false);
                
function uitable_CellSelectionCallback(hObject, eventdata, handles)
    global clickStatus;
    
    guiData=guidata(hObject);
    if(exist('eventdata','var'))
        cellIndices=eventdata.Indices;
        guiData.selectedCellIndices=cellIndices;
        guidata(hObject,guiData);
    else
        cellIndices=guiData.selectedCellIndices;
    end
    
    if(get(guiData.MakingSelectionCB,'Value')==1&&get(guiData.MakingSelectionRB,'Value'))
        
        if(~get(guiData.HoldRB,'Value'))
            inds=cellIndices(:,1);
            guiData.rowSelection(inds)=~guiData.rowSelection(inds);
            guidata(hObject,guiData);
            displaySPMStatSummary(hObject);
        end
        return;
    elseif(get(guiData.MakingSelectionCB,'Value')==2&&get(guiData.MakingSelectionRB,'Value'))
        
        if(~get(guiData.HoldRB,'Value'))
            inds=cellIndices(:,2);
            guiData.colSelection(inds)=~guiData.colSelection(inds);
            guidata(hObject,guiData);
            guidata(hObject,guiData);
            displaySPMStatSummary(hObject);
        end
        return;
    end
    
    %    dispStatImg=clickStatus==2;
    dispStatImg=false;
    if(dispStatImg)
        exceptions={guiData.mainGui};
        if(isfield(guiData,'xjView'))
            exceptions{end+1}=guiData.xjView;
        end
        CommonMethods.closeAll(exceptions);
    end
    argout=guiData.SSHandler.resultsOutput.handleTableSelection(guiData.table, cellIndices,dispStatImg);
    if(isempty(argout))
        return;
    end
    set(guiData.CellMsgStr,'String',argout.msg);
    if(isfield(argout,'xjView'))
        guiData.xjView=argout.xjView;
    end    
    guidata(hObject,guiData);
    
    updateDetailTable(hObject,argout);
    
    function updateDetailTable(hObject, pars)
        guiData=guidata(hObject);
        if(isfield(guiData,'hTable_detailed'))
            if(~isvalid(guiData.hTable_detailed))
                guiData=rmfield(guiData,'hTable_detailed');
            end
        end
        if(~isfield(guiData,'hTable_detailed'))
            quickTabular_spmStatDetails();
            hs=allchild(0);
            for i=1:length(hs)
                h=hs(i);
                try
                    name=get(h,'Name');
                catch err
                    continue;
                end
                if(strcmp(name,'quickTabular_spmStatDetails'))
                    guiData.hTable_detailed=h;
                    guidata(hObject,guiData);
                    break;
                end
            end
        end
        gd=guidata(guiData.hTable_detailed);
        gd.table=pars.table;
        gd.SSHandler=guiData.SSHandler;
        gd.summaryGui=guiData.mainGui;
        guidata(guiData.hTable_detailed,gd);
        gd.updateGUI(guiData.hTable_detailed);
        figure(guiData.hTable_detailed);
        
function updateSPMExt(hObject)
guiData=guidata(hObject);
indeces=guiData.table.indeces;
str=guiData.table.data{indeces(1), indeces(2)};
if(isfield(guiData, 'spmExtension'))
%    if(strcmp(strtrim(guiData.table.header{indeces(2)}), 'MNI coordinates (x y z)'))
    if(indeces(2)==guiData.CoordinateColumn)
        coor=CommonMethods.str2nums(str, ' ')';
        spmExtGuiData=guidata(guiData.spmExtension);
        spmExtGuiData.updateCoordinates_mni(guiData.spmExtension,coor(1:3));
    end
end

% --- Executes on button press in makeHeaderBT.
function makeHeaderBT_Callback(hObject, eventdata, handles)
% hObject    handle to makeHeaderBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function CoordinateColumnTF_Callback(hObject, eventdata, handles)
% hObject    handle to CoordinateColumnTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of CoordinateColumnTF as text
%        str2double(get(hObject,'String')) returns contents of CoordinateColumnTF as a double
guiData=guidata(hObject);
guiData.CoordinateColumn=str2double(get(guiData.CoordinateColumnTF,'String'));
guidata(hObject, guiData);

% --- Executes during object creation, after setting all properties.
function CoordinateColumnTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoordinateColumnTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in SummaryTypeCB.
    
function clickcallback(obj,evt)
    global clickStatus;
persistent chk
if isempty(chk)
      chk = 1;
      pause(0.5); %Add a delay to distinguish single click from a double click
      if chk == 1
          clickStatus=1;
          chk = [];
      end
else
      chk = [];
      clickStatus=2;
end

function keypresscallback(hObject, evnt)
%    fprintf('pressing: %s\n',evnt.Key);
    guiData=guidata(hObject);
    guiData.keyStr=evnt.Key;
    guidata(hObject, guiData);
    
function keyreleasecallback(hObject,evnt)
    guiData=guidata(hObject);
    guiData.keyStr={};
    guidata(hObject, guiData);
%    fprintf('releasing: %s\n',evnt.Key);
    



function pClusterTF_Callback(hObject, eventdata, handles)
% hObject    handle to pClusterTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pClusterTF as text
%        str2double(get(hObject,'String')) returns contents of pClusterTF as a double
displaySPMStatSummary(hObject);


% --- Executes during object creation, after setting all properties.
function pClusterTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pClusterTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pFWETF_Callback(hObject, eventdata, handles)
% hObject    handle to pFWETF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pFWETF as text
%        str2double(get(hObject,'String')) returns contents of pFWETF as a double
displaySPMStatSummary(hObject);

% --- Executes during object creation, after setting all properties.
function pFWETF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pFWETF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ExportTableBT.
function ExportTableBT_Callback(hObject, eventdata, handles)
% hObject    handle to ExportTableBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
if(~isfield(guiData,'SSHandler'))
    return
end
so=guiData.SSHandler;
name=['Stats summary table _ ' so.resultsOutput.logName '.tsv'];
[fname, path]=uiputfile('*.tsv','Specifya tab separated file to export the', fullfile(so.resultsOutput.logFolder,name));
fname=fullfile(path,fname);
fid=fopen(fname,'wt');
data=guiData.table.data;
len=size(data,1);
tab=CommonMethods.tab;
line=CommonMethods.cell2str(guiData.table.header,tab);
fprintf(fid,'%s\n',line);
for l=1:len
    line=CommonMethods.cell2str(data(l,:),tab);
    fprintf(fid,'%s\n',line);
end
fclose(fid);


function pFDRTF_Callback(hObject, eventdata, handles)
% hObject    handle to pFDRTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pFDRTF as text
%        str2double(get(hObject,'String')) returns contents of pFDRTF as a double

displaySPMStatSummary(hObject);

% --- Executes during object creation, after setting all properties.
function pFDRTF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pFDRTF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CopyTableTextBT.
function CopyTableTextBT_Callback(hObject, eventdata, handles)
% hObject    handle to CopyTableTextBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guiData=guidata(hObject);
header=guiData.table.colNames;
data=guiData.table.data;
if(size(header,2)==size(data,2))
    vertcat(header,data);
end
rowNames=vertcat({'contrasts'},guiData.table.rowNames);
data=horzcat(rowNames, vertcat(header,data));
tw=size(data,2);
%line= cell2tableline(strs,pos,tw,del,strSpace)
l1=CommonMethods.cell2tableline(guiData.table.summaryFile,1,tw,CommonMethods.tab,' ');
options=get(guiData.MultCompCorrectionCB, 'String');
multCompCorr=options{get(guiData.MultCompCorrectionCB, 'Value')};
strs=['p_cluster: ' num2str(guiData.table.p_cluster) '; p_v_FWE: ' num2str(guiData.table.p_v_FWE) '; p_v_FDR: ' num2str(guiData.table.p_v_FDR) '; MultCompCorr: ' multCompCorr];
l2=CommonMethods.cell2tableline(strs,1,tw,CommonMethods.tab,' ');
data=vertcat(strsplit(l2,CommonMethods.tab), data);
data=vertcat(strsplit(l1,CommonMethods.tab), data);
ClipboardHandler.copyCellarray(data);
 

% --- Executes on button press in SignificantRowsOnlyRB.
function SignificantRowsOnlyRB_Callback(hObject, eventdata, handles)
% hObject    handle to SignificantRowsOnlyRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SignificantRowsOnlyRB

displaySPMStatSummary(hObject);

% --- Executes on button press in SignificantColumnsOnlyRB.
function SignificantColumnsOnlyRB_Callback(hObject, eventdata, handles)
% hObject    handle to SignificantColumnsOnlyRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SignificantColumnsOnlyRB
displaySPMStatSummary(hObject);

% --- Executes on button press in SetSummaryFileBT.
function SetSummaryFileBT_Callback(hObject, eventdata, handles)
% hObject    handle to SetSummaryFileBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setSummaryFile();
displaySPMStatSummary(hObject);

    function setSummaryFile();
        [summaryFile, initFile]=getSummaryFile();
        if(isempty(summaryFile))
            initdir=pwd;
        else
            initdir=fileparts(summaryFile);
        end
        [name,path]=uigetfile(fullfile(initdir,'*.tsv; *.csv'),'select a spm stat summary file');
        summaryFile=fullfile(path,name);
        
        if(isempty(initFile))
            initFile=uigetfile('*.ini','getting file path for quickTabular_spmSummary.ini', fullfile(CommonEnvMethods.getMatlabProjectPath(), 'Literature Tables','quickTabular_spmSummary.ini'));
            lines={};
            
        else
            lines=CommonMethods.getTextLines(initFile);
        end
        if(ischar(lines))
            lines={lines};
        end
        idx=find(cellfun(@(x)~isempty(strfind(x,'SummaryFile')),lines));
        if(isempty(idx))
            ind=length(lines)+1;
        else
            ind=idx(1);
        end
        lines{ind}=['SummaryFile ' summaryFile];
        CommonMethods.writeTextLines(initFile,lines);
        
    function [SummaryFile, initFile]=getSummaryFile()
        SummaryFile='';
        initFile=which('quickTabular_spmSummary.ini');
        if(isempty(initFile))
            return;
        end
        fid=fopen(initFile);
        lines=textscan(fid,'%s','Delimiter','\n');
        lines=lines{1};
        if(ischar(lines))
            lines={lines};
        end
        idx=find(cellfun(@(x)~isempty(strfind(x,'SummaryFile')),lines));
        if(isempty(idx))
            return;
        end
        SummaryFile=lines{idx(1)};
        strs=strsplit(SummaryFile,' ');
        SummaryFile=strs{2};
        if(~exist(SummaryFile,'file'))
            SummaryFile='';
        end


% --- Executes on button press in SelectingRowsRB.
function SelectingRowsRB_Callback(hObject, eventdata, handles)
% hObject    handle to SelectingRowsRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SelectingRowsRB
displaySPMStatSummary(hObject);

% --- Executes on button press in SelectingColumnsRB.
function SelectingColumnsRB_Callback(hObject, eventdata, handles)
% hObject    handle to SelectingColumnsRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SelectingColumnsRB
displaySPMStatSummary(hObject);

% --- Executes on button press in HoldRB.
function HoldRB_Callback(hObject, eventdata, handles)
% hObject    handle to HoldRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HoldRB
guiData=guidata(hObject);
if(~get(guiData.HoldRB,'Value'))
    uitable_CellSelectionCallback(hObject);
end


% --- Executes on selection change in MakingSelectionCB.
function MakingSelectionCB_Callback(hObject, eventdata, handles)
% hObject    handle to MakingSelectionCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MakingSelectionCB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MakingSelectionCB
guiData=guidata(hObject);
hCB=guiData.MakingSelectionCB;
option=CommonMethods.getSelectedString(hCB);
if(strcmp(option,'Selecting Rows')||strcmp(option,'Selecting Columns'))
    set(guiData.MakingSelectionRB,'Value',1);
elseif(strcmp(option, 'Load Selection Function'))
    funcdir=fullfile(CommonEnvMethods.getMatlabProjectPath(), 'UI Functions', 'Selection Functions');
    if(~exist(funcdir,'file'))
        mkdir(funcdir);
    end
    [name, path]=uigetfile(fullfile(funcdir,'*.m'),'Select a table cell selection function');
    
    if(ischar(name))
        [~,name,ext]=fileparts(name);
        cmd=[name '(hObject);'];
        eval(cmd);
        set(guiData.MakingSelectionRB,'Value',1);
        ind=get(hCB,'Value');
        options=get(hCB,'String');
        options={options{1:ind-1} name options{ind:end}};
        set(hCB,'String',options);
    else
        %%for the convenience of making a new function
        
        cd(funcdir);
    end
else
    cmd=[option '(hObject);'];
    eval(cmd);
    set(guiData.MakingSelectionRB,'Value',1);
end
if(get(guiData.MakingSelectionRB,'Value'))
    displaySPMStatSummary(hObject);
end

% --- Executes during object creation, after setting all properties.
function MakingSelectionCB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MakingSelectionCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MakingSelectionRB.
function MakingSelectionRB_Callback(hObject, eventdata, handles)
% hObject    handle to MakingSelectionRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MakingSelectionRB
displaySPMStatSummary(hObject);


% --- Executes on button press in GoToDirBT.
function GoToDirBT_Callback(hObject, eventdata, handles)
% hObject    handle to GoToDirBT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sns=getSelectedSummaryNodes(hObject);
if(~isempty(sns))
    sn=sns{1};
    cd(fileparts(sn.fnameSPM));
end


% --- Executes on selection change in MultCompCorrectionCB.
function MultCompCorrectionCB_Callback(hObject, eventdata, handles)
% hObject    handle to MultCompCorrectionCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns MultCompCorrectionCB contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MultCompCorrectionCB
displaySPMStatSummary(hObject);


% --- Executes during object creation, after setting all properties.
function MultCompCorrectionCB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MultCompCorrectionCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
