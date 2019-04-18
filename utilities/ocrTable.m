classdef ocrTable
    %OCTTABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rows;
        cols;
        imgdata;
        table;
        txt;
    end
    
    methods
        function obj=ocrTable(rows, cols, noTranspose)
            imgdata=imclipboard('paste');
            if(exist('noTranspose','var'))
                for t=1:noTranspose
                    data0=imgdata;
                    dim=size(imgdata);
                    imgdata=zeros(dim(2),dim(1),dim(3));
                    for i=1:dim(2)
                        for j=1:dim(1)
                            imgdata(i,j,:)=data0(j,i,:);
                        end
                    end
                end
            end
            obj.rows=rows;
            obj.cols=cols;
            obj.imgdata=imgdata;
            obj.txt=ocr(imgdata);
%            obj.txt=OCR(imgdata);
            obj.table=obj.makeTable();
        end
        function table=makeTable(obj)
            txt=obj.txt;
            rows=obj.rows;
            cols=obj.cols;
            table=cell(rows,cols);
            words=txt.Words;
            boxes=txt.WordBoundingBoxes;
           
            len=length(words);
            boxes=boxes+0.01*rand(len,4);%to avoid degeneracey of the coordinate causing problem on sorting.
            
            ys=arrayfun(@(x)boxes(x,2),1:len);%y coordinates of the top of the bounding boxes
            [oys, inds_ys]=sort(ys);
            yns=arrayfun(@(x)boxes(x,2),inds_ys);%y coordinates of the top of the bounding boxes, sorted
            yxs=arrayfun(@(x)(boxes(x,2)+boxes(x,4)),inds_ys);%y coordinates of the bottom of the bounding boxes, sorted
            inds_table=cell(rows,cols);
            unpicked=ones(1,len);
            for r=1:rows
                inds=find(unpicked);%The indexes of the words possibly belong to the row r, index that as sorted according to y
                if(isempty(inds))
                    break;
                end
                first=inds(1);%the index of the word with the smallest y
                set=[first];%collection of indexes of the words in this row, the sorted indexes
                yn=yns(first);
                yx=yxs(first);
                unpicked(first)=0;
                range=[yn yx];
                for i=2:length(inds)
                    ind=inds(i);
                    if(CommonMethods.IsOverlappingRanges(range, [yns(ind) yxs(ind)]))
                        set(end+1)=ind;
                        unpicked(ind)=0;
                    end
                end
                
                len1=length(set);
                inds_row=arrayfun(@(x) inds_ys(x), set);% indexes of the words, the original indexes
                boxes_row=arrayfun(@(x) boxes(x,:), inds_row, 'UniformOutput', false);
                words_row=arrayfun(@(x) words{x}, inds_row, 'UniformOutput', false);%the words in the row r, not sorted according to x position.                 
                xs=arrayfun(@(x)boxes_row{x}(1), 1:len1);%bounding boxes of the row r
                
                [oxs, inds_xs]=sort(xs);%oxs is the sorted x position of the boxes in row r
                ls=arrayfun(@(x) boxes_row{x}(3), inds_xs);% the length of the bonding boxes in assending order of x position
                words_row_sorted=arrayfun(@(x)words_row{x}, inds_xs, 'UniformOutput', false);
                
                delta_xs=arrayfun(@(x)(oxs(x+1)-(oxs(x)+ls(x))), 1:(len1-1));
                [odxs, inds_dxs]=sort(delta_xs);
                inds_end=zeros(1,len1);%finding the indexes in the sorted x positions of the words in the same rows marking the end of each column
                for c=1:(cols-1)
                    inds_end(inds_dxs(len1-c))=1;
                end
                inds_end=[find(inds_end) len1];%the indexes in the sorted x positions of the words in the same rows marking the end of each cell
                
                iF=0;
                for c=1:cols
                    iI=iF+1;
                    iF=inds_end(c);
%                     inds_sorted=inds_xs(iI:iF);
%                     inds_cell=arrayfun(@(x)inds_xs(x),inds_sorted); %the indexes of the words in the same cell of the table in the original order.
                    table{r,c}=arrayfun(@(x) words_row_sorted{x}, iI:iF, 'UniformOUtput', false);
                end
            end
        end
    end
    
end

