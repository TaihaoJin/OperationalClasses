classdef xjviewActionHandler
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        xjviewActionListeners;
    end
    
    methods
        function obj = xjviewActionHander(~)
            xjviewActionListeners = {obj};
        end
        function addActionListener(obj, listener)
            xjviewActionListeners{end+1}=listener;
        end
        function actionFired(obj,action)
            for i=1:length(obj.xjviewActionListeners)
                listener=obj.xjviewActionListeners{i};
                listener.handleAction(action);
            end
        end
        function handleAction(action)
            name=action;
        end
        function updateCoordinates(obj,mmi,cor)
            mmi0=mmi;
            cor0=cor;
        end
    end
    
end

