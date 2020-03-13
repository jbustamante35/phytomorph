classdef linkable < handle

    properties

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set of lists for pointers
        upList;
        downList;
        leftList;
        rightList;

    end


    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % constructor
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = linkable()
            % list connections to layers
            obj.upList = {};
            obj.downList = {};
            obj.leftList = {};
            obj.rightList = {};
        end

        % function to link the device to something above itself
        function [] = linkUp(obj,target)
            obj.upList{end+1} = target;
        end

        % functio to link the device to something below itself
        function [] = linkDown(obj,target)
            obj.downList{end+1} = target;
        end

        % clear all links
        function [] = clearAllLinks(obj)
            obj.upList = {};
            obj.downList = {};
            obj.leftList = {};
            obj.rightList = {};
        end


        % recursive crawl through nodes to build list
        function [nodeList,upid,M,curID] = getNodeList(obj,pid,nodeList,upid,M)

            if nargin == 1
                M = containers.Map('KeyType','char','ValueType','any');
              
                pid = 'null';
                nodeList = {'null'};
                upid = {{},{}};
            end


            
            curID = obj.getUidType();
            M(curID) = obj;
            upid = cat(1,upid,{{curID},{pid}});

            if ~any(contains(nodeList,curID))
                nodeList{end+1} = curID;
                for e = 1:numel(obj.downList)
                    [nodeList,upid,M] = obj.downList{e}.getNodeList(curID,nodeList,upid,M);
                end       
            end
            
            nodeList = unique(nodeList);

            
           
            
        end


    end


end