classdef NeuritesNode
    %Tree node for neurites model
    %   Use to create phytree from Neurolucida tracing data
    
    properties
        id
        countid         %order added - for tree generation
        branchlength    %length of branch as float
        leftnode        %child NeuritesNode obj
        rightnode       %child NeuritesNode obj
        parentnode      %parent NeuritesNode obj
        nodetype        %Indicates whether Branch or END ie no childnodes
        anglearc        %angle of node coords to soma
        nodelevel       %level of branching with 1 as root
        vol             %volume of branch um3
        sa              %sa of branch um2
        points          %xy points
        tree            %tree number
    end
    
    methods
        function obj = NeuritesNode(id,branchlength,nodetype,nodelevel,countid,treenum)
            obj.id = id;
            obj.countid = countid;
            obj.branchlength = branchlength;
            obj.nodetype = nodetype;
            obj.nodelevel = nodelevel;
            obj.tree = treenum;
        end
        function obj = setMeasurements(obj,anglearc, vol, sa, blength, points)
            obj.anglearc = anglearc;
            obj.vol = vol;
            obj.sa = sa;
            obj.branchlength = blength;
            obj.points = points;
        end
            
        %Add child node if empty  - to change, use direct n.left=node
        function [parent,childnode] = addChildNode(obj,childnode)
            %detect nodelevel to add to
            parent = obj;
            
            %parent = obj.findNodeLevel(node.nodelevel-1,1)
            while(parent.nodelevel < childnode.nodelevel-1)
                if obj.isNode(parent.rightnode) && parent.rightnode.isBranchpoint()
                    parent = parent.rightnode;
                elseif obj.isNode(parent.leftnode) && parent.leftnode.isBranchpoint()
                    parent = parent.leftnode;
                else
                    break
                end
            end
            if (childnode.id ~= parent.id)
                childnode.parentnode = parent;   
            end
            % Add childnode to parent on next level up
            if ~obj.isNode(parent.rightnode)|| (parent.rightnode.hasID(childnode.id))
                parent.rightnode = childnode;
                
            elseif ~obj.isNode(parent.leftnode) || (parent.leftnode.hasID(childnode.id))
                parent.leftnode = childnode;
            else 
               disp('Warning:Child node NOT added')
                
            end
        end
        
        
        function bool = isBranchpoint(obj)
            bool = strcmp(obj.nodetype,'BP');
            
        end
        
        function bool = isNode(obj,node)
            bool = isa(node,'NeuritesNode');
        end
        
        function bool = hasID(obj,id)
            bool = (obj.isNode(obj) && obj.id == id);
        end
        
        function bool = hasParent(obj)
            bool = (obj.isNode(obj) && ~isempty(obj.parentnode));
        end
        
        function sibling = getSibling(obj)
            sibling = [];
            if (obj.hasParent())
                if (obj.parentnode.leftnode.id == obj.id)
                    sibling = obj.parentnode.rightnode;
                elseif (obj.parentnode.rightnode.id == obj.id)
                    sibling = obj.parentnode.leftnode;
                end
            end
        end
                
      
        % finds node on level and whether is branchpoint (0|1)
        function childnode = findNodeLevel(obj,level,branchpoint)
            %look along left path first
            childnode = obj;
            while(childnode.nodelevel < level)
                if obj.isNode(childnode.leftnode) && childnode.leftnode.isBranchpoint()==branchpoint
                    childnode = childnode.leftnode;
                elseif obj.isNode(childnode.rightnode) && childnode.rightnode.isBranchpoint()==branchpoint
                    childnode = childnode.rightnode;
                else
                    break
                end
            end
        end
                    
      
        
        %has both Child nodes
        function bool = hasChildNodes(obj)
            if obj.isNode(obj.leftnode) && obj.isNode(obj.rightnode)
                bool = true;
            else
                bool = false;
            end
        end
        
                
    end
    
end

