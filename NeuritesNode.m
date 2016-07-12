classdef NeuritesNode
    %Tree node for neurites model
    %   Use to create phytree from Neurolucida tracing data
    
    properties
        id
        branchlength    %length of branch as float
        leftnode        %child NeuritesNode obj
        rightnode       %child NeuritesNode obj
        nodetype        %Indicates whether Branch or END ie no childnodes
        anglesoma       %angle of node coords to soma
        nodelevel       %level of branching with 1 as root
    end
    
    methods
        function obj = NeuritesNode(id,branchlength,nodetype,nodelevel)
            obj.id = id;
            obj.branchlength = branchlength;
            obj.nodetype = nodetype;
            obj.nodelevel = nodelevel;
        end
        %Add child node if empty  - to change, use direct n.left=node
        function [parent,childnode] = addChildNode(obj,childnode)
            %detect nodelevel to add to
            parent = obj;
            
            %parent = obj.findNodeLevel(node.nodelevel-1,1)
            while(parent.nodelevel < childnode.nodelevel-1)
                if obj.isNode(parent.leftnode) && parent.leftnode.isBranchpoint()
                    parent = parent.leftnode;
                elseif obj.isNode(parent.rightnode) && parent.rightnode.isBranchpoint()
                    parent = parent.rightnode;
                else
                    break
                end
            end
            % Add childnode to parent on next level up
            if ~obj.isNode(parent.leftnode) || (parent.leftnode.hasID(childnode.id))
                parent.leftnode = childnode;
            elseif ~obj.isNode(parent.rightnode)|| (parent.rightnode.hasID(childnode.id))
                parent.rightnode = childnode;
            else 
               disp('Not added')
                
            end
           
            
            
        end
        
        
        function bool = isBranchpoint(obj)
            bool = (strcmp(obj.nodetype,'BP') > 0);
            
        end
        
        function bool = isNode(obj,node)
            bool = isa(node,'NeuritesNode')
        end
        
        function bool = hasID(obj,id)
            bool = (obj.isNode(obj) && obj.id == id)
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

