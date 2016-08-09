function tree = generateTree( csvdata )
%Generate dendrogram style tree from csvdata
%   Uses Matlab Phylogenetic tree: phytree 
% 1. Extract branch data per tree/branch order from csv
% 2. Determine orientation of branch as key (degrees)
% 3. List branches in order of order and location
CSV = readtable(csvdata);
t = unique(CSV.Tree);
u = unique(CSV.Order); %COUNT NUMBER ep AND bp
csvidx = find(CSV.Tree ==1); %indices of matching rows
branches = {}; %cell(t,u);
blength = 0;
bcount = 1;
atypes=[];
n0 = []; %parentnode
for i=1:length(csvidx)
    order = CSV.Order(i);
    btype = char(CSV.PointType(i));
    blength = blength + CSV.Length__m_(i);
    if (strcmp(btype, 'CP')==0)
        n = NeuritesNode(bcount,blength,btype,order)
        if (~isempty(n0) && isa(n0,'NeuritesNode'))
                %%TODO: set root node to parent of this node
                if (n0.nodelevel ~= n.nodelevel -1)
                    n0 = findParent(branches, n.nodelevel-1);
                end
                [n0,n] = n0.addChildNode(n);
                branches{n0.nodelevel}(1,1) = n0; %replace
                if length(branches) >= order
                    branches{order} = cat(1,branches{order},n);
                else
                    branches{order} = [n];
                end
           
        else
            n0 = n; %set root
            branches{n0.nodelevel}(1,1) = n0; %add
        end
        
%         if (length(branches) >= order)
%             branches{n0.nodelevel}(1,1) = n0; %replace
%             branches{order} = cat(1,branches{order},n);
%         else
%             branches{n0.nodelevel}(1,1) = n0;
%         end
        blength = 0;
        bcount = bcount+1;
        
    end
        
end
%load tree vars
% branchdata: contains branch pairs by IDs eg[1 2;3 4];
% branchdistances: 
branchdata = [];
branchdistances = [];% [290.7;125;324.2;139.7;0]
for j=[1:2:bcount-1]
    branchdata = cat(1, branchdata,[j j+1]);    
end
num=0;
tformat = '';
for j=[length(branches):-1:1]
    %Find end branches - first string
    for (m=1:length(branches{1,j}):1)
        N = branches{1,j}(m,1);
        if (N.isBranchpoint())
            branchdistances= cat(1, branchdistances, N.branchlength);
        end
    end
    
    
    if(~isempty(tformat) && ~isempty(bps))
        for k=1:length(bps)
            tformat = sprintf('%s:%0.02f,', tformat,bps(k))
            epformat=[];
            for i=1:length(eps)
                num = num + 1;
                epformat{i}=sprintf('%d:%0.02f',num, eps(i));
            end
            tformat = sprintf('(%s%s)',tformat,strjoin(epformat,','))
        end
    else
        epformat=[];
        
        for i=1:length(eps)
            num = num + 1;
            epformat{i}=sprintf('%d:%0.02f',num, eps(i));%TODO Detect pairs here
        end
        
    end
    
    
    
end
%Check lengths
sprintf('Num branches: %d',length(branches)) 
sprintf('Num nodes: %d',bcount)

tree = phytree(branchdata,branchdistances)
names = get(tree,'LeafNames')

%Show plot
view(tree);
%h = plot(tree,'orient','top');
xlabel('Distance from soma (um)')
%set(h.terminalNodeLabels,'Rotation',65)

end

%move this function to class?
function N = findParent(branches, level)
   N = 0;
   %Find end branches - first string
   for m=1:length(branches{1,level}):1
        N = branches{1,level}(m,1);
        if (N.isBranchpoint() && ~N.hasChildNodes())
            disp('Parent node ')
            N
            break;
        end
   end
end
   
            
            
