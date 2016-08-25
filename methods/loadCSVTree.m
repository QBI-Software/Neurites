function neuron = loadCSVTree( csvdata, soma, scale )
%Generate dendrogram style tree from csvdata
%   Uses Matlab Phylogenetic tree: phytree 
% 1. Extract branch data per tree/branch order from csv
% 2. Determine orientation of branch as key (degrees)
% 3. List branches in order of order and location
CSV = readtable(csvdata);
t = unique(CSV.Tree);
u = unique(CSV.Order); %COUNT NUMBER ep AND bp
neuron = {};
for j=1:length(t)
    csvidx = find(CSV.Tree ==j); %indices of matching rows
    branches = {}; %cell(t,u);
    neuron{j} = branches;
    blength = 0;
    bcount = 1;
    bvol = 0;
    bsa = 0;
    points =[];
    atypes=[];
    n0 = []; %parentnode
    idx = 1; %initial 
    for i=1:length(csvidx)
        order = CSV.Order(i);
        btype = char(CSV.PointType(i));
        blength = blength + CSV.Length__m_(i);
        bvol = bvol + CSV.Volume__m__(i);
        bsa = bsa + CSV.SurfaceArea__m__(i);
        if (strcmp(btype, 'CP')==0)
            bradians = findAngleSoma(CSV.StartX(i),CSV.StartY(i),soma,scale);
            bangle = radtodeg(abs(bradians));
            n = NeuritesNode(bcount,blength,btype,order);
            n = n.setMeasurements(bangle,bvol,bsa,blength,points)
            if (~isempty(n0) && isa(n0,'NeuritesNode'))
                    %%TODO: set root node to parent of this node
                    if (n0.nodelevel ~= n.nodelevel -1)
                        [n0,idx] = findParent(branches, n.nodelevel-1);
                    end
                    [n0,n] = n0.addChildNode(n);
                    branches{n0.nodelevel}(idx,1) = n0; %replace
                    if length(branches) >= order
                        branches{order} = cat(1,branches{order},n);
                    else
                        branches{order} = [n];
                    end

            else
                n0 = n; %set root
                branches{n0.nodelevel}(1,1) = n0; %add
            end

            blength = 0;
            bcount = bcount+1;

        end

    end
    
end


end

%move this function to class?
function [N,m] = findParent(branches, level)
   N = 0;
   %Find end branches - first string
   for m=1:length(branches{1,level})
        N = branches{1,level}(m,1);
        if (N.isBranchpoint() && ~N.hasChildNodes())
            disp('Parent node ')
            N
            break;
        end
   end
end
   
            
            
