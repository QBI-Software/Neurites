function neuron = loadCSVTree( csvdata, soma, scale )
%Generate dendrogram style tree from csvdata
%   Uses Matlab Phylogenetic tree: phytree 
% 1. Extract branch data per tree/branch order from csv
% 2. Determine orientation of branch as key (degrees)
% 3. List branches in order of order and location
if (isstr(csvdata))
    CSV = readtable(csvdata);
else
    CSV = csvdata
end
t = unique(CSV.Tree);
%u = unique(CSV.Order); %COUNT NUMBER ep AND bp
neuron = cell(length(t),1);
for j=1:length(t)
    csvidx = find(CSV.Tree ==j); %indices of matching rows
    branches = {};
    %neuron{j} = {};
    %branches = neuron{j};
    blength = 0;
    bcount = 1; 
    bvol = 0;
    bsa = 0;
    points =[];
    xpoints=[];
    ypoints=[];
    atypes=[];
    n0 = []; %parentnode
    idx = 1; %initial 
    for k=1:length(csvidx)
        i = csvidx(k); %get row idx
        order = CSV.Order(i);
        btype = char(CSV.PointType(i));
        bid = i;
        blength = blength + CSV.Length__m_(i);
        bvol = bvol + CSV.Volume__m__(i);
        bsa = bsa + CSV.SurfaceArea__m__(i);
        xpoints(end+1)=CSV.StartX(i);
        ypoints(end+1)=CSV.StartY(i);
        if (strcmp(btype, 'CP')==0)
            bradians = findAngleSoma(CSV.StartX(i),CSV.StartY(i),soma,scale);
            bangle = radtodeg(abs(bradians));
            xpoints(end+1)=CSV.EndX(i);
            ypoints(end+1)=CSV.EndY(i);
            points =cat(2,xpoints',ypoints');
            n = NeuritesNode(bid,blength,btype,order,bcount,CSV.Tree(i));
            n = n.setMeasurements(bangle,bvol,bsa,blength,points);
            points =[];
            xpoints=[];
            ypoints=[];
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
            bvol = 0;
            bsa = 0;
            blength = 0;
            bcount = bcount+1;

        end

    end
    neuron{j} = branches;
end


end

%move this function to class?
function [N,m] = findParent(branches, level)
   N = 0;
   sprintf('Level %d', level)
   %Find end branches - first string
   for m=1:length(branches{1,level})
        N = branches{1,level}(m,1);
        if (N.isBranchpoint() && ~N.hasChildNodes())
            sprintf('Parent node %d',N.id);
            %N
            break;
        else
            N = 0;
        end
   end
end
   
            
            
