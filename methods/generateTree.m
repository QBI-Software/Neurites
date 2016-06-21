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
bcount = 0;
atypes=[];
for i=1:length(csvidx)
    order = CSV.Order(i);
    btype = CSV.PointType(i);
    blength = blength + CSV.Length__m_(i);
    if (strcmp(btype, 'CP')==0)
        if (strcmp(btype, 'BP')>0)
            stype = 1; % BP
        else
            stype = 2; % EP
        end
        if (length(branches) >= order)
            branches{order}{stype} = cat(1,branches{order}{stype},blength);
        else
            branches{order}{stype} = [blength];
        end
        blength = 0;
        bcount = bcount+1;
        
    end
        
end
%load tree vars
branchdata = [];%[1 2;3 4];
branchdistances = [];% [290.7;125;324.2;139.7;0]
for j=[1:2:bcount-1]
    branchdata = cat(1, branchdata,[j j+1]);    
end
num=0;
tformat = '';
for j=[length(branches):-1:1]
    eps = branches{1,j}{2};
    bps = branches{1,j}{1};
    branchdistances= cat(1, branchdistances, eps);
    
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

