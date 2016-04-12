%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
% swc_converter.m
%
% Anthony Kellems, 11/7/07
%
% Reads in a SWC (Duke Southampton standard) morphology file.
%
% INPUTS:   inputfile = name of the .swc file
%           ploton    = 1 for plotting the cell in a window, 0 for not
%
% OUTPUT:  nc = structure of neuron characteristics
%
% Elements of this code were borrowed from Nan Xiao's neur_converterh.m
% file and its associated codes for reading .ASC files (now antiquated)

function nc = swc_converter(inputfile,ploton)

fin = fopen(inputfile, 'r');
newLine = fgetl(fin);

if nargin < 2, ploton = 1; end

fprintf('Reading .swc data for %s...',inputfile)
nc.swc = [];
j = 1;
info.shrinkage_correction = [1 1 1]';

while ischar(newLine)

    if numel(newLine) == 0              % check for blank line

    elseif newLine(1) == '#'            % check for comments that are relevant
        comm = strread(newLine(2:end),'%s','delimiter',' ');

        if isempty(comm)
        elseif strcmp(comm{1},'SOMA_AREA')
            info.soma_area = str2num(comm{2});
        elseif strcmp(comm{1},'ORIGINAL_SOURCE')
            info.original_source = makesentence(comm,1);
        elseif strcmp(comm{1},'CREATURE')
            info.creature = makesentence(comm,1);
        elseif strcmp(comm{1},'REGION')
            info.traced_region = makesentence(comm,1);
        elseif strcmp(comm{1},'FIELD/LAYER')
            info.traced_layer = makesentence(comm,1);
        elseif strcmp(comm{1},'TYPE')
            info.traced_type = makesentence(comm,1);
        elseif strcmp(comm{1},'CONTRIBUTOR')
            info.contributor = makesentence(comm,1);
        elseif strcmp(comm{1},'REFERENCE')
            info.traced_reference = makesentence(comm,1);
        elseif strcmp(comm{1},'RAW')
            info.raw_file = makesentence(comm,1);
        elseif strcmp(comm{1},'EXTRAS')
            info.traced_extras = makesentence(comm,1);
        elseif strcmp(comm{1},'SHRINKAGE_CORRECTION')
            info.shrinkage_correction(1,1) = str2num(comm{2});
            info.shrinkage_correction(2,1) = str2num(comm{3});
            info.shrinkage_correction(3,1) = str2num(comm{4});
        elseif strcmp(comm{1},'VERSION_NUMBER')
            info.traced_version_number = makesentence(comm,1);
        elseif strcmp(comm{1},'VERSION_DATE')
            info.traced_version_date = makesentence(comm,1);
        elseif strcmp(comm{1},'SCALE')
            info.scale(1,1) = str2num(comm{2});
            info.scale(2,1) = str2num(comm{3});
            info.scale(3,1) = str2num(comm{4});
        end
    else
        nextLine = sscanf(newLine,'%f')';
        if length(nextLine) ~= 7
            fprintf('\n    ERROR READING .SWC FILE:\n    Improper formatting of data (missing data).')
            fprintf('\n    The bad line was:    %s\n\n',newLine)
            nc = [];
            return
        end
        nc.swc = [nc.swc; nextLine];
    end
    newLine = fgetl(fin);
    
    if mod(j,100) == 0
        fprintf(1,'.');
    end
    if mod(j,4000) == 0
        fprintf('\n                              .')
    end
    j = j+1;
end
fclose(fin);

npts = size(nc.swc,1);

% apply correction factor
% nc.swc(:,3:5) = [nc.swc(:,3)*info.shrinkage_correction(1) ...
%                  nc.swc(:,4)*info.shrinkage_correction(2) ...
%                  nc.swc(:,5)*info.shrinkage_correction(3)];

somaind = find(nc.swc(:,2) == 1);     % obtain soma data
nc.somadata.x = nc.swc(somaind,3);
nc.somadata.y = nc.swc(somaind,4);
nc.somadata.z = nc.swc(somaind,5);
nc.somadata.r = nc.swc(somaind,6);
if isfield(info,'soma_area')
    nc.somadata.As = info.soma_area;
else
    % note: this method can be very bad because of potentially intersecting
    % polygonal edges. Also, one hopes that each cell really does have the
    % soma area listed as a comment in the header. But this is here just in
    % case.
    fprintf('\nNo soma area found in header\n')
    area = 0;                           % from TK's CAAM 420 Hw 4 code
    for j = 1:length(nc.somadata.x)
        if j == length(nc.somadata.x)
            k = 1;
        else
            k = j+1;
        end
        x1 = nc.somadata.x(j);
        x2 = nc.somadata.x(k);
        y1 = nc.somadata.y(j);
        y2 = nc.somadata.y(k);
        norm_n = sqrt( (y2-y1)^2 + (x1-x2)^2 );
        
        n1 = (y2-y1)/norm_n;
        n2 = (x1-x2)/norm_n;
        
        xmid = x1 + (x2-x1)/2;
        ymid = y1 + (y2-y1)/2;
        
        h = sqrt( (x2-x1)^2 + (y2-y1)^2 );
        
        area = area + h*(xmid*n1 + ymid*n2)/2;
    end
    
    nc.somadata.As = area;
end

fprintf('\nStoring coordinates')

branchpts = find(nc.swc(:,1) ~= (nc.swc(:,7)+1));             % find where new parents occur

changepts = find(nc.swc(1:end-1,2) ~= nc.swc(2:end,2)) + 1;   % find where different types change

chunkpts = union(branchpts,changepts);           % points where new chunks begin (including soma)
totalchunks = length(chunkpts);

chunkpts(end+1) = npts+1;

maxrad = 0;
numdends = 0;
numroots = 0;

for j = 1:totalchunks
    jj = chunkpts(j);
    jjend = chunkpts(j+1) - 1;
    
    if ismember(nc.swc(jj,1),somaind)     % skip over soma chunks
        continue
    else
        numdends = numdends + 1;                % store all coordinate data for this branch
        startpt = nc.swc(jj,7);
        nc.celldata{numdends}(:,1) = nc.swc([startpt jj:jjend],3);  % x
        nc.celldata{numdends}(:,2) = nc.swc([startpt jj:jjend],4);  % y
        nc.celldata{numdends}(:,3) = nc.swc([startpt jj:jjend],5);  % z coords
        nc.celldata{numdends}(:,4) = nc.swc([startpt jj:jjend],6);  % radius
        nc.celldata{numdends}(:,5) = nc.swc([startpt jj:jjend],1);  % swc segment ID
        nc.celldata{numdends}(:,6) = nc.swc([startpt jj:jjend],2);  % segment type
        
        bigrad = max(nc.celldata{numdends}(:,4));              % check for new max radius
        if bigrad > maxrad, maxrad = bigrad; end
        
        rootcheck = ismember(nc.swc(jj:jjend,7),somaind);   % check if branch is a root    
        if sum(rootcheck) > 0
            numroots = numroots + 1;
            nc.roots(numroots) = numdends;
        end
        
        nc.parent(numdends) = startpt;          % segment index of parent (not branch number!)
    end

    if mod(j,5) == 0
        fprintf(1,'.');
    end
end

% find parents of each branch
for j = 1:numdends
    if ismember(nc.parent(j),somaind)           % if soma is parent
        nc.parent(j) = 0;
    else
        for k = 1:numdends
            if ismember(nc.parent(j),nc.celldata{k}(:,5))  % if branch k is parent
                nc.parent(j) = k;
                break
            end
        end
    end
end

% Check for branches that are "hidden" in the input file. These branches
% occur because the user tracing the actual neuron went in a successive
% line past a fork point, meaning that multiple branches are interpreted as
% one long branch.
% This is corrected by checking all fork points and creating branches as
% needed.
fprintf('\nFinding all branches...')
for j = 1:numdends
    if nc.parent(j) ~= 0          % if soma is NOT parent
        par = nc.parent(j);
        for k = 2:size(nc.celldata{par},1)-1                   % only check between endpoints
            if nc.celldata{par}(k,1:3) == nc.celldata{j}(1,1:3)
                numdends = numdends + 1;
                nc.celldata{numdends} = nc.celldata{par}(k:end,:);    % add new branch
                nc.celldata{par}(k+1:end,:) = [];                  % remove old points from parent
                nc.parent(numdends) = par;
                break
            end
        end
    end
    if mod(j,ceil(numdends/10)) == 0
        fprintf(1,'.');
    end
end

% Now go back and find right parents, as indices have changed with addition
% of new branches.
fprintf('\nGetting parents...')
for j = 1:numdends
    if nc.parent(j) == 0          % if soma is parent
        nc.parent(j) = 0;
    else
        for k = 1:numdends
            if nc.celldata{k}(end,1:3) == nc.celldata{j}(1,1:3)  % if branch k is parent
                nc.parent(j) = k;
                break
            end
        end
    end
    if mod(j,ceil(numdends/10)) == 0
        fprintf(1,'.');
    end
end

% find children of each branch
for j = 1:numdends
    children = find(nc.parent == j);
    if isempty(children)
        nc.treedata{j} = [];
    else
        nc.treedata{j} = children;
    end
end

% Remove possible redundant points
for j = 1:numdends
    numslink = 0;
    k = 1;
    while k < size(nc.celldata{j},1)-numslink 
        
        if (nc.celldata{j}(k,1) == nc.celldata{j}(k+1,1) && ...
            nc.celldata{j}(k,2) == nc.celldata{j}(k+1,2) && ...
            nc.celldata{j}(k,3) == nc.celldata{j}(k+1,3)),
            
            nc.celldata{j}(k,:) = [];
            numslink = numslink + 1;
        end
        k = k + 1;
    end
end

% Check for "lone child" branches 
% (weird case, but it's in some of the traces)
numsingle = 0;

% index shows whether or not the segment is a lone child
singleindex = zeros(1,numdends);        

% combines child and parent into one segment
for j = 1:numdends                      
    
    % checks for lone children
    if length(nc.treedata{j}) == 1         

        % if so, store number of points in parent
        parentlength = size(nc.celldata{j},1);                 
        
        % store number of points in child
        childlength = size(nc.celldata{nc.treedata{j}(1)},1);     
        
        % merge parent with child
        nc.celldata{j}(parentlength+1:parentlength+childlength-1,:) ...    
                     = nc.celldata{nc.treedata{j}(1)}(2:childlength,:);
                 
        % memorize that this segment is a lone child
        singleindex(nc.treedata{j}(1)) = 1;

        % increment number of lone children
        numsingle = numsingle + 1;              
    end
end

% re-adjust non-single child segments
if numsingle > 0                      

    newlength = size(nc.celldata,2) - numsingle;
    
    % initialize new cell array
    newnc.celldata{newlength} = [];                
    
    % initialize new tree array
    newnc.treedata{newlength} = [];                
    
    % For each lone child (that will be deleted from the cell array) ...
    % every index in tree data that is greater than the index of the lone...
    % child must be decremented by one.
   
    % traverse the segments
    for j = 1:length(singleindex)         
        
        % check if the branch is a lone child
        if singleindex(j) == 1                  
            
            % traverse the tree data
            for k = 1:size(nc.treedata,2)          
                
                % check if segment has 2 children (normal case)
                if length(nc.treedata{k}) == 2     
                    
                    % if so, check if children index are larger than ... 
                    % that of the current segment
                    if nc.treedata{k}(1) > j       
                        nc.treedata{k}(1) = nc.treedata{k}(1) - 1;    
                    end
                    
                    if nc.treedata{k}(2) > j
                        nc.treedata{k}(2) = nc.treedata{k}(2) - 1;
                    end
                end
            end
            
            % traverse nc.roots
            for i = 1:length(nc.roots)             
                % decrement root indices along the same lines
                if nc.roots(i) > j                 
                    nc.roots(i) = nc.roots(i) - 1;    
                end
            end
            
            nc.parent(j) = [];
        end
    end
    
    k = 1;
    i = 1;
    
    % Remove the lone child from the original arrays by copying ...
    % the normal segments into new arrays
    
    % traverse segments
    for j = 1:numdends                      
        % copy only if segment is not a lone child
        if singleindex(j) == 0              
            newnc.celldata{k} = nc.celldata{j};
            k = k + 1;
        end

        % copy only if tree data does not point to lone child
        if length(nc.treedata{j}) ~= 1         
            newnc.treedata{i} = nc.treedata{j};
            i = i + 1;
        end
    end    
    
    % new fixed arrays  
    nc.celldata = newnc.celldata;                    
    nc.treedata = newnc.treedata;     
    
    % new fixed numdends   
    numdends = size(nc.celldata,2);            
    
    clear newnc.celldata;
    clear newnc.treedata;   
end

% Now go back and find right parents, as indices have changed again
fprintf('\nGetting parents...')
for j = 1:numdends
    if nc.parent(j) == 0          % if soma is parent
        nc.parent(j) = 0;
    else
        for k = 1:numdends
            if nc.celldata{k}(end,1:3) == nc.celldata{j}(1,1:3)  % if branch k is parent
                nc.parent(j) = k;
                break
            end
        end
    end
    if mod(j,ceil(numdends/10)) == 0
        fprintf(1,'.');
    end
end

% find children of each branch
for j = 1:numdends
    children = find(nc.parent == j);
    if isempty(children)
        nc.treedata{j} = [];
    else
        nc.treedata{j} = children;
    end
end

% fudge with trifurcating branches
% Add a minute segment (one point from a trifurcation child)
% Change connections so that new segment connects to two of the...
% parent's children; the parent still is connected to its first ...
% child but is now also connected to the new segment.
for j = 1:numdends
    
    if length(nc.treedata{j}) == 3
        newindex = numdends + 1;

        % check for s segment with only two points, and interpolate to add
        % a third point
        if size(nc.celldata{j},1) == 2
            xtmp = linspace(nc.celldata{j}(1,1),nc.celldata{j}(2,1),3)';
            ytmp = linspace(nc.celldata{j}(1,2),nc.celldata{j}(2,2),3)';
            ztmp = linspace(nc.celldata{j}(1,3),nc.celldata{j}(2,3),3)';
            rtmp = linspace(nc.celldata{j}(1,4),nc.celldata{j}(2,4),3)';
            segtmp = [nc.celldata{j}(:,5); npts+1];
            typetmp = nc.celldata{j}(1,6)*ones(3,1);
            
            npts = npts+1;
            
            nc.celldata{j} = [xtmp ytmp ztmp rtmp segtmp typetmp];
        end
        
        nc.celldata{newindex}(1:2,:) = nc.celldata{j}(end-1:end,:);       % create new segment
        nc.celldata{j}(end,:) = [];                                   % delete point from parent
        nc.celldata{nc.treedata{j}(1)}(1,:) = nc.celldata{j}(end,:);
        
        nc.treedata{newindex}(1) = nc.treedata{j}(2);
        nc.treedata{newindex}(2) = nc.treedata{j}(3);
        nc.treedata{j}(2) = newindex;
        nc.treedata{j}(3) = [];
        
        nc.parent(newindex) = j;
        nc.parent(nc.treedata{newindex}(1)) = newindex;
        nc.parent(nc.treedata{newindex}(2)) = newindex;
        
        numdends = numdends + 1;
    end
end

nc = neur_reorder(nc);
nc.info = info;
nc.type = 'swc';
nc.filename = inputfile;
nc.maxrad = maxrad;

% Write down branch lengths
fprintf('\nWriting branch lengths...')

nc.lendata{numdends} = [];
for j = 1:numdends
    
    nc.lendata{j} = ...
         getlength(nc.celldata{j}(:,1),...
                   nc.celldata{j}(:,2),...
                   nc.celldata{j}(:,3));
               
    if mod(j,ceil(numdends/10)) == 0  
        fprintf(1,'.');
    end
    
    nc.segsize(j) = size(nc.celldata{j},1);
end

rootcount = length(nc.roots);
fprintf('\n');
fprintf('\n');
disp(['Roots: ' num2str(rootcount)]);
disp(['Segments: ' num2str(numdends)]);

if ploton == 1
    swc_plotcell(nc)
end




% plot the cell
function swc_plotcell(nc)
figure( ...
        'Name','Cell Viewer', ...
        'NumberTitle','off', ...
        'Visible','on', ...
        'Color','black',...
        'BackingStore','off'); hold on
set(gca,'Color','black')
    
for ii = 1:length(nc.lendata)
    plot3(nc.celldata{ii}(:,1),...
            nc.celldata{ii}(:,2),...
            nc.celldata{ii}(:,3),...
            'b',...
            'LineWidth',ceil(2*nc.celldata{ii}(1,4)/nc.maxrad));
end

patch(nc.somadata.x,...
        nc.somadata.y,...
        nc.somadata.z,'r','EdgeColor','r')
    
% plot branch numbers
segcount = length(nc.celldata);

for i = 1:segcount
    hold on
    h = text(nc.celldata{i}(1,1),...
        nc.celldata{i}(1,2),...
        nc.celldata{i}(1,3),num2str(i));

    set(h,...
        'Color','cyan',...
        'FontSize',8);
end

xlim = get(gca,'XLim');
ylim = get(gca,'YLim');

mlen = round(1/8*(xlim(2)-xlim(1))*1e-2)*1e2;

xp = [xlim(1) xlim(1)+mlen];
yp1 = [ylim(2)-1/8*(ylim(2)-ylim(1)) ylim(2)-1/8*(ylim(2)-ylim(1))];
yp2 = [ylim(2)-3/16*(ylim(2)-ylim(1)) ylim(2)-3/16*(ylim(2)-ylim(1))];

plot(xp,yp1,'c')
mscale = text(xp(1),yp2(1),[num2str(mlen) '\mum']);
set(mscale,'Color','Cyan','Fontsize',12)



% Takes a cell array of strings as input and returns one string with spaces
% separating each string in the original array.
% Leaves out the first string if 1 is given as argument flag.
function sentence = makesentence(cellinput,flag)

if nargin < 2, flag = 0; end

if flag == 1
    start = 2;
else
    start = 1;
end

sentence = [];
for j = start:length(cellinput)
    sentence = [sentence ' ' cellinput{j}];
end




% perform a Hines reordering of the data
function newnc = neur_reorder(nc)
nc
k = 1;
parent = 0;
children = 1;

for ii = 1:length(nc.treedata)
    if size(nc.treedata{ii},2) > 0
    % find parents
        
        parent(k) = ii;
        children(k) = nc.treedata{ii}(1);
        k = k + 1;
        parent(k) = ii;
        children(k) = nc.treedata{ii}(2);
        k = k + 1;      
    end
end

pardex = zeros(1,length(nc.treedata));
pardex(children) = parent;

% loop through branches to find depth
for ii = 1:length(nc.treedata)
    
    jj = ii;
    depth(ii) = 1;
    jj = pardex(jj);
    
    while jj
        jj = pardex(jj);
        depth(ii) = depth(ii) + 1;
    end
end

% rearrange the branches by depth (highest to lowest)
m = 1;
maxdepth = max(depth);
newnc.somadata = nc.somadata;

for jj = maxdepth:-1:1

    cdex = find(depth == jj);
    
    for kk = 1:length(cdex)
    
        ii = cdex(kk);
        
        newnc.celldata{m} = nc.celldata{ii};
        newnc.treedata{m} = nc.treedata{ii};
        
        newdex(ii) = m;
        m = m + 1;
    end    
end

% take care of new connections and reorder the roots
for ii = 1:length(nc.treedata)
    
    if size(newnc.treedata{ii},2) > 0
        newnc.treedata{ii}(1) = newdex(newnc.treedata{ii}(1));
        newnc.treedata{ii}(2) = newdex(newnc.treedata{ii}(2));
    end
end

for ii = 1:length(nc.roots)
    newnc.roots(ii) = newdex(nc.roots(ii));
end

% reversal of fortune
for ii = 1:length(nc.treedata)
    newnc.celldata{ii} = flipud(newnc.celldata{ii});
end   

% store depths, leaves, and parents
m = 1;
newnc.parent = zeros(1,length(newnc.treedata));
for ii = 1:length(depth)

    newnc.depth(newdex(ii)) = depth(ii);

    if(length(newnc.treedata{ii}) == 0)     % leaf means it has no children
        newnc.leafs(m) = ii;
        m = m + 1;
    else
        newnc.parent(newnc.treedata{ii}) = ii;
    end
end



% computes the cumulative length of a given branch
function arc = getlength(x,y,z)
sum = sqrt((x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2 + (z(2:end)-z(1:end-1)).^2);
arc = [0 cumsum(sum)'];
