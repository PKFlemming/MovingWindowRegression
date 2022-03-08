%%%%%%
% Purpose: assess spatial variation in regression coefficients.
% Mechanism:
% MWR (Moving Window Regression) is like a simpler version of GWR
% (Geographically Weighted Regression). It takes map data for a dependent
% variable Y and some number of independent variables X1, X2...Xn. Data
% format is described below. The arrays in have the shape of the map data,
% i.e. each array could be mapped directly onto the map using its i,j
% coordinates in the array as lon and lat (with appropriate scaling and 
% translation constants).
% User specifies a bandwidth bw and a spacing sp. The code then scrolls
% through the Y and X1...Xn arrays, starting in the top left corner and taking 
% for each variable a square subarray corresponding to a map "window" of 
% edge length bw. It runs a linear regression on these subarrays and writes 
% the regression data (r2, pvalue, coefficients etc.) out to csv. It then 
% moves to the right by a distance equal to the specified spacing and takes 
% a new set of subarrays (the next window). It continues along the row 
% until it reaches the end, then moves down by sp and starts again at the 
% beginning of the next row.
%%%%%%

%%%% USER INPUT %%%%
% enter file path info for input/output data
path = 'D:\UKuni\3rdYr\Project\GIS_data\DissPy\Regression\resampled_asciis\snap_15_500_full\';

% filenamearray is array of input file names
% Inputs are .txt files written out from ArcMap with Raster to ASCII tool
% First one needs to be the y vals (i.e., SEG or log(SEG) vals)
% Rest are x vals (independent variables)
% Data is set up such that the filenames have the form variablename +
% batchtag + .txt, so the filename without batchtag can be used for
% labelling columns in the output data
% The filename array contains all the possible variables. For each run, 
% variables not to be included are commented out. So this run modelled
% log(SEG density) in terms of GTI, AODb, slope, powerline density and
% (log(population density))_(averaged at 30km)
filenamearray = {
    'pv'  % original var name "PV", not "SEG". Check back-compatibility
%     'L10pv' 
    'GTI'
    'aodb'          
    'T'
    'slope'          
    'load'
    'pl_dens'
    'L10PCI'
%     'L10TDP'
%     'ed_rds_dv'
    'logpop_fs20'
%     'logpop_fs30'
%     'dummy_in_R1'    % used to check windowing against python translation
%     'dummy_in_R2'
%     'dummy_in_R3'
    };

windowtest = 0 % used to toggle values for python translation testing

batchtag = '_15_500_full.txt'; % invariant for each dataset
if windowtest == 1
    batchtag = '.txt'; % used to check windowing against python translation
end
ncheckval = 2; % which of the vars you want exported as a checkval
               % to doublecheck the imported points are in the right place
               % on the map. Y is 1, X1 is 2, Xn is n+1
% path to write the data out
outpath = 'D:\UKuni\3rdYr\Project\GIS_data\DissPy\Regression\resampled_asciis\snap_15_500_full\GDR_logpv';
% first bit of the outfile name, specified by user at each run to identify 
% output files:
outfilestem = 'linPV_GATSLPldlPcilPFS20';
% outfiletag could just be .csv, or can include an additional version 
% identifier ('b.csv'; 'version_2.csv') if running slightly modified
% versions of the same basic regression
outfiletag = '_ml_LonLat1.csv';      
% dimensions, in meters, of map cells:
[cellx, celly] = deal(500,500);

if windowtest == 1
    [cellx, celly] = deal(1000,1000);
end

% coords of top left corner of map (get from Arc):
[refx, refy] = deal(398334.8956, 2032383.292); 

% we're going to create subarrays by slicing into our main arrays
% first we specify the bandwidth (edge length) of our subarrays, in km:
bw = 200;
% next we specify the spacing, in km, that we want between our focal points 
% (the subarray centres). If sp = bw, there will be no overlap, and each 
% cell will be counted only once. sp << bw -> lots of overlap -> focus will 
% move more slowly over the map, producing a more smoothly varying output.
sp = 50;

if windowtest == 1
    bw = 4;
    sp = 4;
end

% Now we set the maximum fraction (0-1) of a subarray that can be NaN 
% before we reject the subarray. Tried playing with this but for non-log 
% SEG 0.5 seems best. Otherwise you get centres plotted that are off the 
% edge of the mapped data.
% With logPV data it needs to be higher though- there's gaps in the data bc
% all 0 vals -> NoData. NaNmax = 0.7 works at bw = 150. Presumably the
% further one zooms in the lower NaNmax can be. 0.6 works for bw = 150.
NaNmax = 0.7;

%%%% END OF USER INPUT %%%%

t_start = now; % timer. Used in dev to compare speed of diff approaches.

%%%% READ BLOCK %%%%
dimsfilenamearray = size(filenamearray); % dimensions of array of filenames
nvars = dimsfilenamearray(1)   % n variables (incl. y) = n files in array
YXarr = {};                     % empty cell array to populate w data
for n = 1 : nvars
    filename = filenamearray{n};                   % get nth filename
    filename = strcat(filename, batchtag);         % append file suffix
    vals = dlmread(fullfile(path, filename));      % read data
    dimsvals = size(vals);
    vdimvals = dimsvals(1); % vertical dimension of input array
    hdimvals = dimsvals(2); % horizontal dimension of input array
    formatSpec = 'vertical x horizontal dimensions of input data for nvar = %d: %dx%d';
%     str = sprintf(formatSpec,n,vdimvals,hdimvals)
    YXarr{n} = vals;                               % add data to array
end
disp("size of YXarr:")
disp(size(YXarr))
% YXarr now holds all the data. YXarr has n cells at the top level, where n
% = nvars = number of variables (Y and X). Each one of those cells
% contains an i by j array, where i and j are the dimensions of the map (in
% gridcells, not metres)

%%%% RUN BLOCK %%%%
% First we make the output filepath using the kilometre measures of bw, sp
% This makes it easier to see which data is from which run
outname = [outfilestem, '_bw', num2str(bw), 's', num2str(sp), outfiletag];
outfile = fullfile(outpath, outname);

% Create a csv with this name and write out the column headers
% makeheaders is a helper function that concatenates the variable name with
% the statistic type (coefficient, p value etc.) to make the headers for
% the output csv. It also writes the fixed headers for model stats.
makeheaders(filenamearray, ncheckval, outfile);

% next we convert bandwidth and spacing from km to number of cells:
bw = bw * 1000 / cellx; % cellx = celly, so can use either
sp = sp * 1000 / cellx;

% Next we see how many subarrays will fit into our data arrays on each axis
% This tells how many rows and columns of subarrays we'll have
% Any data points in a strip of width bw-1 cells along the right and bottom 
% of the main arrays will be ignored. In practice, the impact of this will 
% be very small- the dataset on the map is an irregular shape, so most of 
% the edge of the map is NoData anyway.

[gridY, gridX] = size(YXarr{1}); % dimensions of the data arrays
n_cols = floor((gridX-bw)/sp)+1 % n subarrays in x direction
n_rows = floor((gridY-bw)/sp)+1 % n subarrays in y direction
% Note that for these n subarray calcs it's necessary first to subtract the
% bandwidth from the grid dimensions- otherwise, when the top left of the
% subarray gets close to the bottom/right of the map, the bottom/right of 
% the focal window is off the edge.

% Now we make counters to help us iterate
i = 1; % counter for x direction
j = 1; % counter for y direction
while j <= n_rows % go row by row till we hit the target number of rows
    % subx_, suby_, sp and bw are all numbers of grid cells on the map, and
    % our arrays of datapoints directly reflect the grids of cells, so we
    % can now use subx_, suby_ etc. to index into these arrays.
    suby_top = (j-1)*sp + 1;       % get index of top and bottom of row
    suby_bottom = suby_top + bw - 1;
    while i <= n_cols % go col by col up to target number of cols in row
        subx_left = (i-1)* sp + 1; % get left and right limits of column
        subx_right = subx_left + bw - 1;
        % get subarrays from each main array:
        subarrs = {nvars}; % make cell array to hold subarrays
        for n = 1:nvars
            vararr = YXarr{n}; % get main array for 1st,2nd,...,nvars-th 
            % variable. Index into it to get nth subarray (corresponding 
            % to nth variable, where n=1 gives Y array, n=1 gives first X
            % array, etc. etc.) and add that to the cell array:
            subarrs{n} = vararr(suby_top:suby_bottom,subx_left:subx_right);
            if n==1
                outnameM = ['j', num2str(j-1), 'i', num2str(i-1), '_PVbar_ml1.csv'];
                outpathM = fullfile(outpath, outnameM);
                csvwrite(outpathM,subarrs{n})
            end
            
            if windowtest == 1
                formatSpec = 'printing subarr for j (row) = %d; i (col) = %d; nvar = %d';
                str = sprintf(formatSpec,j,i,n)
                subarrs{n}
            end
        end
        
%       Before processing the subarrays to run the regression, we check
%       how many NoData values (-9999) there are. Because the map data is
%       preprocessed to make all arrays perfectly coterminous, we only 
%       have to check one, except when using log SEG data- that array is 
%       more limited, as explained above, so for log runs we have to use
%       the first (Y) array.
        nNoData = sum(subarrs{1}(:) == -9999); % number of NoData vals
        if nNoData > bw * bw * NaNmax          % acceptable limit
            i = i + 1;
%             disp("too much NoData; skipping window")
            continue
        end
        % If the window was too NoData-heavy, it will now have been
        % discounted. If it's still being processed, we know there aren't
        % too many NoData values, but there may still be some, so we check.
        % If we find any, we replace them with NaN in all subarrs:
        
        if nNoData > 0                      % check for NoData vals
            for n = 1:nvars
                subarr = subarrs{n};        % get the nth subarray
                subarr(subarr==-9999)=NaN;  % set -9999 vals to NaN
                subarrs{n} = subarr;        % put the subarray back in
            end
        end
        % Now we turn our subarray into a vertical array with one column
        % per variable and one row per mapcell
        % stackcols takes an n by m array and stacks all the columns on top 
        % of each other, giving an nm by 1 array. We do this with the
        % subarray of Y data first:
        subarrs_in = stackcols(subarrs{1});   % array with col of Y vals
        % and then for each remaining var...
        for n = 2:nvars                     
            subarr = subarrs{n};              % ...we take the subarray...
            vsubarr = stackcols(subarr);      % ...stack it vertically...
            subarrs_in = horzcat(subarrs_in, vsubarr); % ...add col to array
        end
        % subarrs_in is now nm high and nvars wide
        
        if windowtest == 1
            formatSpec = 'printing subarrs_in (equiv. subarrsDF) for j (row) = %d; i (col) = %d';
            str = sprintf(formatSpec,j,i)
            subarrs_in
%             comparisonOutfile = 'D:\UKuni\3rdYr\Project\GIS_data\DissPy\Regression\MWRtoPy\window_rvals_NRrows_3.csv';
%             dlmwrite(comparisonOutfile,subarrs_in,'delimiter',',','-append','precision',8)
%             subarrs_in = rowShuffle(subarrs_in);
%             comparisonOutfile = 'D:\UKuni\3rdYr\Project\GIS_data\DissPy\Regression\MWRtoPy\window_rvals_Rrows_3.csv';
%             dlmwrite(comparisonOutfile,subarrs_in,'delimiter',',','-append','precision',8)
        end
        
        % Now we get the lat and lon coords of the centre of our subarray 
        % so we can locate it on the map again. Mean of edge vals gets us 
        % the centre; -1 corrects for MatLab's 1-indexing; timesing by cell 
        % dims converts to distance on map. 
        fprintf('i (col) = %d; j (row) = %d\n',i,j)
       
        meanLeftRight = mean([subx_left, subx_right]);        
        fprintf('subx_left = %d; subx_right = %d. meanLeftRight = %f.\n',subx_left,subx_right,meanLeftRight)
            
        lonToAdd = (mean([subx_left, subx_right]) - 1) * cellx;
        lon = refx + (mean([subx_left, subx_right]) - 1) * cellx;
        fprintf('lonToAdd = %d; lon = %f\n',lonToAdd,lon)
     
        meanTopBottom = mean([suby_top, suby_bottom]);
        fprintf('suby_top = %d; suby_bottom = %d. meanTopBottom = %f.\n',suby_top,suby_bottom,meanTopBottom)
        
        latToSubtract = (mean([suby_top, suby_bottom]) - 1) * celly;
        lat = refy - (mean([suby_top, suby_bottom]) - 1) * celly;
        fprintf('latToSubtract = %d; lat = %f\n',latToSubtract,lat)
        

        % checkval is a reference value used during dev to check the mapped 
        % values are lining up correctly. If they are, this reference value
        % should be the same as the original value in the layer on the 
        % map in the cell where the generated datapoint lands.
        ref_in_sub = floor(bw/2);
        subtocheck = subarrs{ncheckval};
        checkval = subtocheck(ref_in_sub,ref_in_sub);
        
        % Now we run the regression on our selected data by calling the 
        % regressout function that runs the regression and writes the data 
        % to csv for mapping (see below)
        regressout(checkval, lon, lat, subarrs_in, outfile);
        i = i + 1; % increment x counter to move to next column within row
    end
    % When above loop finishes, we're done with one row & move to the next
    j = j + 1; % increment y counter to move to next row
    i = 1;     % reset x counter to move to beginning of row
end
t_end = now;
% print how long it took (86400 converts from days to seconds)
runtime = (t_end-t_start) * 86400
