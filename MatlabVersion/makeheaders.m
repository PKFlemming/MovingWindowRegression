function [headers] = makeheaders(vars, ncheckval, outfile)
%%%
% Function to generate headers for regression data csv written by
% regressout function. It creates the file and writes the headers to the
% first line.
%%%

PVtype = vars{1};       % first var is dependent variable Y, (log)PV density
xvars = {vars{2:end}};  % rest are independent vars X1, X2...Xn
dims_xs = size(xvars);  % dimensions of array of x vars
nxs = dims_xs(2);       % number of x vars

coordheads = {'xcoord','ycoord'}; % headers for coordinate columns
checkval = vars{ncheckval};       % get name of checkval var
checkvalhead = strcat(checkval, '_ref'); % make header for checkval col
% statsheads are the headers for the columns that contain the general model
% statistics (as opposed to those that pertain to individual y and x vars)
statsheads = {'RMSE' 'r2' 'adjr2' 'mod_fstat' 'mod_pval' 'mod_ssresid'};
% make headers for columns containing data on y var
pvmeanhead = strcat(PVtype, '_mean');
pvstdevhead = strcat(PVtype, '_stdev');
pvheads = {pvmeanhead, pvstdevhead};
% fixed headers are all those about coords, the model, the checkval, and
% the dependent (PV) var- everything except the x vars.
fixedheads = [coordheads checkvalhead statsheads pvheads];

% xtags are descriptive suffixes, one for each piece of informtion we're
% going to write out about the x vars. For each x var we take (in same
% order as in xtags cell array below) mean, standard deviation, estimated
% regression coefficient, standard error, t statistic, p value,
% standardised regression coefficient b*, and absolute value of b*.
xtags = {'_mean' '_stdevX' '_b' '_st_error' '_t_stat' '_pval' '_b_star' '_b_star_abs'};
dims_xtags = size(xtags);   % dimensions of xtags cell array
nxtags = dims_xtags(2);     % number of xtags per x var
nxheads = nxs * nxtags;     % number of x column headers = n(x vars) * 
                            % n(x headers per x var)
xheads = {nxheads};         % empty cell array, one space for each x header

for n = 1:nxs        
    xvar = xvars{n};                    % for each x var...
    index_n = (n-1) * nxtags;
    for m = 1:nxtags        
        xtag = xtags{m};                % ...for each xtag...
        xhead = strcat(xvar, xtag);     % concatenate var + tag -> col head
        index_head = index_n + m;       % get its position
        xheads{index_head} = xhead;     % and add it to the cell array of 
                                        % column headers
    end
end

headers = [fixedheads xheads]; % combine fixed headers and x var headers
commaHeader = [headers;repmat({','},1,numel(headers))]; % insert commaas
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader); % cHeader in text with commas for csv
% write headers to file
fid = fopen(outfile,'w'); 
fprintf(fid,'%s\n',textHeader);
fclose(fid);