% This is the backend MWR code: the bit that actually runs the regression
% on the data selected by the driver code, and writes the output to csv.

function [data_out] = regressout(checkval, i, j, subarrs_in, outfile)

pthresh = 0.05; % critical p value
% subarrs_in contains all the input data to be regressed
Ys = subarrs_in(:,1);                    % first column is Y vals
dims_subarrs_in = size(subarrs_in);      % get dimensions
ncols = dims_subarrs_in(2);              % get number of columns
Xs = subarrs_in(:,2:ncols);   % remaining columns are X1, X2... Xn vals

% run regression on Xs and Y
[b,se,pval,inmodel,stats] = stepwisefit(Xs, Ys);
% each of these outputs is an array
% coefficients, standard errors, p value, input (which we ignore), model
% stats

% make an array to hold the data we're going to write out
nxs = ncols - 1;            % 1 y col, rest are x cols
ncols_out = 11 + 8*nxs;     % 11 cols of model data; 8 cols on each x var
data_out = NaN(1,ncols_out);% make array of NaNs with this shape

% Write out coords and checkval
% Coords are in metres, calculated to coincide with Arcmap's values
% Grid system in use here is Kalianpur 1975 (UTM 43N). See driver code for
% explanation of checkval
data_out(1) = i;    % lon
data_out(2) = j;    % lat
data_out(3) = checkval;

% Now we want to write out the model stats for this subarray regression,
% starting with the general model stats. We'll do stats for each ind var
% later.
% root mean square error
rmse = stats.rmse;  
data_out(4) = rmse;

% r^2 and adjusted r^2
% source of formulae is https://uk.mathworks.com/matlabcentral/answers/
% 93200-how-can-i-obtain-the-r-squared-and-adjusted-r-squared-values-from-
% stepwisefit-in-the-statistics-tool
nY = length(Ys);
varY = var(Ys, 'omitnan');
degfree = stats.df0;
r2 = 1 - (rmse^2/varY) * ((nY-1-degfree)/(nY-1));
adjr2 = 1 - rmse^2/varY;
data_out(5) = r2;
data_out(6) = adjr2;

% assorted other model stats
data_out(7) = stats.fstat;      % F statistic
data_out(8) = stats.pval;       % p value
data_out(9) = stats.SSresid;    % sum of squares of residuals

% calculate and write out mean, stdev of Y data
meanY = mean(Ys, 'omitnan');
stdevY = std(Ys, 'omitnan');
data_out(10) = meanY;
data_out(11) = stdevY;

% now we loop through the model info for the independent variables

for n = 1:nxs
    xvals = Xs(:,n);
    N = 12 + 8*(n-1); % first 11 cols are holding general model info
    
    data_out(N) = mean(xvals, 'omitnan'); % mean of Xn in subarray
    stdevX = std(xvals, 'omitnan');
        % note: I've used sample std here, rather than population std, but
        % the difference is negligible. For a given set of input values,
        % stdev(population)/stdev(sample) = sqrt((N-1)/N), which tends to 1
        % as N becomes large.
    data_out(N+1) = stdevX;

    pvalue = pval(n);
    data_out(N+5) = pvalue; % it's a bit ugly having the N+5 term 
                            % here, but my output format was already
                            % pretty standardised before I realised I
                            % needed to exclude data based on pval
    % if the pvalue is greater than threshold, the ind var will not be
    % included in the final model, and we don't want to plot info about
    % its coefficient
    if pvalue > pthresh
        continue
    end
    data_out(N+2) = b(n);           % estimated coefficient
    data_out(N+3) = se(n);          % standard error
    data_out(N+4) = stats.TSTAT(n); % T statistic
    bstar = b(n) * stdevX/stdevY;   % standardised regression coefficient
    data_out(N+6) = bstar;
    data_out(N+7) = abs(bstar);     % abs of b* for comparing magnitudes
end
% write data out to csv
% dlmwrite(outfile,data_out,'delimiter',',','-append','precision',8);
