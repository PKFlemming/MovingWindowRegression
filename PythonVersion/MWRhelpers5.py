import statsmodels.api as sm
import numpy as np
import csv
from math import sqrt, floor

def dictrounder(d, roundto):
    if roundto == None: # if there's no rounding to apply, just return the original dict
        outdict = d
    else:
        outdict = {k : f"{d[k]:.{roundto}g}" for k in d}
    return(outdict)

def getDfSegment(df, segment):
    nRows = len(df.index)
    firstHalfN = floor(nRows / 2)
    secondHalfN = nRows - firstHalfN
    if segment.lower() in ["n", "north"]:
        dfSegment = df.head(firstHalfN)
    elif segment.lower() in ["s", "south"]:
        dfSegment = df.tail(secondHalfN)
    elif segment.lower() in ["a", "all"]:
        dfSegment = df                      # this shouldn't come up bc if segment == all we don't need to segment,
                                            # but just in case...
    print(f"using {segment} segment of data")
    return dfSegment

# to make this usable for full, segment and MW regression, it needs to take df, not path+x+y+segment.
# assume first column is y
# segment df before inputting, if needed
def backStep(y, X, thold=0.05, roundto=None,verbosity=None):
    # roundto sets rounding for display of coefs and pvals. If roundto == None, no rounding is applied
    # rounding is done with string formatting, so rounded vals will display as 'strings'

    # df = pd.read_csv(dataPath)
    # if segment.lower() in ["n", "s", "north", "south"]:
    #     df = getDfSegment(df, segment)

    varsToRemove = True # we assume there are some variables whose pvals > thold (up until we see that there are not)
    anyVarsRemoved = False      # to enable us to carry things over between loops we initialise this marker
                                # we'll flip it if we remove some variables
    xsdictToCarryOver = {}      # These are the dicts we'll use to carry things over between loops
    xsdictCarriedOver = {}
    while varsToRemove:
        Xwconst = sm.add_constant(X) # TODO: this goes above loop? No, bc needs to change each time as xs removed?
        # adds a constant (column of 1s) to the array of xs
        xvars = list(X.columns)      # make list of X column names (x vars)
        xswconst = ["const"] + xvars # prepend "const" to list of var names

        OLSregression = sm.OLS(y, Xwconst) # these two lines run the OLS regression and get the results
        OLSresults = OLSregression.fit()

        # get the coefficients, standard errors, t statistics and pvalues out of the OLS results
        bs, ses, ts, pvals = OLSresults.params, OLSresults.bse, OLSresults.tvalues, OLSresults.pvalues
        if verbosity == "verbose":
            print(OLSresults.summary())

        if anyVarsRemoved == True:
            xsdict = xsdictCarriedOver
            # iff we're on nth loop, n>1 (1-indexed), we've removed some vars, and anyVarsRemoved == True. This means
            # that we've removed at least one xvar from the input to the regression, but we still want to know the pval
            # for that xvar in the regression just to confirm that it's been removed intentionally. For that reason we
            # initialise the xdict as already containing that data on the removed variable(s)- the data carried over
            # from the previous loop. If anyVarsRemoved == False, we haven't removed any vars and are necessarily on the
            # first loop, so we just want an empty dict as our xsdict.
        else:
            xsdict = {}  # initialise empty dict to hold info on x vars
            # the keys of xsdict will be the xkeys; the values will be- for each x variable- a dict of statistics
        for xidx, xkey in enumerate(xswconst):
            b = bs[xidx]
            # populate dict of info on this x var
            xsdict[xkey] = {"b": b,                     # regression coefficient
                            "se": ses[xidx],            # standard error
                            "tval" : ts[xidx],          # t value
                            "pval": pvals[xidx],
                            }
        if verbosity == "verbose":
            for x in xsdict:
                print(x, ":", dictrounder(xsdict[x],roundto))       # TODO: make (recursive?) function?
            # print(OLSresults.summary())

        # We don't want to remove the constant from our regression, even if its pvalue is high- removing it throws off
        # everything else and gives unrealistic r values for the overall regression. For this reason we make a "pvalues
        # without constant" object over which to iterate. Const,pvalue is the first member of pvals, so we remove it:
        pvalsWOConst = pvals[1:]

        # if any remaining pval is greater than our threshold value, we remove the x variable with the highest pvalue
        # and go again
        if any([pval > thold for pval in pvalsWOConst]):  # TODO: all->any; modeldict-> end?
            for x in xsdict:
                # print(f"printing x: {x}")
                xparams = xsdict[x]
                xpval = xparams["pval"]
                if xpval == max(pvalsWOConst): # TODO: surely more efficient way to do this? "if ismax..." or sth. idx max? argmax?
                    X = X.drop(columns=x)
                    if verbosity == "verbose":
                        print(f"\nremoving variable '{x}' because its pval = {round(xpval,3)} > {thold}\n")
                    anyVarsRemoved = True
                    # now we populate the carry-over dict with vals we want to carry over to the next loop. We're
                    # removing an xvar bc its pval is above our threshold, so we want to know that pval, but we don't
                    # want to retain other info on it. Therefore other dict vals are set to nan.
                    xsdictToCarryOver[x]={"b" : np.nan,
                                          "se" : np.nan,
                                          "tval" : np.nan,
                                          "pval" : xpval}
                    xsdictCarriedOver = xsdictToCarryOver # could be condensed, but this makes the var names more explicit
        # if no pvals are greater than threshold, there are no variables to remove and the bs, pvals etc. we've generated
        # are our final model params
        else:
            varsToRemove = False
            modeldict = {"r2": OLSresults.rsquared,
                         "r2adj": OLSresults.rsquared_adj,
                         "fstat": OLSresults.fvalue,
                         "pval": OLSresults.f_pvalue,
                         "rmse": sqrt(OLSresults.mse_resid),
                         # of mse_model, _total, and _resid, this is the one that matches the MatLab output
                         }

    return(xsdict, modeldict)

###

def regressout(checkval, i, j, yX, outfile, thold = 0.05):
    # get y and X variables
    y = yX.iloc[:, 0]  # get y data from df (first column)
    X = yX.iloc[:, 1:]
    # run backwards stepwise regression. y is dependent variable, X is set of independent variables
    xsdict, modeldict = backStep(y, X, thold=thold, roundto=None,verbosity="None") #TODO: add toggle in driver?

    # Calculate how many columns we need to have for the data we're going to write out:
    nxs = len(X.columns)
    ncols_out = 11 + 8 * nxs # 11 cols of model data + 8 cols for each x var

    # Make array with this many columns (each window will be described by one row in the output csv)
    data_out = np.empty(ncols_out)

    # Write out coords and checkval
    # Coords are in metres, calculated to coincide with Arcmap's values
    # Grid system used here is Kalianpur 1975 (UTM 43N). See driver code for explanation of checkval
    data_out[0] = i  # lon
    data_out[1] = j  # lat
    data_out[2] = checkval

    ## Now we want to write out the model stats for this subarray regression, starting with the general model stats.
    ## We'll do stats for each independent variable later.
    data_out[3] = modeldict["rmse"]
    data_out[4] = modeldict["r2"]
    data_out[5] = modeldict["r2adj"]
    data_out[6] = modeldict["fstat"]
    data_out[7] = modeldict["pval"]
    data_out[8] = np.nan                        # SSresids goes here

    # Calculate and write out mean and stdev of y data
    data_out[9] = np.nanmean(y)
    print(y)
    print(np.nanmean(y))
    print(np.mean(y))
    stdevY = np.nanstd(y, ddof=1)
    # ddof is degrees of freedom. ddof=1 sets divisor to N-1 in stdev formula, as appropriate for stdev of a sample
    # rather than a whole population. In practice the difference is negligible: bandwidth used was generally of the
    # order of 100, meaning N >> 1000. At N = 1000, 1 > stdevP/stdevS > 0.999
    data_out[10] = stdevY

    # Now we loop through the model info for the independent variables
    # For the nth X variable, we want the nth member of each array of model stats
    xvars = list(X.columns)
    # print(f"xvars (X.columns) = {xvars}")

    for n in range(nxs):
        # cols 0-10 are holding general model info, so we start writing at column 11
        N = 11 + 8 * n
        xvals = X.iloc[:, n]                # x_n
        xvar = xvars[n]
        # print(f"xvar = {xvar}")
        xdict = xsdict[xvar]
        # print(f"xdict = {xdict}")
        b, se, tval, pval = xdict["b"], xdict["se"], xdict["tval"], xdict["pval"]
        stdevX = np.nanstd(xvals, ddof=1)

        data_out[N] = np.nanmean(xvals)  # line 74 of regressout_clean.m
        data_out[N+1] = stdevX
        data_out[N+2] = b           # estimated coefficient
        data_out[N+3] = se          # standard error
        data_out[N+4] = tval        # T statistic
        data_out[N+5] = pval
        bstar = b * stdevX/stdevY   # standardised regression coefficient
        data_out[N+6] = bstar
        data_out[N+7] = abs(bstar)     # abs of b* for comparing magnitudes

    # write data out to csv
    with open(outfile, 'a', newline='') as csvfile:
        rowwriter = csv.writer(csvfile, delimiter=',')
        rowwriter.writerow(data_out)


def txtToArray(txtfile):
    with open(txtfile) as file_name:
        txtAsArray = np.loadtxt(file_name)
    return txtAsArray

def makeheaders(vars, ncheckvar, outfile):
    # Function to generate headers for regression data csv written by regressout function. It creates the file and
    # writes the headers to the first line.

    PVtype = vars[0]  # first var is dependent variable Y, (log)PV density
    xvars = vars[1:]  # rest are independent vars X1, X2...Xn
    nxs = len(xvars)  # number of x vars

    coordheads = ["xcoord", "ycoord"]  # headers for coordinate columns
    checkvar = vars[ncheckvar]  # get name of checkvar (check variable)
    checkvarhead = checkvar + "_ref"  # make header for checkvar col

    # statsheads are the headers for the columns that contain the general model
    # statistics (as opposed to those that pertain to individual y and x vars)
    statsheads = ["RMSE", "r2", "adjr2", "mod_fstat", "mod_pval", "mod_ssresid"]
    # make headers for columns containing data on y var
    pvmeanhead = PVtype + "_mean"
    pvstdevhead = PVtype + "_stdev"
    pvheads = [pvmeanhead, pvstdevhead]
    # fixed headers are those about coords, the model, the checkval, & the dependent (PV) var- everything except x vars.
    fixedheads = coordheads + [checkvarhead] + statsheads + pvheads # all lists except cvh, so [cvh]

    # xtags are descriptive suffixes, one for each piece of informtion we're
    # going to write out about the x vars. For each x var we take (in same
    # order as in xtags cell array below) mean, standard deviation, estimated
    # regression coefficient, standard error, t statistic, p value,
    # standardised regression coefficient b*, and absolute value of b*.
    xtags = ["_mean", "_stdevX", "_b", "_st_error", "_t_stat", "_pval", "_b_star", "_b_star_abs"]
    nxtags = len(xtags)         # number of xtags per x var
    xheads = []                 # empty list

    for n in range(nxs):
        xvar = xvars[n];                    # for each x var...
        for m in range(nxtags):
            xtag = xtags[m]                 # ...for each xtag...
            xhead = xvar + xtag             # concatenate var + tag -> col head
            xheads.append(xhead)            # and add it to the array of column headers

    headers = fixedheads + xheads           # combine fixed headers and x var headers
    with open(outfile, 'w', newline='') as csvfile:
        headerwriter = csv.writer(csvfile, delimiter=',')
        headerwriter.writerow(headers)
