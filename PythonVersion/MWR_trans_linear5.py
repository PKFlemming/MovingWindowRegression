# Front end of MWR. Driver code. Equivalent to GDR_2_2_clean (internal dev reference)

"""
Purpose:
Assess spatial variation in regression coefficients. Originally the regressions explored the relationship between
density of solar energy generation (MW of installed capacity per unit area) and various possible siting criteria such as
solar irradiance, population density and atmospheric dust levels
Mechanism:
MWR (Moving Window Regression) is like a simpler version of GWR (Geographically Weighted Regression). It takes map data
for a dependent variable Y and some number of independent variables X1, X2...Xn. Data format is described below. The
arrays in have the shape of the map data, i.e. each array could be mapped directly onto the map using its i,j
coordinates in the array as lon and lat (with appropriate scaling and translation constants).
User specifies a bandwidth
bw and a spacing sp. The code then scrolls through the Y and X1...Xn arrays, starting in the top left corner and taking
for each variable a square subarray corresponding to a map "window" of edge length bw. It runs a linear regression on
these subarrays and writes the regression data (r2, pvalue, coefficients etc.) out to csv. It then moves to the right by
a distance equal to the specified spacing and takes a new set of subarrays (the next window). It continues along the row
until it reaches the end, then moves down by sp and starts again at the beginning of the next row.

Note: Originally written in MatLab, and I'm still in the process of tweaking the translation, so there are a bunch of print
calls and #TODOs that can be removed. I'll upload a clean version when it's ready.
"""


from MWRhelpers5 import txtToArray, makeheaders, regressout
from math import floor
import numpy as np
from os.path import join
import pandas as pd

### INPUT BLOCK ###
# enter file path info for input/output data
inOutPath = r"D:\UKuni\3rdYr\Project\GIS_data\DissPy\Regression\resampled_asciis\snap_15_500_full"

# filenamearray is an array of the files to use as input for the regressions. The first is the y (dependent) variable.
# all the rest are candidate x (independent) variables. Files are all in the folder pointed to by inOutPath.
# Variables to be excluded are commented out.
filenamearray = [
    'pv',  # first
    # 'L10pv',
    'GTI',
    'aodb',
    # 'T',
    'slope',
    'load',
    'pl_dens',
    'L10PCI',
#     'L10TDP',
#     'ed_rds_dv',
    'logpop_fs20',
    # 'logpop_fs30',
#     'dummy_in_R1', # used to check windowing during translation from matlab code
#     'dummy_in_R2',
#     'dummy_in_R3'
    ]

windowtest = 0 # used to toggle values for python translation testing

batchtag = "_15_500_full.txt" # invariant tag for each group of input files
if windowtest == 1:
    batchtag = '.txt'; # used to check windowing during translation from matlab code

ncheckvar = 1 # which of the vars you want exported as a checkval to doublecheck the imported points are in the right
              # place on the map. Y is 0, X1 is 1, Xn is n

# path to write the data out
outpath = r"D:\UKuni\3rdYr\Project\GIS_data\DissPy\Regression\resampled_asciis\snap_15_500_full\GDR_logpv"

# first bit of the outfile name, specified by user at each run to identify output files:
outfilestem = "linPVb_GASLPldlPcilPFS20"
if windowtest == 1:
    outfilestem = "WT"

# outfiletag could just be .csv, or can include an additional version identifier ('b.csv'; 'version_2.csv') if running
# slightly modified versions of the same basic regression
outfiletag = '_py6_kk_OW_LLx.csv'

# dimensions, in meters, of map cells:
[cellx, celly] = [500,500]
if windowtest == 1:
    [cellx, celly] = [1000, 1000]

# coords of top left corner of map (get from Arc):
[refx, refy] = [398334.8956, 2032383.292]

# we're going to create subarrays (windows) by slicing into our main arrays
# first we specify the bandwidth (edge length) of our subarrays, in km:
bw = 200
# next we specify the spacing, in km, that we want between our focal points (the subarray centres). If sp = bw, there
# will be no overlap, and each map cell will be counted only once. sp << bw -> lots of overlap -> focus will move more
# slowly over the map, producing a more smoothly varying output.
sp = 50

if windowtest == 1:
    bw = 4
    sp = 4

# Now we set the maximum fraction (0-1) of a subarray that can be NaN before we reject the subarray. Tried playing with
# this but for non-log SEG 0.5 seems best. If it's higher you get centres plotted that are off the edge of the mapped data.
# With logPV data it *needs* to be higher though- there's gaps in the data bc all 0 vals -> NoData. NaNmax = 0.6~0.7 works at
# bw = 150. Presumably the further one zooms in the lower NaNmax can be.
NaNmax = 0.7

#### END OF USER INPUT ####

## this is where the timer went; haven't included that in the translation yet

#### READ BLOCK ####

# txtToArray takes a text file of the type outputted by Arc's raster-to-ascii function (element delimiter = space,
# line delimiter = newline) and turns it into an array
# this line runs a list comprehension: it takes each file, as pointed to by inOutPath+filename+batchtag, turns it into an
# array, and makes a list of all of the (2d) arrays
separateArrays = [txtToArray(join(inOutPath, filename + batchtag)) for filename in filenamearray]
# we then turn that back into a (3d) array
YXarr = np.array(separateArrays) # debug comment: array has same dims in py and matlab

# YXarr now holds all the data. YXarr has shape n by i by j, where n = nvars = number of variables (both Y and X) while
# i and j are the dimensions of the map (in gridcells, not metres)

#### RUN BLOCK ####
# First we make the output filepath using the kilometre measures of bw, sp
# This makes it easier to see which data is from which run
outname = outfilestem + "_bw" + str(bw) + "s" + str(sp) + outfiletag
outfile = join(outpath, outname)

# Create a csv with this name and write out the column headers
# makeheaders is a helper function that concatenates the variable name with
# the statistic type (coefficient, p value etc.) to make the headers for
# the output csv. It also writes the fixed headers for model stats.
makeheaders(filenamearray, ncheckvar, outfile)

### this (including the makeheaders helper) gets us up to line 113 in GDR_2_2_clean.m (internal dev ref) ###

# next we convert bandwidth and spacing from km to number of cells:
[bw, sp] = [val * 1000 / cellx for val in [bw, sp]]          # cellx = celly, so can use either
                                                             # debug comment: bw and sp same as matlab
# Next we see how many subarrays will fit into our data arrays on each axis
# This tells how many rows and columns of subarrays we'll have
# Any data points in a strip of width (gridX-bw)%sp cells along the right and height (gridY-bw)%sp cells along the
# bottom of the main arrays will be ignored. If the area of mapped data were rectangular this could introduce a bias, as we
# always have full data along the left and top and variably partial data along the right and bottom. In practice, the
# impact of will be small, as the dataset on the map is an irregular shape, so most of the edge of the map is
# NoData anyway. Could be slightly improved by centring.

YXarrDims = YXarr.shape
# print(f"\nprinting YXarrDims: \n{YXarrDims} \nend of YXarrDims\n")

[nvars, gridY, gridX] = YXarrDims[:] # dimensions of the data arrays
n_cols = floor((gridX - bw) / sp) + 1      # n subarrays in x direction # TODO: review this, incl. for odd bw
n_rows = floor((gridY - bw) / sp) + 1      # n subarrays in y direction

print(f"n_cols, n_rows = {n_cols, n_rows}") # debug comment: same n_cols, n_rows as matlab
                                            # also makes sense on paper, incl. for (gridX - bw) % sp == 0
# Note that for these n subarray calcs it's necessary first to subtract the bandwidth from the grid dimensions-
# otherwise, when the top left of the subarray gets close to the bottom/right of the map, the bottom/right of the focal
# window is off the edge.

# subx_, suby_, sp and bw are all numbers of grid cells on the map. Our arrays of datapoints directly reflect the grids
# of cells, so we can now use subx_, suby_ etc. to index into these arrays as if we were moving a window around the map
for j in range(n_rows):                      # go row by row till we hit the target number of rows TODO: n_rows+1?
    print(f"j={j}")
    suby_top = int(j * sp)                   # get index of top and bottom of row TODO: check for OBOEs
    suby_bottom = int(suby_top + bw)

    for i in range(n_cols):                  # go col by col up to target number of cols in row TODO: n_cols+1?
        print(f"i={i}")
        subx_left = int(i * sp)              # get left and right limits of column
        subx_right = int(subx_left + bw)
        # get subarrays from each main array:
        subarrs = []                    # empty list
        for n in range(nvars):
            vararr = YXarr[n]           # get main array for 1st,2nd,...,nth variable
            # index into it to get subarray and add that to the 3d array:
            subarr = vararr[suby_top:suby_bottom,subx_left:subx_right]
            # if n == 0:
            #     outnameM = "j"+str(j)+"i"+str(i)+"_PVbar_py.csv"
            #     outfileM = join(outpath, outnameM)
            #     np.savetxt(outfileM, subarr, delimiter=',')

            # print(f"n={n}")
            # print(f"subs={suby_top, suby_bottom, subx_left, subx_right}")
            subarrs.append(subarr)
            arrOfSubarrs = np.array(subarrs)
            if windowtest == 1:
                print(f"printing subarray for j (row) = {j}; i (col) = {i}; nvar = {n} (matlab: {j+1},{i+1}, {n+1})")
                print(subarr.astype(int))
            # print(arrOfSubarrs.shape)

        #### this gets us up to line 155 in GDR_2_2_clean.m (internal dev ref)####

        # Before processing the subarrays to run the regression, we check how many NoData values (-9999) there are.
        # Because the map data is preprocessed to make all arrays perfectly coterminous, we only have to check one,
        # except when using log SEG data- that array is more limited, as explained above, so for log runs we have to use
        # the first (Y) subarray. Default is therefore to use the y subarray.
        ySubarr = subarrs[0]
        nNoData = np.count_nonzero(ySubarr == -9999) # count no. of NoData values in window (= nNoData in subarray)
        if nNoData > bw * bw * NaNmax:               # acceptable limit
            continue
        # If the window was too NoData-heavy, it will now have been discounted. If it's still being processed, we know
        # there aren't too many NoData values, but there may still be some, so we check. If we find any, we replace them
        # with NaN in all subarrs:
        if nNoData > 0:                       # check for NoData vals
            for n in range(nvars):
                subarr = subarrs[n]                              # get the nth subarray
                subarr= np.where(subarr==-9999, np.nan, subarr)  # set -9999 vals to NaN
                subarrs[n] = subarr                              # put the subarray back in

        # Now we turn our 3d array (n by i by j) into a 2d array with one column per variable and one row per mapcell
        # first we take each of the n 2d subarrays (i by j) and turn it into a 2d array with 1 column (ij by 1)
        subarrs_in = [subarr.reshape(-1, 1) for subarr in subarrs]
        # TODO: speedcheck (orig. had subarrs overwritten here, but changed this bc need to get checkval out of subarr.)
        # might mean a lot of additional storage in memory
        # now we stack these horizontally into a 2d array (ij by n)
        hstackedSubarrs_in = np.hstack(subarrs_in)

        # note that, bc the windows are square, i = j = bw/(linear resolution in km)
        #### this gets us to line 191 in in GDR_2_2_clean.m (internal dev ref)####

        # now we make this a dataframe, bc currently that's what the backstep function is set up to accept as input
        # TODO: make more flexible/consistent/simple
        subarrsDF = pd.DataFrame(data = hstackedSubarrs_in,
                                 columns = filenamearray)
        # subarrsDF has been checked manually against subarrs_in in MatLab using dummy data, and is identical

        # Now we get the lat and lon coords of the centre of our subarray so we can locate it on the map again. Mean of
        # edge vals gets us central cell; -0.5 gets us to centre of cell; timesing by cell dims converts to distance on map.

        lon = refx + (np.mean([subx_left, subx_right]) - 0.5) * cellx  # TODO: check for odd bw
        lat = refy - (np.mean([suby_top, suby_bottom]) - 0.5) * celly

        print(f"j={j}, i={i} (1-indexed j={j+1}, i={i+1})")

        print(f"suby_top={suby_top}, suby_bottom={suby_bottom}. meanTopBottom={np.mean([suby_top, suby_bottom])}")
        print(f"latToSubtract = {(np.mean([suby_top, suby_bottom])-0.5) * celly}; lat = {lat}")

        print(f"subx_left={subx_left}, subx_right={subx_right}. meanLeftRight={np.mean([subx_left, subx_right])}")
        print(f"lonToAdd = {(np.mean([subx_left, subx_right])-0.5) * cellx}; lon={lon}")
        
        # checkval is a reference value used during dev to check the mapped values are lining up correctly. If they are,
        # this reference value should be the same as the original value in the layer on the map in the cell where the
        # generated datapoint lands. Checkval = check value; checkvar = check variable
        ref_in_sub = floor(bw/2)-1 # TODO: check for odd bw
        subtocheck = subarrs[ncheckvar]
        checkval = subtocheck[ref_in_sub,ref_in_sub]

        # remove nan values from dataframe
        subarrsDF = subarrsDF.dropna()      # TODO: is this ok???

        # print(f"\nprinting NaNs in y (I think?) \n{y.isna().sum()} \nend of NaNs in y \n")
        regressout(checkval, lon, lat, subarrsDF, outfile)
