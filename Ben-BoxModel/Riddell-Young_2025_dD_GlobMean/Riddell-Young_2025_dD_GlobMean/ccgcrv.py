#!/usr/bin/env python
"""
Program for applying curve fitting/filtering of time series data and
printing results.

Takes an input file containing two columns of data, a decimal date and value.
Applies curve fitting algorithm to the data, and depending on options,
prints results to stdout (default) or to files if specified.

General python requirements:
	numpy, scipy, dateutil

ccgg python requirements:
	ccg_filter - curve fitting/filtering
	ccg_dates - date conversion routines
"""
from __future__ import print_function

import sys
import os
#import argparse

from dateutil.parser import parse
from dateutil.rrule import rrule, DAILY

import ccg_filter
from ccg_dates import calendarDate, decimalDateFromDatetime, datetimeFromDecimalDate

##########################################################################
def check_export(options):
	""" Check if any options are set for filtered curve data. """

	a = options.orig \
               or options.func \
               or options.poly \
               or options.smooth \
               or options.trend \
               or options.detrend \
               or options.smcycle \
               or options.harm \
               or options.res \
               or options.smres \
               or options.trres \
               or options.ressm \
               or options.gr

	return a


##########################################################################
def export_data(options, filt):
	""" print out the curve data """

	if options.sample:
		if options.samplefile:
			try:
				fp = open(options.samplefile, "w")
			except IOError as e:
				sys.exit("Can't open file for writing. %s" % e)

		else:
			fp = sys.stdout

		export_dates(options, fp, filt, filt.xp)

	if options.equal or options.user:
		if options.file:
			try:
				fp = open(options.file, "w")
			except IOError as e:
				sys.exit("Can't open file for writing. %s" % e)

		else:
			fp = sys.stdout

		if options.equal:
			# Create a new list of dates at sample interval to give to export_dates().
			# Not quite the same as filt.xinterp because it takes into account leap years
			dates = rrule(DAILY, interval=int(filt.sampleinterval), dtstart=options.startdate, until=options.lastdate)
			xdates = [decimalDateFromDatetime(dt) for dt in dates]
			if xdates[-1] > filt.xp[-1]:
				xdates[-1] = filt.xp[-1]  # avoid problems with rounding and interpolation in ccgfilt
			if xdates[0] < filt.xp[0]:
				xdates[0] = filt.xp[0]  # avoid problems with rounding and interpolation in ccgfilt

		else:
			f = open(options.user)
			xdates = [float(line.split()[0]) for line in f]
			f.close()


		export_dates(options, fp, filt, xdates)


##########################################################################
def export_dates(options, fp, filt, x):
	""" Export (print) data to file pointer fp at dates given by x.
	The values to print are given as boolean flags in the export class.
	Some values can only be printed at sample dates, i.e. original data and residuals
	"""

	if options.showheader:
		export_header(options, fp)

	frmt = "%13.6e"

	h = filt.getHarmonicValue(x)	# harmonics
	p = filt.getPolyValue(x)	# poly
	s = filt.getSmoothValue(x)	# function + short term smoothing
	t = filt.getTrendValue(x)	# poly + long term smoothing
	g = filt.getGrowthRateValue(x)	# growth rate, derivative of trend
	f = filt.getFunctionValue(x)    # function, poly + harmonics

	for i, xp in enumerate(x):
		if options.cal:
			(yr, mon, dy, hr, mn, sec) = calendarDate(xp)
			if options.hour:
				print("%4d %02d %02d %2d" % (yr, mon, dy, hr), end='', file=fp)
			else:
				print("%4d %02d %02d" % (yr, mon, dy), end='', file=fp)
		else:
			print("%13.8f" % xp, end='', file=fp)

		if options.sample and options.orig:    print(frmt % filt.yp[i], end='', file=fp)
		if options.func:                       print(frmt % f[i], end='', file=fp)
		if options.poly:                       print(frmt % p[i], end='', file=fp)
		if options.smooth:                     print(frmt % s[i], end='', file=fp)
		if options.trend:                      print(frmt % (t[i]), end='', file=fp)
		if options.sample and options.detrend: print(frmt % (filt.yp[i] - t[i]), end='', file=fp)
		if options.smcycle:                    print(frmt % (s[i] - t[i]), end='', file=fp)
		if options.harm:                       print(frmt % (h[i]), end='', file=fp)
		if options.sample and options.res:     print(frmt % (filt.yp[i] - f[i]), end='', file=fp)
		if options.smres:                      print(frmt % (s[i] - f[i]), end='', file=fp)
		if options.trres:                      print(frmt % (t[i] - p[i]), end='', file=fp)
		if options.sample and options.ressm:   print(frmt % (filt.yp[i] - s[i]), end='', file=fp)
		if options.gr:                         print(frmt % (g[i]), end='', file=fp)

		print(file=fp)


##########################################################################
def export_header(options, fp):
	""" Export a line with column header names to file pointer fp.
	"""

	frmt = "%-13s"

	print(frmt % "date", end='', file=fp)

	# make sure these are in same order as in export_dates()
	if options.sample and options.orig:    print(frmt % "value", end='', file=fp)
	if options.func:                       print(frmt % "function", end='', file=fp)
	if options.poly:                       print(frmt % "polynomial", end='', file=fp)
	if options.smooth:                     print(frmt % "smooth", end='', file=fp)
	if options.trend:                      print(frmt % "trend", end='', file=fp)
	if options.sample and options.detrend: print(frmt % "detrended", end='', file=fp)
	if options.smcycle:                    print(frmt % "smooth_cycle", end='', file=fp)
	if options.harm:                       print(frmt % "harmonics", end='', file=fp)
	if options.sample and options.res:     print(frmt % "residuals", end='', file=fp)
	if options.smres:                      print(frmt % "smooth_resid", end='', file=fp)
	if options.trres:                      print(frmt % "trend_resid", end='', file=fp)
	if options.sample and options.ressm:   print(frmt % "resid_smooth", end='', file=fp)
	if options.gr:                         print(frmt % "growth_rate", end='', file=fp)

	print(file=fp)

#########################################################################
def read_data(filename=None):
	"""
	# Read in the input data file.
	# Format is always two columns,
	# the first column a decimal date value, (e.g. 2010.5 is halfway through 2010)
	# the second column is the corrsponding measurement value.
	"""

	if filename is None:
		fp = sys.stdin
	else:
		try:
			fp = open(filename)
		except IOError as e:
			sys.exit("Cannot open input file. %s" % e)


	x = []
	y = []
	for line in fp:
		(xv, yv) = line.split()
		x.append(float(xv))
		y.append(float(yv))

	fp.close()


	return x, y


#########################################################################

startdate = None

# manually set variables here instead of using argparse
class Options:
    def __init__(self):
        self.npoly = 5
        self.nharm = 4
        self.interv = 0
        self.short = 80
        self.long = 667
        self.gap = 0
        self.gain = False
        self.timez = None
        self.file = None
        self.samplefile = None
        self.equal = False
        self.sample = True
        self.cal = False
        self.hour = False
        self.date = None
        self.user = None
        self.showheader = False
        self.orig = False
        self.func = False
        self.poly = False
        self.smooth = False
        self.trend = False
        self.detrend = False
        self.smcycle = False
        self.harm = False
        self.res = False
        self.smres = False
        self.trres = False
        self.ressm = False
        self.gr = False
        self.coef = None
        self.stats = False
        self.amp = False
        self.mm = False
        self.annual = False
        self.args = ['data/ALT_MPI_dD.txt']

options = Options()



if options.npoly < 0 or options.npoly > 10:
	sys.exit("Error in --npoly argument: value out of range (0-10) %s" % options.npoly)

if options.nharm < 0 or options.nharm > 10:
	sys.exit("Error in --npoly argument: value out of range (0-10) %s" % options.nharm)

if options.interv < 0:
	sys.exit("Error in --interval argument: value out of range (>=0) %s" % options.interv)

if options.short < 0:
	sys.exit("Error in --short argument: value out of range (must be >=0) %s" % options.short)

if options.long < 0:
	sys.exit("Error in --long argument: value out of range ( must be >=0) %s" % options.long)

if options.gap < 0:
	sys.exit("Error in --gap argument: value out of range (must be >=0) %s" % options.gap)

if options.date:
	try:
		startdate = parse(options.date)
	except ValueError as err:
		sys.exit("Can not get valid date from --date argument '%s': %s" % (options.date, err))



if options.coef:
	try:
		(begcoef, endcoef) = options.coef.split(",")
		begcoef = int(begcoef)
		endcoef = int(endcoef)
	except ValueError:
		sys.exit("Cannot get coefficient range.")


args = options.args
if not len(args):
	xp, yp = read_data()
else:
	inputfile = args[0]
	xp, yp = read_data(inputfile)


# if user dates or equal spaced dates aren't specified, use sample dates as default
if not options.user and not options.equal: options.sample = True

# Compute the filtered data
if options.timez is None: options.timez = int(xp[0])
filt = ccg_filter.ccgFilter(xp, yp, options.short, options.long, options.interv, options.npoly, options.nharm, options.timez, options.gap, options.gain)


# If starting date is not specified, set it to the date of the first data point
# Set ending date to date of last data point
if startdate is None:
	options.startdate = datetimeFromDecimalDate(filt.xp[0])
else:
	options.startdate = startdate
options.lastdate = datetimeFromDecimalDate(filt.xp[-1])

if check_export(options):
	export_data(options, filt)

if options.amp:
	months = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
	amps = filt.getAmplitudes()

	print(" *****  Seasonal Cycle Statistics.  *****")
	print(" Year      Amplitude     Maximum   Date     Minimum   Date")
	print("-----------------------------------------------------------")

	frmt = "%5.0f %12.2f %12.2f   %3s %2d %9.2f   %3s %2d"
	for (year, amp, maxdate, maxval, mindate, minval) in amps:

		(yr, mnmax, dmax, hr, mn, sec) = calendarDate(maxdate)
		(yr, mnmin, dmin, hr, mn, sec) = calendarDate(mindate)

		print(frmt % (year, amp, maxval, months[mnmax], dmax, minval, months[mnmin], dmin))

if options.stats:
	print(filt.stats())

if options.mm:
	mm = filt.getMonthlyMeans()
	for (year, month, val, std, n) in mm:
		print("%4d %02d %7.2f %5.2f %2d" % (year, month, val, std, n))

if options.annual:
	am = filt.getAnnualMeans()
	for (year, val, std, n) in am:
		print("%4d %7.2f %5.2f %2d" % (year, val, std, n))

if options.coef:
	for i in range(begcoef, min(filt.numpm, endcoef+1)):
		print(" %.6f" % filt.params[i], end='')
#		print("%d %.6f" % (i, filt.params[i]))
	print()


# Save output

# specify the output directory and file name in a writable location
output_path = os.path.expanduser('output/filename.txt')  # saves the file to your home directory

# create the directory if it doesn't exist
os.makedirs(os.path.dirname(output_path), exist_ok=True)

# open the file for writing
with open(output_path, 'w') as f:
    sys.stdout = f  # redirect stdout to the file

    # Perform all logic inside the context where stdout is redirected

    if options.timez is None: 
        options.timez = int(xp[0])
    filt = ccg_filter.ccgFilter(xp, yp, options.short, options.long, options.interv, 
                                options.npoly, options.nharm, options.timez, options.gap, options.gain)

    # Call export functions within the redirected stdout context
    if check_export(options):
        export_data(options, filt)

    if options.amp:
        months = ["", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
        amps = filt.getAmplitudes()

        print(" *****  Seasonal Cycle Statistics.  *****")
        print(" Year      Amplitude     Maximum   Date     Minimum   Date")
        print("-----------------------------------------------------------")

        frmt = "%5.0f %12.2f %12.2f   %3s %2d %9.2f   %3s %2d"
        for (year, amp, maxdate, maxval, mindate, minval) in amps:
            (yr, mnmax, dmax, hr, mn, sec) = calendarDate(maxdate)
            (yr, mnmin, dmin, hr, mn, sec) = calendarDate(mindate)
            print(frmt % (year, amp, maxval, months[mnmax], dmax, minval, months[mnmin], dmin))

    if options.stats:
        print(filt.stats())

    if options.mm:
        mm = filt.getMonthlyMeans()
        for (year, month, val, std, n) in mm:
            print("%4d %02d %7.2f %5.2f %2d" % (year, month, val, std, n))

    if options.annual:
        am = filt.getAnnualMeans()
        for (year, val, std, n) in am:
            print("%4d %7.2f %5.2f %2d" % (year, val, std, n))

    if options.coef:
        for i in range(begcoef, min(filt.numpm, endcoef + 1)):
            print(" %.6f" % filt.params[i], end='')
        print()

# restore stdout to terminal
sys.stdout = sys.__stdout__

