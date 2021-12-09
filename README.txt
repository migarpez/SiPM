This set of ROOT macros has been developed to characterize SiPMs.
They expect lists of root files or csv files as inputs.

///////////////////////////////
inverse.C
///////////////////////////////

This macro computes the breakdown voltage for each inverse IV curve
provided in .csv format. The input file is a list of .csv files (lists/
inverse.list).

It computes the (1/I)*(dI/dV) function for each curve at fits its maximum
to a landau function. The mean value of the fit is assumed to be the
breakdown voltage.

Change values VMIN, VMAX and RANGE to set a proper fit interval.

From time to time the fit doesn't work properly, so I suggest to review
each fit and do it manually using the fit panel if it did not perform
properly.

Final part of the macro is for drawing results, it should be modified
as needed.

///////////////////////////////
forward.C
///////////////////////////////

This macro computes the quenching resistance for each forward IV curve
provided in .csv format. The input file is a list of .csv files (lists/
forward.list).

The fit is performed between IMIN and IMAX values.

///////////////////////////////
wf_tt_tree.C
///////////////////////////////

This macro generates the desired root files from the originals .csv files.
The created root files containts a TTree called wf. Each entry of the tree
corresponds to a waveform and has three variables:

    -V[WFLENGTH]: a double array of size WFLENGTH with the recorded voltages
    of each waveform.
    -time[WFLENGTH]: a double array of size WFLENGTH with the recorded times
    of each waveform. Ploting V[]:time[] provides the waveforms.
    -wftime: a double value which represents the time delay of each waveform
    with respect to the previous waveform.

It expects a list of .csv files (lists/wf_tt.list) of 'wf' type. For each
file, e.g. 02_01_wf.csv, it looks for another file of name 02_01_tt.csv in
the same location. If it is not found, a default wftime of 2 us is assumed
for each waveform.

Finally, the root file is created in the same location as the original csv
files.

///////////////////////////////
cnoise_from_tree.C
///////////////////////////////

This macro computes all correlated noise variables (DCR, AP and XT) for each
.root file provided in the lists/cnoise.list file. It is quite chaotic.

Before a waveform in analyzed, it is determined if it is electronic noise or
it isn't. There is no general way to determine a waveform as noise, so the
algorithm should be modified depending on the situation.
