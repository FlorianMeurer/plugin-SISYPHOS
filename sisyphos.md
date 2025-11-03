# Welcome to SISYPHOS

<font color='red'>`**THIS MODULE IS IN DEVELOPMENT, ATM YOU CANNOT STOP A RUNNING PROCEDURE OTHER THAN KILLING OLEX2**`</font>
<br>
<br>

This tool is COOL AND  meant to provide means of serially refine multiple models against multiple datasets of the same structure.
In other words, it solves Sisyphos dilemma for you, repeating the same task for the same structure over and over again.
It can be used for finding out which crystal was the best out of many (looking at synchrotron users!), or for benchmarking theoretical calculations.

<br>

NoSpherA2 (* Kleemiss et al., Chem. Sci., 2021, 12, 1675
&nbsp; URL[https://pubs.rsc.org/en/content/articlehtml/2021/sc/d0sc05526c,PAPER] ) is implemented and can be used in this procedure.

<br>
<br>

Anomalous dispersion refinement (* Meurer et al., IUCrJ, 2022, 9
&nbsp; URL[https://journals.iucr.org/m/issues/2022/05/00/lt5050/,PAPER]) can be also used, especially when in the vicinity of absorption edges (see Anomalous Dispersion in Info tab).

# Inputs

## Set working dir

Hallo
Here you can show the program the location of your measured data. SISYPHOS will look for any .hkl that is in any subfolder here, except for output folders of SISYPHOS.

## Set solution ins

Here you can give a starting model in form of an .ins file. Further filetypes such as .cif and .res will be implemented in the future. The model WILL be refined for every dataset, so this is really just the starting point.

# Adjustment

## Adjustment / eV

This setting is only for synchrotron users who know an offset between the enery / wavelength in their headers and the exact energy / wavelength of their experiment. Leave this empty if you measured on an in-house difffractometer or you don't care about minute differences in the wavelength.

# Resolution

## Resolution / Angstroms

Here you can apply a general cutoff in d-spacing for all of your datasets.

# Modes

## Use NoSpherA2

Decide if you want to use NoSpherA2 to perform a non-spherical atomic formfactor calculation. Further details see NoSpherA2 (* Kleemiss et al., Chem. Sci., 2021, 12, 1675
&nbsp; URL[https://pubs.rsc.org/en/content/articlehtml/2021/sc/d0sc05526c,PAPER] ).

## Benchmark mode

This mode is can be used to test e.g. the optimal settings for HAR for your structure. You have to prepare a benchmarkfile, which contains the settings for NoSpherA2 using the ORCA5 quantum chemical software. Each line in the benchmark file will create a new model, with the calculation taking semi-colon seperated keywords (keyword:setting;next_keyword:nextsetting) for the calculation. Optional keywords are: "basis_name","method",
                  "ncpus", "mem", "charge",
                  "multiplicity", "full_HAR",
                  "Max_HAR_Cycles","becke_accuracy",
                  "Relativistic", "h_aniso",
                  "h_afix",
                  "cluster_radius", "DIIS",
                  "cluster_grow", "ORCA_SCF_Conv",
                  "ORCA_SCF_Strategy", "ORCA_Solvation", "ORCA_DAMP"
                  and "pySCF_Damping".

## DISP refinement

Decide if you want to refine anomalous dispersion parameters for a present element / position. Further details see DISP Refinement (* Meurer et al., IUCrJ, 2022, 9
&nbsp; URL[https://journals.iucr.org/m/issues/2022/05/00/lt5050/,PAPER]).

## Update weight

This is a switch to automatically update the weights used.

# Elements

Decide which Elements the Anomalous Dispersion values should be refined for.

# disp

If you have more than one position for one of the Elements given, you can choose if they should have the same DISP value (constrained) or an individual one for each site.

# Energy and wavelength sources

Here you can set up, where SISYPHOS will search for alternative wavelengths for your experiments. You can choose the file header with format NAME_ENERGY or if you already have an .ins file containing this information (has to have same name as the data block) you can also use this instead.
