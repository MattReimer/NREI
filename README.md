## Introduction

Drift‐foraging, a behavior by which fish capture drift‐born prey items from the water column, is a commonly observed feeding behavior among salmonid species. Salmonids are often observed holding positions in slower velocities while making feeding forays into faster waters with relatively higher prey deliver rates. It is believed that this behavior minimizes energetic expenditures due to swimming (by favoring positions sheltered from the fastest velocities when not actively capturing prey items) while simultaneously maximizing potential energy intake (by capturing and ingesting prey items from the higher prey delivery rates of adjacent faster waters).

Fish ecologists have attempted to model, quantitatively, the energetic tradeoffs of drift‐foraging. Such efforts typically estimate the energy gain from foraging and the energetic costs of swimming, and then they difference the two to estimate net rate of energy intake (NREI). NREI is usually expressed in units of energy per unit time, and it represents the net energy balance a drift‐forager could hope to obtain at a given foraging location per unit time of foraging effort. Fausch {}, Hughes and Dill {}, and Hill and Grossman {} pioneered studies linking NREI to foraging position selection in a variety of lotic settings, and others have since linked NREI to both habitat quality {} and capacity {}.

NREI models typically integrate multiple sources of information to produce NREI estimates. Characterization of depths and velocities near the foraging point of interest, information about the sizes and concentration (individuals/m3) of drifting invertebrates, and an estimate of fish length all serve to inform a foraging model. Velocity at the foraging point, water temperature, and an estimate of fish weight serve as inputs to a swim costs model that predicts the energetic costs of swimming at the foraging location. Because NREI modeling integrates a template of physical habitat (i.e., a map of depths and velocities derived from hydraulic modeling and a CHaMP DEM) with important biological variables (e.g., water temperature, prey concentration, and fish size) to characterize the energetic quality of habitat, NREI has potential as a tool to integrate multiple CHaMP data products in order to produce metrics related to growth potential and energetic habitat quality.

The Columbia Habitat Monitoring Program (CHaMP){} and The Integrated Status and Effectiveness Monitoring Program (ISEMP){} have funded an implementation of NREI modeling for use under the greater CHaMP/ISEMP research umbrella. The model is implemented in R {}, a free and open‐source statistical scripting language. The model uses CHaMP‐standard Delft3D outputs to characterize flow fields near potential foraging points throughout CHaMP reaches, and it has been customized to work with other CHaMP data products. The model allows multiple z foraging positions throughout the water column, which may be important in streams having deep pools {}, and it estimates NREI for many combinations of fish size, water temperature, and prey concentration from a single simulation (rather than one simulation for each unique combination of those inputs).

This document is intended as a user’s guide for distribution with the CHaMP/ISEMP NREI model files.

## Pre­simulation preparation 

### R and Rtools installation

You will need the free statistical software R installed on your computer to use the NREI model scripts. To install R, go to <https://cran.r‐project.org/> and install the latest version. You will also need to install Rtools, a set of tools designed for R developers. This is necessary so the NREI scripts can compress (i.e., zip) some of the larger model outputs to save space on your computer. To install R tools, go to <https://cran.r‐project.org/bin/windows/Rtools/> and follow the installation instructions. Other R packages necessary for NREI simulations (i.e., extra tool sets not included in the base distribution of R) will be installed when you run the model.

## Hard drive setup

1. Unzip nreiExampleFilesStart.zip
2. Paste the code, data_in, data_out, and simulations folders from the unzipped nreiExampleFilesStart folder onto your
C:/ drive (or another similarly high level location). Note that:
    1. The code/nrei directory is intended to contain the most up‐to‐date copies of the `runNreiMultiFish.R` and
    `nreiFunctionsMultiFish.R` scripts. If you wanted to permanently implement a change to the code for all future simulations, you would open the scripts here, make the change, and then save the new versions with an appropriate filename. Similar (superior) functionality could also be achieved using git or another version tracking software.
    2. The data_in directory is intended to house input data for the NREI model and other CHaMP/ISEMP analysis tools. Hydraulic model outputs that serve as inputs to NREI can be downloaded in compressed form (i.e., zipped) and later unzipped for use in NREI modeling. You will notice that the data_in folder contains a folder called champ, which contains folders for both zipped and unzipped input files. The provided data_in/champ/unzipData folder contains example inputs for your first NREI simulation. Information on populating the data_in folder with input data for other CHaMP visits can be found in Appendix XX.
    Note: If you already have a folder on your hard drive intended to hold
    Input data for similar CHaMP/ISEMP models (e.g. the HSI habitat model), you may be able to continue using your current folder structure. For now, however, please follow these instructions and use the provided example input files. By going through this user’s guide as‐written, you should be able to tell if your existing folder will meet the needs of the NREI model. You can later adjust the appropriate sections of `runNreiMultiFish.R` to work with your existing file structure.
    3. The data_out/nrei directory is intended to contain outputs from NREI simulations.
    4. The simulations/nrei directory is intended to hold simulation‐specific information for each of your NREI
    simulations, simulation sets, or other analyses. Typical NREI modeling process is to create a new, appropriately named folder inside the simulations/nrei directory each time you begin a new NREI analysis. Inside the new folder, you should also place an application‐specific version of the `runNreiMultiFish.R` script (copy it from the code folder), an input lookup table (this will be explained later), and any other simulation‐related notes that may be needed for later reference (e.g., methodology of your simulations). In depth analyses of the model itself may, from time to time, require copying the `nreiFunctionsMultiFish.R` script to this folder as well in order to make changes to the model’s functions, but this not necessary for routine NREI analyses.
    5. The remainder of this user’s guide will walk you through simulating NREI with these file structure guidelines in mind.


------------------

## Introduction to model scripts and input files

Before your first NREI simulation, it is a good idea to briefly familiarize yourself with the model’s scripts and input files.

### R scripts

The NREI model has two main scripts. The first script, `nreiFunctionsMultiFish.R`, contains the functions necessary for calculating NREI. The second script, `runNreiMultiFish.R`, reads inputs, estimates NREI, and then writes model outputs to file.

#### Input files

Each NREI simulation requires an input table. Each row of the required input table contains the inputs for a single NREI simulation. To see an example input lookup table, navigate to the directory simulations/firstNreiSim_26May2016 and open exampleInlookTab.csv in Excel. You will notice that exampleInlookTab.csv only has a single row (i.e., NREI will only be simulated for one CHaMP visit using this table), but users can batch‐process multiple CHaMP visits by creating similar input lookup tables with multiple rows. The columns VisitID, WatershedName, SiteName, VisitYear, Grad, Area_Wet, DpthBf_Avg, and SubD50 can be found in standard data exports available from champmonitoring.org. The columns est.grid.size, est.sim.duration.h, kroughness, grdreductionfactor, DZ, zreductionfactor are calculated using the Grad, Area_Wet, DpthBf_Avg, and SubD50 columns (described in more detail in Appendix XX). The kSpecies column is used to create a folder in the output file structure that corresponds to the species being modeled. Note that changing this value simply changes the name of a folder that will be inserted in the output file structure. Currently, it does not influence parameters or equations used in the model, although this functionality could be added in the future. The column kPreyLength contains the length of a representative drifting prey item (in meters) and, generally speaking, should not be changed. If you do need to change the kPreyLength value, you will also need to change corresponding lines of `runNreiMultiFish.R` (see Appendix XX). Further information on column names and their units as well as creating input lookup tables can be found in Appendix XX.

#### Hydraulic inputs

NREI simulations require depth‐averaged 2D hydraulic model output from Delft3D to characterize depths and velocities throughout CHaMP reaches. To see an example, navigate to the directory data_in/champ/unzipData/JohnDay/CBW05583‐ 240498/2012/VISIT_960 and open the dem_grid_results.csv file in Excel. Each row of dem_grid_results.csv contains select hydraulic model outputs for one wetted raster cell. During NREI simulations, the NREI model will read the contents of dem_grid_results.csv and use its data to interpolate depths and velocities near simulated foraging locations throughout the CHaMP reach being modeled. For information on how to populate the data_in folder with appropriate hydraulic outputs for other CHaMP visits, see Appendix XX.

## Simulating NREI

Once you have an appropriate input lookup table (provided with the example files), you are ready to begin using the model. To simulate NREI for the provided example files:

1. Copy and paste `runNreiMultiFish.R` from the code/nrei folder to the simulations/firstNreiSim_26May2016 folder.
2. Open `runNreiMultiFish.R`. Preferably, you’ll want a text editor with syntax highlighting for R scripts (e.g.,
Notepad++), but you could also open the script in the R console.
3. Make sure location pointers in the script match the file structure of your hard drive. If you pasted the code, data_in,
data_out, and simulations folders directly to your C:/ drive then no changes should be necessary (but you should still read the points below).
a. The variable inputs.dir (line 30) should point to the high level directory containing the unzipped Delft3D solution files (e.g., C:/data_in/champ/unzipData).
b. The variable outputs.dir (line 31) should point to the directory where model outputs are to be written. Note that if you provide a directory that doesn’t exist, it will be created for you upon running the NREI scripts.
c. The variable nrei.func.fn (line 32) should be the full file path of the `nreiFunctionsMultiFish.R` script on your hard drive.
d. The variable in.look.tab (line 54) should be the full file path of the input lookup table you wish to use. For the example files, it should point to exampleInLookTab.csv.
4. Adjust, if necessary, the number of processor cores to be used during simulation.
a. The variable sim.cores (line 88) indicates the number of processor cores to be used during simulation. While, in general, distributing the workload across multiple processor cores can decrease the time needed for simulation, there can be a point of diminishing returns and you cannot utilize more cores than your computer has. You can find instructions online to identify the number of processor cores your computer has. Choose a number based on your goals for simulation (i.e., speed vs. using your computer for other things while simulations are running).
5. Save any changes to `runNreiMultiFish.R`, and then open the script in R (if you haven’t already).
6. Press ctrl + a to select the entire script, then press ctrl + r to send the highlighted text to the R interpreter and begin
the simulation. Note: the script can also be run block‐by‐block if desired.
Note: Users may also wish to activate/deactivate certain model outputs or alter the default fish size, water temperature, or drifting prey concentration values in the model. See below for more information.

