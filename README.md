About the necessary Python 3.6+ compiler:

On Windows its possible to run python programms on Visual Studio Code (see   https://code.visualstudio.com/docs/setup/windows).

On Linux S.O. its easier to run python as on the lastest versions com Ubuntu, python is pre-installed and you can check its version with:

$ python --version

If you don't have python on your system its quite easy to install with the following:

$ sudo apt update
$ sudo apt install build-essential
$ sudo apt install python3.8

If you are having problems, there are multiple easy to follow guides like these: https://phoenixnap.com/kb/how-to-install-python-3-ubuntu#:~:text=How%20to%20Install%20Python%203%20on%20Ubuntu%2018.04,Code%20%28Latest%20Version%29%20Step%201%3A%20Update%20Local%20Repositories

https://blog.eldernode.com/install-python-3-ubuntu-20/

------------------------------------------------------------------------------------------

#Molecular-Dynamics related files.

The following was tested on both Windows and Linux S.O.:

Analysis_Plots_4.0.py - RMSD, RMSF and Radgyr data analysis. This code needs numpy and matplotlib modules with an python 3.6+ compiler (if you dont usually use Terminal, install Visual Studio Code, which runs on Windows and Linux). 

--------------------------------------
About Analysis_Plots_4.0.py execution:

$ python3 Analysis_Plots_4.0.py -h

    Welcome to Analysis plot 4.0:

    Copyright (C) 2021  Braga, B. C.
    This program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions; use option '-v' for details.


    Usage:
        -v or --version         prints current version, the python libraries needed and the data format expected.       

        -h or --help            prints this message.

        -type   Quantity of files used, for data comparison.
                one: Normal plot with one data file.
                two: Plots two data files with the same axis information (ex:RMSD for pH7 and pH8).
                four: Plots four data files with the same axis information (ex:RMSD for pH5, pH6, pH7, pH8).

        -anatp          Analysis type:
                rmsd;
                rmsf;
                radgyr.

        -i              input data file(s) with a name for the plot (separated by space).
                        Ex: -i ph7.00_rmsd.dat pH=7.00

        -stitle         (Valid only for type four) Title for comparison plot.

        -lblcrd (Valid only for type four) Label coords.

        -fram2time      (For anatp = rmsd or radgyr) Sets X-axis will change from frame to time (pico seconds). Must inform frame-step conversion.
                        Ex:-fram2time 10000

        -nanosec                (For anatp = rmsd or radgyr) Sets time intervals on X-axis to nanoseconds.

        -dpi            Sets the plot quality (recommended to use only with type=one). Default: 100


    Examples:
        $ python3 Analysis_Plots_4.0.py -type one -anatp rmsd -i ph7.00_rmsd.dat Production pH=7.00 -fram2time 10000 -dpi 120


        $ python3 Analysis_Plots_4.0.py -type two -anatp rmsf -i ph7.00_rmsf.dat Production pH=7.00 ph8.00_rmsf.dat Production pH=8.00


        $ python3 Analysis_Plots_4.0.py -type four -stitle Production -anatp radgyr -lblcrd (16.61,200) -i ph7.00_radgyr.dat pH=7.00 ph8.00_radgyr.dat pH=8.00 ph9.00_radgyr.dat pH=9.00 ph10.00_radgyr.dat pH=10.00 -fram2time 10000 -nanosec 
        

The following was tested only on Linux S.O.:

ASM.py - Molecular Dynamics simulation manager (tested on AMBER18 and AMBERtools18+). Automates the process of performing a molecular dynamics simulation using the AMBER/AMBERtools computational packages.

  Ps: AMBERtools is an open source package of AMBER and both can be found on https://ambermd.org/
  
This programm needs both Python 3.6+ and AMBER/AMBERtools installed on the running machine.

-----------------------
About ASM.py execution:

$ python3 ASM.py -h

    Welcome to Amber Simulation Manager 1.1:

    Copyright (C) 2021  Braga, B. C.
    This program comes with ABSOLUTELY NO WARRANTY; This is free software, and you are welcome to redistribute it under certain conditions; use option '-v' for details.


    Usage:
        -v or --version         Prints current version and its corrections relating previous versions.

        -h or --help            Prints this message.

        -i              Input PDB file.

        -g              Goal as MD or CpHMD. If CpHMD was chosen, one pH must be given right after this flag. Default: MD.
                        OBS: If you choose CpHMD and don't use flag 'res' the code will choose by default to tritate AS4 SER HIP.

        -phdset         (This option has priority over the pH of option -g) Defines the pH list (only if you'll do multiple pH CpHMD production). After this flag an initial pH, a final pH and an interval unit must be given in this order (Ex: -phdset 4.0 7.0 0.5). Which represents pH:[4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0].

        -res            Residues to tritate, for CpHMD (which must be informed right after this flag).

        -rdmut          Random active site mutation.
                        1:random aminoacid mutates to a similar one (eg. SER_160-THR). 
                        2:random aminoacid mutates to one not so similar (eg. SER_160-MET).

        -atsite         (Auxiliary option for -rdmut) Active site - for enzymes only. Must be informed right after this flag as eg: -atsite SER_160 HIS_237 ASP_206. Necessary only if a mutation will be done. Defaut: set as petase's active site.

        -rdydone        (Auxiliary option for -rdmut) Informs which mutation you're restricting from the "random" choice. Eg. -rdydone SER_160-MET ASP_206-GLY ASP_206-ILE.

        -mut            Mutation of choice. Eg: -mut SER_160-MET.
                                Obs: For the mutation options (-mut and -rdmut) a solvent must be used. If none given, water will be used by default.

        -mode           Duration of the simulation's stages (in time steps).
                        low: Annealing/Equilibration: 10**4 |  Production: 10**5
                        med: Annealing/Equilibration: 10**5 |  Production: 10**6
                        gpu-med: Annealing/Equilibration: 10**6 |  Production: 5*10**6
                        gpu-high: Annealing/Equilibration: 5*10**6 |  Production: 10**7
                        gpu-ultra: Annealing/Equilibration: 10**7 |  Production: 5*10**7

        -explicit       Explicit solvent will be used. One of the following options must be given after this flag:
                        water
                        methanol
                        chloroform
                        N-methyacetamide
                        urea

        -icyc           (Integer between 50 and 1000) Number of times information regarding energy, rmsd, rmsf, coordinates, etc, are saved for each simulation stage (Default: 200).
                                Obs: For -mode gpu-ultra, this is set between 800 and 1000.

        -arq            Choice in architecture
                        gpu (if chosen, the GPU id must be informed right after this flag)
                        sander
                        pmemd

        -mpi            Multicore run (the number of cores must be informed right after this flag). Default compiler: mpiexec.
                                Obs: Not implemented for -arq gpu.

        -cut            Nonbonded cutoff in angstrom (Default: 12).

        -prepstp        Stops the run after all preparations are complete (right before minimization).
                                Obs: The shell script, to run the simulation, will still be created for you.


    Examples:
        $ python3 ASM.py -i 1UBQ.pdb -arq sander -mpi 4 -mode low -g CpHMD 6 -phdset 6.0 7.0 0.5


        $ python3 ASM.py -i 6eqe.pdb -arq pmemd -mpi 4 -mode low -g CpHMD 6.5 -explicit water


        $ python3 ASM.py -i 6eqe.pdb -arq gpu 0 -mode gpu-med -g MD -explicit water -atsite SER_160 HIS_237 ASP_206 -rdmut 1


        $ python3 ASM.py -i 6eqe.pdb -arq gpu 0 -mode gpu-med -g MD -explicit water -atsite SER_160 HIS_237 ASP_206 -mut SER_160-MET


        $ python3 ASM.py -i 6eqe.pdb -arq gpu 0 -mode gpu-med -g MD -explicit water -atsite SER_160 HIS_237 ASP_206 -mut SER_160-MET HIS_237-ALA


        $ python3 ASM.py -i paracetamol.pdb -arq sander -mode med


        $ python3 ASM.py -i 6eqe.pdb -arq sander -mpi 4 -mode low -g MD -explicit water
