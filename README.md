# MCPU
//start	

MCPU - Monte-Carlo protein simulation program


											Directories
											
mcpu_prep - the directory containing the code to create input files
sim - the directory with files prepared for simulations of any specific protein. This is 
also where output is stored
src_mpi_umbrella - the directory with source code
src_mpi_umbrella/cfg - the configuration file
config_files - the directory with parameters





Creates a grid of simulations with multiple temperatures. 
At each temperature, multiple cores can be run. If desired, umbrella sampling can be 
implemented such that a harmonic term with respect to native contacts is added to the 
energy of the form:

U_umbrella = 1/2*K_BIAS*(N - S)^2 (Eq. 1)

Where N is the number of native contacts for a proposed configuration, S is the set point,
and K_BIAS is the spring constant. Umbrella biasing can also be turned off, in which case
each core at a given temperature has the same conditions. There is also the option of
implementing replica exchange, where a core can exchange with its neighbors in the grid 
 along either the temperature and set point directions.




NOTE: At the moment, this code works for PDB files with up to 8000 atoms and 
1000 residues. To increase this, change the values for MAX_ATPMS and MAXSEQUENCE, 
respectively, in define.h and recompile


1. Create necessary input files: 
	<PDB_ID>.triple
	<PDB_ID>.sctorsion
	<PDB_ID>.sec_str
To create the first two files, run save_triple.c (in the mcpu_prep directory): 
	./save_triple <PDB_ID>
with triple.energy, sct.energy, and <PDB_ID>.fasta in the directory. This may take a few minutes to run.
NOTE: It is important that the FASTA file have 80 characters per line

Create <PDB_ID>.sec_str manually. File contains secondary structure assignment for each protein residue (see publication [1]).
first line: use input secondary structure? (9/0 = yes/no)
second line: secondary structure type (H/E/C = helix/sheet/coil)
For most applications, the first line is entirely 0's (no input secondary structure) and the second line is set to entirely C's 

Place input files, along with the pdb file, in the directory sim/DHFR/files/
(currently contains sample input files for DHFR)


2. Edit configuration options in cfg file. The most relevant options (without changing the potential) are:

									NATIVE PROTEIN DATA
	NATIVE_DIRECTORY -- Contains a set of PDB files that will be used to initialize simulation for each respective core. The PDB files in this directory should be named 0.pdb, 1.pdb, 2.pdb, etc., and the ith core will initialized with PDB file i.pdb. If you wish to initialize all cores with the same input file, set this option to None, and simply edit NATIVE_FILE and STRUCTURE_FILE below 
	NATIVE_FILE and STRUCTURE_FILE -- input PDB file for simulations (folded structure, single chain, no hydrogens). Even if the NATIVE_DIRECTORY option above is set to something other than None, these options should still be specified to allow for RMSD computation
	PDB_OUT_FILE -- path to simulation output. Should be formatted as {path to output directory}/{protein name}, such that all simulation output will be saved in {path to output directory} and the output files will all incorporate {protein name} in their name.
									
									MONTE-CARLO PARAMETERS

	MC_STEPS -- length of the simulation
	MC_PDB_PRINT_STEPS -- frequency of outputting coordinates to a pdb file
	MC_PRINT_STEPS -- frequency of outputting energies to log file

									Replica Exchange Parameter
	MC_REPLICA_STEPS -- frequency of replica exchange. To turn off exchange, set to a value greater than MC_STEPS.


									SIMULATION PARAMETERS
	MC_TEMP_MIN -- Lowest simulation temperature in grid
	TEMP_STEP -- Spacing between successive temperatures in grid
	NODES_PER_TEMP -- How many cores are assigned to each temperature
	USE_CLUSTER -- Gives the frequency at which knowledge-based moves are attempted, given that the function LoopBackboneMove has been called to make a move (see [4] for description of knowledge-based moves). The function LoopBackboneMove is only called with probability 0.33, so the overall chance of attempting a knowledge move is 0.33*USE_CLUSTER. At MC steps above MAX_CLUSTERSTEP (see below), this value is automatically set to 0 and knowledge-based moves are no longer used.
		NOTE: Knowledge-based moves violate detailed balance and thus, steps that incorporate them should not be used to compute thermodynamic properties. But these moves may nonetheless be useful to incorporate at the beginning of a simulation to help the simulation find energy minima at intermediate numbers of native contacts (as in [5])
		NOTE: At the moment, if the secondary structure file has any helices (H characters), then the value in the cfg file will be ignored and USE_CLUSTER will be set to 0.5 by default. This can be changed if desired by modifying the function LoopBackboneMove (in Move.h) and recompiling.
	MAX_CLUSTERSTEP -- Largest step at which knowledge-based moves (see [4]) are to be used. For all MC steps beyond this, these moves will be turned off. As above, this requires the secondary structure file to have no H characters.
	



									Umbrella parameters

	UMBRELLA--indicates whether or not umbrella sampling is to be used. 1 if so, 0 if not. All subsequent parameters in this section are moot if set to 0
	K_BIAS -- Spring constant for umbrella biasing. See equation (Eq. 1) above
	NUMBER_OF_CONTACTS_MAX	-- Highest set point to be used in umbrella biasing
	CONTACTS_STEP -- Separation between set points. These parameters set up a simulation grid with nodes_per_temp cores at each temperature, whose set points range from NUMBER_OF_CONTACTS_MAX and descending in increments of CONTACTS_STEP
	MIN_SEQ_SEP -- Minimum separation in sequence between residues for pairs of residues in the native PDB file to define a contact. Setting to values above 4 ensures that short range contacts in alpha helices are not included in contacts definition.
	CONTACT_CALPHA_CUTOFF --Maximum distance, in angstroms, for two alpha carbons to be considered in contact
	
					




	 								PARAMETER FILES
	-- direct these to the correct input file in the sim folder. 
		- TRIPLET_ENERGY_FILE is <PDB_ID>.triple (see step 1)
		- SIDECHAIN_TORSION_FILE is <PDB_ID>.sctorsion
		- SECONDARY_STRUCTURE_FILE is <PDB_ID>.sec_str
	
			



3. Additional parameters related to the energy function can be modified in define.h if desired, but this will require the program to be recompiled.
Contains weights for different energy terms (see publications [1], [3]): 
POTNTL_WEIGHT -- contact potential
HBOND_WEIGHT -- hydrogen bonding
TOR_WEIGHT -- torsional energy for amino acid triplets
SCT_WEIGHT -- side chain torsional energy
ARO_WEIGHT -- relative orientations of aromatic residues


4. Compile (if necessary) and run
If it is necessary to re-compile the code, one can do so from src_mpi_umbrella directory by simply typing ./compile (which runs the commands  gcc -c rng.c and mpicc -O3 -o fold_potential_mpi backbone.c -lm rng.o)

To run:
mpiexec -n <# of procs> ./fold_potential_mpi <cfg filename>	




5. Data analysis
A separate readme for data analysis method used in [5] and [6] is forthcoming


Publications:
[1] J.S. Yang et al., Structure 15, 53 (2007)
[2] J. Tian et al., PLOS Comp. Bio., in press
[3] J. Xu, L. Huang, E. I. Shakhnovich, Proteins 79, 1704 (2011)
[4] W. Chen, J.S. Yang, E. I. Shakhnovich, Proteins 66, 682 (2007)
[5] A. Bitran, W. M. Jacobs, X. Zhai, E. I Shakhnovich, PNAS (2020)
[6]A. Bitran, W. M. Jacobs, E.I. Shakhnovich, in press


//end