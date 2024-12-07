The folder contains

1.	PDBs.zip folder contains all the PDB files used in the work. The "_A" suffix corresponds to active, while "_I" corresponds to inactive conformational state of the kinase.

2.	cmapCalcElecBlock.m computes contact map, pairwise electrostatic energy and block approximation using PDB file name, pH, distance and sequence cut-off. (Requires user input)
		1.	The user has to give input file name without file extension (if "test.pdb", then only "test" as input)
		2.	System pH. Can only run for 7, 5, 3.5 and 2 (default is 7)
		3.	srcutoff which is a distance cut-off in Angstrom (default is 5 Angstrom excluding H atoms)
		4.	BlockSize for block approximation
		5.	Secondary structure file generated from STRIDE webserver in the 'structfile' variable and must be present in the same folder. The naming of the structure file 				must be prefixed by "struct". For the test.pdb the structure file should be named as structtest.txt.In Linux systems, uncomment the lines 12-14 and comment and 			comment 15. The system must have stride downloaded in it with stride run file in the same folder.
	
	The output (o) will be 
		1.	residue-wise contact map ("contactmapmatElec" prefix), 
		2.	residue-wise electrostatic energy with distance between the residue pair ("contdistElec" prefix)
		3.	after block approximation, we will have a block contact map ("contactmapmatElecB" prefix), 
		4.	block electrostatic energy with distance between the block pair ("contdistElecB" prefix). 
		5.	"BlockDet" prefixed matrix mapping the residue number (column 1) with its respective block number (column 2), 
		6.	"BlockSize" prefixed variable containing block size,
		7.	"disr" prefixed variable contains the position of non-helical, non-strand residues (identified by STRIDE webserver) along with Glycine residues at any 					secondary structure region,
		8.	"ppos" prefixed variable containing proline residue position at any secondary structure region.

3.	FE_Coupling.m file for generating folding and thermodynamic coupling free energy (Requires user input)
		1.	The user must run cmapCalcElecBlock.m before FE_Coupling.m and its output must be in the same folder as FE_Coupling.m unless mentioned in the *.m file. 
		2.	The user has to give input file name without file extension (if "test.pdb", then only "test" as input)
		3.	Under ene variable, vdW interaction energy in kJ/mol should be provided
		4.	Under DS variable, entropic cost in kJ/(mol.K) should be provided
		5.	Under DCp variable, heat capacity change in kJ/(mol.K) per native contact should be provided
		6.	Under T variable, temperature in Kelvin should be provided
		7.	Under IS variable, ionic strength in Molar unit (M) should be provided
	o	The output will be MATLAB data-file (*mat extension) containing all the variables generated in the run.

4.	There are two excel files (*xlsx extension) containing
	o	data of 274 individual kinases in KinaseDatabase_All.xlsx and 104 pairs in KinaseDatabase_ActiveInactivePairs.xlsx file
	o	Each excel file contains two sheets suffixed "Database" and "Residues".
		1.	The sheet with "Database" suffix contains kinase family, name ("_2" suffixed for second-catalytic domain), HGNC name, Uniprot id, Domain boundary, PDB id, 				conformational state, vdW interaction energy, number of microstates, complementary PDB id and state (only in *Pair.xlsx), Intermediates with cut-off 1RT, 2RT, 				3RT, full name and sequence.
		2.	The sheet with "Residues" suffix contains kinase family, name ("_2" suffixed for second-catalytic domain), HGNC name, Uniprot id, PDB id, complementary PDB id 				(only in *Pair.xlsx), orthosteric and allosteric (MT3, AAS, PDIG, PIF, CMP, MPP and DRS) residues at structurally aligned identical positions.

5.	The two *plot.m and *data.mat files for plotting the data provided in research article
	
	o	Requirements for plotting the data present in a research article
		1.	MATLAB version 21 and after
		2.	Download AlaScanData.mat and FullVarData.mat in the same folder where DataPlot.m and AlaScanPlot.m are present
	
	o	Instructions for plotting the 1D- or 2D-free energy profiles/landscapes, residue probability along reaction coordinate, global residue folding probability and 				stability, positive/negative and effective thermodynamic coupling free energies.
		1.	Use DataPlot.m, which utilizes FullVarData.mat file to plot the files.
		2.	The input can be kinase name or PDB id. The name and id should be same as that provided in supplementary table S3 or in KinaseDatabase_ActiveInactivePairs.xlsx 			or in KinaseDatabase_All.xlsx file.
		3.	User can give both PDB id as well as name, both individual or together (name and pdb ids must match, else error will be shown).
		4.	In case the kinase's active-inactive pair is used. There will be prompt on the MATLAB command window, asking the conformational state of the kinase. Type 				A/Active or I/Inactive for active and inactive conformational state, respectively.
		5.	In case of ABL1 kinase, the prompt will show three states- A for active, I1 for inactive 1 and I2 for inactive. User has to type one of it to see the results.
		6.	The FullVarData.mat contains 379 rows with information on 379 kinase (170 + 104 pairs*2 + inactive 1 ABL1) in 14 columns as a cell vector. The data stored in 				columns are
			•	Column 1: Family
			•	Column 2: Kinase name ("_2" suffix is added for second catalytic kinase domain)
			•	Column 3: PDB id
			•	Column 4: Conformation state
			•	Column 5: vdW interaction energy (in kJ/mol)
			•	Column 6: Two column vector, where 1st column contains residue number and corresponding value in 2nd column represent the block number in which the 					respective residue lie.
			•	Column 7: One-dimension free-energy profile
			•	Column 8: Residue folding probability
			•	Column 9: Residue folding probability at each reaction coordinate (number of structured blocks)
			•	Column 10: Free-energy landscape 
			•	Column 11: Positive thermodynamic coupling free energy
			•	Column 12: Negative thermodynamic coupling free energy
			•	Column 13: Effective thermodynamic coupling free energy
			•	Column 14: Total number of sampled microstates
		7.	The DataPlot.m plots
			•	one-dimension free-energy profile with total number of sampled microstates mentioned in the plot (plot 1)
			•	residue folding probability across reacting coordinate (plot 2)
			•	residue folding stability (plot 3)
			•	free-energy landscape (plot 4)
			•	effective thermodynamic coupling free energy (plot 5)
			•	mean per residue effective thermodynamic coupling free energy (plot 6)
		8.	User can use the FullVarData.mat to plot positive/negative thermodynamic coupling free-energy and residue folding probability.
	
	o	Instructions for plotting the data of Alanine scan use AlaScanPlot.m
		1.	Use AlaScanPlot.m, which utilizes AlaScanData.mat file to plot the files.
		2.	Users can give kinase family, name and PDB id either, individually or together (family, name and pdb ids must match, else error will be shown).
		3.	The residue to be investigated should be provided by the user in "mutpos" variable. If the residue is either Alanine, Glycine or Proline, an error will be 				shown, as these residues aren’t mutated during Alanine scan analysis.
		4.	The AlaScanData.mat contains data variable with 9 rows containing information about 9 representative kinases (one from each family except RGC) in 18 columns as 			a cell vector. The data stored in columns are
			•	Column 1: Family
			•	Column 2: Kinase name
			•	Column 3: PDB id
			•	Column 4: Residue undergoing mutation
			•	Column 5: Residue number undergoing mutation
			•	Column 6: Vector with each column containing Cα-Cα distance from the mutation site in the row.
			•	Column 7: Wild-type one-dimension free-energy profile (1D array of size: total number of blocks i.e.NumBlock)
			•	Column 8: One-dimension free-energy profile of each mutated position (2D vector of size: NumBlock*mutated residue i.e., mut_res)
			•	Column 9: Wild-type positive thermodynamic coupling free energy (2D vector of size: Total residue i.e., Tot_res*Tot_res)
			•	Column 10: Positive thermodynamic coupling free energy of each mutated position (3D vector: Tot_res*Tot_res*mut_res)
			•	Column 11: absolute differential coupling indices (DCI) obtained using differential coupling matrix (2D vector: Tot_res*mut_res)
			•	Column 12: Parameters obtained after fitting (2D array of size: mut_res*3, where columns 1, 2 and 3 contain amplitude, coupling distance and shift 					obtained after fitting)
			•	Column 13: Confidence interval of fitting parameters (1D cell array of size: mut_res with each cell containing 3*2 matrix where rows 1, 2 and 3 contain 				lower (column 1) and upper (column 2) limit of amplitude, coupling distance and shift, respectively) 
			•	Column 14: Fitted y-values at a distance interval of 0.2 nm using fit parameters, given in column 13 (1D cell array of size: mut_res) 
			•	Column 15: Standard deviation of Cα-Cα distance from the mutation site at an interval of 0.2 nm
			•	Column 16: Mean of Cα-Cα distance from the mutation site at an interval of 0.2 nm
			•	Column 17: Mean of absolute DCI at an interval of 0.2 nm
			•	Column 18: Standard deviation of absolute DCI at an interval of 0.2 nm
		5.	The AlaScanPlot.m plots
			•	one-dimension free-energy profile of both wild-type and Alanine mutated residue (plot 1)
			•	absolute mutational response vs. mutated residue distance, with standard deviations as bars (black for distance and red for response), and exponential 					fitted line with its equation (plot 2)
		6.	Users can use the AlaScanData.mat to get confidence intervals for the fit parameters. 
