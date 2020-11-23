This an example of applying PiOS for multiple pi-subsystems in a system of interest:
anthracene-phenol-pyridine molecular triads studied in :
E. R. Sayfutyarova and S. Hammes-Schiffe, "Substituent Effects on Photochemistry of Anthracene-Phenol-Pyridine Triads Revealed by Multireference Calculations,
 J. Am. Chem. Soc., 142, 487 (2020).

1) pay attention to the AssignTag function, which was used to adjust the number of electrons on proton donor and acceptor atoms, 
as we studied the electronic states along the proton transfer coordinates. We run calculations on a grid with 24 points. 
The input files were different only in terms of geometries.
2) you have to modify the ActiveSpaceMaker.py for this application: define separate lists of atoms for each of three pi-systems, 
and then we selected only two HOMO and two LUMO type orbitals for each of them
 
