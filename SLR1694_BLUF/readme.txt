This an example of applying PiOS for multiple pi-subsystems in the active site model of the SLR1694 BLUF photoreceptor studied in :
E. R. Sayfutyarova, J. J. Goings and S. Hammes-Schiffer, "Electron-Coupled Double Proton Transfer in the Slr1694 BLUF Photoreceptor:
 A Multireference Electronic Structure Study, J. Phys. Chem. B, 123, 439 (2019).

With a simple script you can generate all other input files (they will differ simply by xyz- and chk-files names in the input files).
Other geometries differ by positions of two hydrogens being moved between oxygen and nitrogen atoms.
1) pay attention to the AssignTag function, which was used to adjust the number of electrons on proton donor and acceptor atoms,
as we studied the electronic states along the proton transfer coordinates. We run calculations on a 2D grid along O-O and N-N.
The input files were different only in terms of geometries.
2) you have to modify the ActiveSpaceMaker.py for this application: define separate lists of atoms for each of four pi-systems,
and the transferring hydrogen atoms for the AVAS.

Note: In this work I did not include the O into the conjugated system of the Tyr, but it would be better if I did as its lone pair is clearly delocalized into Tyr ring. 
It became clear after later calculations with electrostatic embedding, affecting this even more.

