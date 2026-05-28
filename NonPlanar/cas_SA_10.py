import numpy
from pyscf import gto
from pyscf import scf
from pyscf import mcscf
from NonPlanarPiOS import MakePiOS
from pyscf import fci
from pyscf import mrpt


# this is how we initialize a molecule: info about the total spin, charge, choice of a basis set and molecular geometry
mol=gto.Mole()
mol.atom="TH2S_10.xyz"
mol.basis = '631G*'
mol.verbose=4  # this is parameter for setting-up the level of detalization of the input file
mol.charge=0
mol.spin =0
mol.build()

#===================================Running RHF================================================

mf = scf.RHF(mol)
mf.chkfile = 'hf_10.chk'
mf.conv_tol = 1.e-11
mf.kernel()


#chkfile = 'hf_10.chk'
#mol = scf.chkfile.load_mol(chkfile)
#mf.__dict__.update(scf.chkfile.load(chkfile, 'scf'))

#========================================Active space orbital construction======================

PiAtoms = range (1,11)
N_Core,N_Act, N_Virt,nelec,coeff=MakePiOS(mol,mf,PiAtoms)

print(" # of core orbs    ", N_Core)
print(" # of active orbs    ", N_Act)
print(" # of virtual orbs    ", N_Virt)

#============Optional: dump the molden file to vizualize the constructed orbitals=======

HF_en =mf.mo_energy
HF_MOcoeff =mf.mo_coeff
HF_occ = mf.mo_occ
ova=mol.intor_symmetric("cint1e_ovlp_sph")

e_d=numpy.diag(mf.mo_energy)
cect=numpy.dot(mf.mo_coeff,numpy.dot(e_d,mf.mo_coeff.T))
f=numpy.dot(ova,numpy.dot(cect,ova))
en = numpy.diag(numpy.dot(coeff.T,numpy.dot(f,coeff)))

fname = 'constructed_orbitals_10'
from pyscf.tools import molden
with open(fname+'.molden','w') as thefile:
    molden.header(mol,thefile)
    molden.orbital_coeff(mol,thefile,coeff,ene=en, occ=mf.mo_occ)

#==========================CASSCF with two FCI solvers for two spin values===================

weights = numpy.ones(9)/9
solver1 = fci.addons.fix_spin(fci.direct_spin1.FCI(mol), ss=2) # <S^2> = 2 for Triplet
solver1.spin = 2
solver1.nroots = 4
solver2 = fci.addons.fix_spin(fci.direct_spin1.FCI(mol), ss=0) # <S^2> = 0 for Singlet
solver2.spin = 0
solver2.nroots = 5


mycas = mcscf.CASSCF(mf, N_Act, nelec)
AS=range(N_Core,N_Core+N_Act)
activeMO = mcscf.sort_mo(mycas,coeff,AS,base=0)
mycas.chkfile ='cas_10_SA.chk'
mycas.verbose = 5
mcscf.state_average_mix_(mycas, [solver1, solver2], weights)
mycas.max_memory=120000
mycas.max_cycle_macro = 50
mycas.kernel(activeMO)

#=====================Running CASCI and NEVPT2 for singlets using saved MOs from CASSCF=======

mc1 = mcscf.CASCI(mf, N_Act, nelec)
mc1.__dict__.update(scf.chkfile.load('cas_10_SA.chk', 'mcscf'))
mc1.fcisolver.nroots = 5
mc1.fcisolver.spin =0
mc1.fix_spin_(ss=0)
mc1.verbose = 6
mc1.kernel()


ci_nevpt_e1 = mrpt.NEVPT(mc1, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mc1, root=1).kernel()
ci_nevpt_e3 = mrpt.NEVPT(mc1, root=2).kernel()
ci_nevpt_e4 = mrpt.NEVPT(mc1, root=3).kernel()
ci_nevpt_e5 = mrpt.NEVPT(mc1, root=4).kernel()

#=====================Running CASCI and NEVPT2 for triplets using saved MOs from CASSCF=======
na=nelec//2+1
nb=nelec//2-1
nelec=(na,nb)
mol.spin=2

mc2 = mcscf.CASCI(mol, N_Act, nelec)
mc2.__dict__.update(scf.chkfile.load('cas_10_SA.chk', 'mcscf'))
mc2.fcisolver.nroots = 4
mc2.fcisolver.spin=2
mc2.fix_spin_(ss=2)
mc2.verbose = 6
mc2.kernel()


ci_nevpt_e1 = mrpt.NEVPT(mc2, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mc2, root=1).kernel()
ci_nevpt_e3 = mrpt.NEVPT(mc2, root=2).kernel()
ci_nevpt_e4 = mrpt.NEVPT(mc2, root=3).kernel()
