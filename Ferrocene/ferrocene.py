import numpy
from pyscf import gto
from pyscf.gto import mole
from pyscf.gto import moleintor
from pyscf import scf
from pyscf import ao2mo
from pyscf import mcscf
from functools import reduce
from pyscf import fci
from pyscf.tools import fcidump
from ActiveSpaceMaker import MakeActiveSpace
from pyscf import mrpt

mol=gto.Mole()
mol.atom='''
Fe         0.00002        0.00035        0.00009
C         -0.23129        1.19815        1.71660
C         -1.21118        0.15030        1.71647
C         -0.51743       -1.10540        1.71638
C          0.89119       -0.83370        1.71638
C          1.06805        0.58996        1.71655
H         -0.43481        2.25429        1.70913
H         -2.27856        0.28292        1.70890
H         -0.97354       -2.07950        1.70826
H          1.67643       -1.56875        1.70790
H          2.00976        1.10960        1.70855
C         -1.06805        0.58996       -1.71666
C          0.23125        1.19809       -1.71675
C         -0.89115       -0.83367       -1.71646
H         -2.00970        1.10971       -1.70788
C          1.21111        0.15027       -1.71667
H          0.43470        2.25423       -1.70862
C          0.51744       -1.10539       -1.71649
H         -1.67648       -1.56861       -1.70723
H          2.27846        0.28305       -1.70847
H          0.97356       -2.07949       -1.70764
'''
mol.basis = 'cc-pvtz'
mol.verbose=4
mol.spin =0
mol.build()


mf=scf.ROHF(mol)
mf.chkfile = 'hf_FeCp2.chk'
mf.kernel()
#chkfile = 'hf_FeCp2.chk'     you can use this to restart from previous HF calculations if you saved the chk file
#mol = scf.chkfile.load_mol(chkfile)
#mf.__dict__.update(scf.chkfile.load(chkfile, 'scf'))


N_Core,N_Act, N_Virt,nelec,coeff=MakeActiveSpace(mol,mf)
print " # of core orbs    ", N_Core
print " # of active orbs  ", N_Act
print " # of virtual orbs ", N_Virt
norb=N_Act
nalpha=(nelec+mol.spin)/2
nbeta=(nelec-mol.spin)/2


#===================this section is optional==============
# if you would like to check the generated orbitals
HF_en =mf.mo_energy
HF_MOcoeff =mf.mo_coeff
HF_occ = mf.mo_occ
ova=mol.intor_symmetric("cint1e_ovlp_sph")
e_d=numpy.diag(HF_en)
cect=numpy.dot(HF_MOcoeff,numpy.dot(e_d,HF_MOcoeff.T))
f=numpy.dot(ova,numpy.dot(cect,ova))
en = numpy.diag(numpy.dot(coeff.T,numpy.dot(f,coeff)))
fname = 'generated_AS_beforeCAS'
from pyscf.tools import molden
with open(fname+'.molden','w') as thefile:
    molden.header(mol,thefile)
    molden.orbital_coeff(mol,thefile,coeff,ene=en, occ=mf.mo_occ)
#==========================================================


mycas = mcscf.CASSCF(mf, norb, [nalpha,nbeta])
AS=range(N_Core,N_Core+N_Act)
mycas.chkfile ='cas_FeCp2.chk'   # it is useful to keep the chk file for CASSCF in case you want to run some subsequent CASCI and NEVPT2 calculations
mycas.fcisolver.nroots = 1
mycas.fix_spin_(ss=0)  # in newer PYSCF it is better to specify ss which is usually your mol.spin
activeMO = mcscf.sort_mo(mycas,coeff,AS,base=0)
mycas.verbose = 5
mycas.max_cycle_macro = 150
mycas.kernel(activeMO)

mc = mcscf.CASCI(mf, norb, [nalpha,nbeta])
mc.__dict__.update(scf.chkfile.load('cas_FeCp2.chk', 'mcscf'))
mc.fcisolver.nroots = 1
mc.fix_spin_(ss=0)
mc.verbose = 6
mc.kernel()


ci_nevpt_e1 = mrpt.NEVPT(mc, root=0).kernel() 



