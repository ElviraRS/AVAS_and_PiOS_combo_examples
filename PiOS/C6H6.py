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
from pyscf import mrpt
from pyscf.mcscf.PiOS import MakePiOS

mol=gto.Mole()
mol.atom='''
 C      0.000000    0.001262    0.000000
 C      0.000000    2.782738    0.000000
 C      1.204419    0.696631    0.000000
 C     -1.204419    0.696631    0.000000
 C      1.204419    2.087369    0.000000
 C     -1.204419    2.087369    0.000000
 H     -2.141507    2.628401    0.000000
 H      2.141507    2.628401    0.000000
 H     -2.141507    0.155599    0.000000
 H      2.141507    0.155599    0.000000
 H      0.000000   -1.080796    0.000000
 H      0.000000    3.864796    0.000000
'''
mol.basis ='aug-cc-pvtz'
mol.verbose=5
mol.spin =0
mol.build()


mf=scf.RHF(mol)
mf.chkfile = 'hf_c6h6.chk'
mf.kernel()

PiAtoms = [1,2,3,4,5,6]
N_Core,N_Act, N_Virt,nelec,coeff=MakePiOS(mol,mf,PiAtoms)
norb=N_Act
nalpha=nelec/2
nbeta=nelec/2

mycas = mcscf.CASSCF(mf, norb, [nalpha,nbeta])
AS=range(N_Core,N_Core+N_Act)
mycas =mycas.state_average_([0.1432,0.1428,0.1428,0.1428,0.1428,0.1428,0.1428])
mycas.chkfile ='cas_c6h6.chk'
mycas.fcisolver.nroots = 7
mycas.fix_spin_(ss=0)
activeMO = mcscf.sort_mo(mycas,coeff,AS,base=0)
mycas.verbose = 5
mycas.max_cycle_macro = 150
mycas.kernel(activeMO)

mycas = mcscf.CASCI(mf, norb, [nalpha,nbeta])
mycas.__dict__.update(scf.chkfile.load('cas_c6h6.chk', 'mcscf'))
mycas.fcisolver.nroots = 7
mycas.fix_spin_(ss=0)
mycas.verbose = 5
mycas.kernel()

ci_nevpt_e1 = mrpt.NEVPT(mycas, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mycas, root=1).kernel()
ci_nevpt_e3 = mrpt.NEVPT(mycas, root=2).kernel()
ci_nevpt_e4 = mrpt.NEVPT(mycas, root=3).kernel()
ci_nevpt_e5 = mrpt.NEVPT(mycas, root=4).kernel()
ci_nevpt_e6 = mrpt.NEVPT(mycas, root=5).kernel()
ci_nevpt_e6 = mrpt.NEVPT(mycas, root=6).kernel()

