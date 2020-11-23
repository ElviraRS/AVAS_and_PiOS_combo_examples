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



def readgeom(xyzfile):
 input=open(xyzfile,'r')
 lines=input.readlines()
 input.close()
 N_lines=len(lines)
 AtomsPyscf = []
 for n in range (2,N_lines):
  AtomsPyscf.append(lines[n])
 AtomsPyscf = "\n".join(AtomsPyscf)
 return AtomsPyscf

mol=gto.Mole()
xyz_input= "triad01_i.xyz"
mol.atom=readgeom(xyz_input)
mol.basis = '6-31g**'
mol.verbose=5
mol.spin =0
mol.build()


mf=scf.RHF(mol)
mf.chkfile = 'hf_triad01_i.chk'
mf.kernel()

N_Core,N_Act, N_Virt,nelec,coeff=MakeActiveSpace(mol,mf)
print " # of core orbs    ", N_Core
print " # of active orbs  ", N_Act
print " # of virtual orbs ", N_Virt
norb=N_Act
nalpha=nelec/2
nbeta=nelec/2


mycas = mcscf.CASSCF(mf, norb, [nalpha,nbeta])
AS=range(N_Core,N_Core+N_Act)
mycas =mycas.state_average_([0.1428,0.1428,0.1428,0.1428,0.1432,0.1428,0.1428])
mycas.chkfile ='cas_triad01_i.chk'
mycas.fcisolver.nroots = 7
mycas.fix_spin_(ss=0)
activeMO = mcscf.sort_mo(mycas,coeff,AS,base=0)
mycas.verbose = 5
mycas.max_cycle_macro = 150
mycas.kernel(activeMO)

mc = mcscf.CASCI(mf, norb, [nalpha,nbeta])
mc.__dict__.update(scf.chkfile.load('cas_triad01_i.chk', 'mcscf'))
mc.fcisolver.nroots = 7
mc.fix_spin_(ss=0)
mc.verbose = 6
mc.kernel()


def analyze_charges(Charges,n_excited_states):
 # here for charge analysis I did not includ H atoms, since the main redistribution of charges happens on pi-system C, N, O atoms

  Anthracene = [28,29,30,32, 34,36,38, 39,40,41, 43,45,47,49] 
  Phenol = [1,7,9,10,12,14,15] #here when O loses H, O involved in the pi-conjugated systems of the phenolic ring
  Pyridine = [3,16,17,19,20,22]

  chg_A=[]
  chg_F=[]
  chg_P=[]
  n_A=len(Anthracene)
  n_F=len(Phenol)
  n_P=len(Pyridine)
  for state in range(1,n_excited_states+1):
    dchg_A=0.0
    dchg_F=0.0
    dchg_P=0.0
    for k in range(n_A):
     ia=Anthracene[k]-1
     dchg_A+=(Charges[state][ia]-Charges[0][ia])
    chg_A.append(dchg_A)

    for k in range(n_F):
     ia=Phenol[k]-1
     dchg_F+=(Charges[state][ia]-Charges[0][ia])
    chg_F.append(dchg_F)

    for k in range(n_P):
     ia=Pyridine[k]-1
     dchg_P+=(Charges[state][ia]-Charges[0][ia])
    chg_P.append(dchg_P)

  for i in range(n_excited_states):
      print ('charge differences for Anthracene, Phenol and Pyridine units for excited state w.r.t. GS %d  %5.5f    %5.5f   %5.5f   '% (i+1,chg_A[i],chg_F[i],chg_P[i]))

def read_charges(casci_object,n_excited_states):
  chgs=[]
  for i in range(1+n_excited_states):
    ncore = mc.ncore
    ncas = mc.ncas
    nocc = ncore + ncas
    nelecas = mc.nelecas
    casdm1= mc.fcisolver.make_rdm1(mc.ci[i], ncas, nelecas)
    mo_coeff = mc.mo_coeff
    mocore = mo_coeff[:,:ncore]
    mocas = mo_coeff[:,ncore:nocc]
    dm1 =(numpy.dot(mocore, mocore.T) * 2
             + reduce(numpy.dot, (mocas, casdm1, mocas.T)))
    ovlp_ao = mc._scf.get_ovlp()
    pop,chg= mc._scf.mulliken_meta(mc.mol, dm1, s=ovlp_ao)
    chgs.append(chg)
  return chgs

charges=read_charges(mc,6)
analyze_charges(charges,6)


ci_nevpt_e1 = mrpt.NEVPT(mc, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mc, root=1).kernel()
ci_nevpt_e3 = mrpt.NEVPT(mc, root=2).kernel()
ci_nevpt_e4 = mrpt.NEVPT(mc, root=3).kernel()
ci_nevpt_e5 = mrpt.NEVPT(mc, root=4).kernel()
ci_nevpt_e6 = mrpt.NEVPT(mc, root=5).kernel()
ci_nevpt_e7 = mrpt.NEVPT(mc, root=6).kernel()



