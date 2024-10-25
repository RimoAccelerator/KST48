# KST48: A Powerful Tool for MECP locating

KST48 is a purely-python program to locate Minimum Energy Crossing Points (MECPs), based on the searching algorithm by Bearpark, Robb and Schlegal in 1994 (Chemical physics letters, **1994**, 223(3): 269-274.).  An early implement of this algorithm is Harvey's pioneering work in 1998 (Theoretical Chemistry Accounts, **1998**, 99(2): 95-99.), which remains the most popular program to look for MECPs. Harvey's Fortran-based MECP program works in combination with Gaussian: by invoking Gaussian to calculate the forces and energies of the two states, a BFGS minimization algorithm was performed and gives an MECP starting from the input geometry.

KST48 provides a much more user-friendly and extensible solution for MECP locating, with a much more simplified preparation procedure. The optimization is based on the GDIIS algorithm, and is generally as effective. In addition to the locating of MECPs between different spin states by ground-state DFT methods, as implemented in Harvey's program, KST48 also supports the crossing between excited states with TD-DFT methods, or any other situations where the gradient can be read from a standard quantum chemical program.

In addition to the locating of MECP without any geometrical constrains, KST48 supports constrained optimization of bond lengths or angles, as well as the 1D- or 2D- scanning in the 3N-7 dimension of space with E1=E2 satisfied.

![logo](https://github.com/RimoAccelerator/KST48/blob/main/logo_KST48.png)

# Requirements
Python 3

Numpy
(KST48 is written with Numpy 1.19.3. However, it only ultilizes the most basic features of Python and Numpy, and should be compatible with any popular versions)

A quantum chemical program (currently supported are Gaussian, BAGEL and ORCA)

# Usage
1. Put kst48.py in any folder you want, together with a folder names JOBS. The task jobs will be placed here.
2. Prepare an input file:
```
#This subset is required. It controls your quantum chemistry tasks.
nprocs = 8 #processors
mem = 8GB # memory to be used. change this into the maxcore value if you want to use ORCA
method = wb97xd def2sv # your keywords line. It will be presented in the job files. Don't write guess=mix or stable=opt; they will be added automatically.
td1 =  # keywords for TD-DFT of state A (only for Gaussian; please write it in the tail part for ORCA)
td2 = # keywords for TD-DFT of state B (only for Gaussian; please write it in the tail part for ORCA)
mp2 = false #set true for MP2 or doubly hybrid calculation in Gaussian
charge = 0
mult1 = 1 # multiplicity of state A
mult2 = 3 # multiplicity of state B
mode = stable #normal; stable; read; inter_read; noread

#This subset is optional. It controls the convergence threshols, and the details of the GDIIS algorithm. Shown here are the default values.
dE_thresh = 0.000050
rms_thresh = 0.0025
max_dis_thresh = 0.004
max_g_thresh = 0.0007
rms_g_thresh = 0.0005
max_steps = 100
max_step_size = 0.1
reduced_factor = 0.5 # the gdiis stepsize will be reduced by this factor when rms_gradient is close to converge

# This subset controls which program you are using, and how to call them
program = gaussian  #gaussian or orca
gau_comm = 'g16'
orca_comm = '/opt/orca5/orca'
xtb_comm = 'xtb'

#Between *geom and *, write the cartesian coordinate of your initial geometry (in angstrom)
*geom
 C                 -0.03266394   -2.07149616    0.00000000
 H                  0.76176530   -1.26928217    0.00000000
 H                 -0.91673236   -1.36926867    0.00000000
 *
#If you have anything to be put at the end of the input file, write it here. This part is especially useful for ORCA: you can add anything related to your keywords here.
*tail1
#-C -H 0
#6-31G(d)
#****
*
*tail2
*

#This subset controls the constraints. R 1 2 1.0 means to fix distance between atom 1 and 2 (start from 1) to be 1.0 angstrom.
*constr
#R 1 2 1.0
#A 1 2 3 100.0 # to fix angle 1-2-3 to be 100 degree
#S R 1 2 1.0 10 0.1 # to run a scan of R(1,2) starting from 1.0 angstrom, with 10 steps of 0.1 angstrom
#S R 2 3 1.5 10 0.1 # you can at most set a 2D-scan
*
```
3. If the input file is called inp, then simply run `python3 kst48.py inp`. The output will be printed to the screen.
You will see the convergence of each step. After the MECP is converged, you will find its input and output files at current folder.

Since Oct. 2023, the geometry can also be defined by the following way, using an external link:
*geom
@a.log # Gaussian .log, .gjf, and .xyz files are supported.
*
It will read the last geometry that appears in a.log, and use to be the the initial geometry.
The atom fixation can be also defined, as shown in the input_example.txt file.

# Modes
There are several modes you can choose:
1. normal
  No special treatment.
2. stable
  Sometimes you might concern whether you are using stable wavefunctions in the MECP locating procedure. By setting mode=stable, KST48 will automatically run a wavefunction stability calculation of the initial geometry, and the so-generated wavefunction will be read for the following optimization.
  Note: please explicitly write UKS when you are using ORCA with the stable mode, to calculate a singlet state.
3. read
  If you have already had the wavefunction file (.chk for Gaussian or .gbw for ORCA), you could put them together with kst48.py, named as a.chk and b.chk (or .gbw for ORCA) respectively. By setting mode=read, KST48 will read them and check the stability, then do the optimization.
4. inter_read
  In some very hard cases, a stable wavefunction may not be correct. For example, for some molecule, you can find a stable RHF wavefunction for its singlet state, but there is a UHF wavefunction significantly lower in energy, and only can be obtained by SCF reading the triplet wavefunction (together with careful convergence controlling, such as guess=mix or even scf=vshift=-1000 in Gaussian). By setting mode=inter_read, state B will firstly be calculated in the first step, and then state A reads its wavefunction. Stability will be automatically checked. It is your work to write any convergence controlling keywords.
5. noread
  In normal mode, each step in the optimization reads the wavefunction from the previous step. In noread mode, no wavefunction is read.

# BAGEL Support
Since Feb. 2022, interface to BAGEL has been added to KST48. In order to invoke BAGEL, you need to set three options in the input file:
```
program = bagel  
bagel_comm = mpirun -np 36 /opt/bagel/bin/BAGEL
bagel_model = model.inp 
state1 = 0 #only set it for the multireference calculation using BAGEL
state2 = 1 #only set it for the multireference calculation using BAGEL
```
As you see, now a input template (here named model.inp) is required. Necessary keywords for a force calculation in BAGEL should be properly set in this file, and leave the geometry part as \*geometry\*. Then KST48 will automatically change the spin multiplicity and "target" option according to your mult1, mult2, state1 and state2. An example model.inp is attached along with the code (example_bagel_model_file.inp).
(The old version has a bug when reading CASPT2 energy. It has been fixed since 2023Oct.)

# Gaussian ONIOM Support
In ths 2023Oct version, the ONIOM calculation in Gaussian is supported. In order to use this feature, please set the following options:
```
isONIOM = true # if you are using the ONIOM feature of Gaussian, set this true.
chargeAndMultForONIOM1 = 0 1 0 1 0 1 0 1 0 1 0 1 #only useful for ONIOM calculation
chargeAndMultForONIOM2 = 0 1 0 1 0 1 0 3 0 3 0 3 #only useful for ONIOM calculation
```
Once isONIOM is set to be true, KST48 will not read the charge, mult1 and mult2 options, and turns to the chargeAndMultForONIOM1 and chargeAndMultForONIOM2 instead. In the geometry section, the cartesian coordinates with the layer info should be given in the following format:   C        2.02254521   -0.36259181   -0.01365118 M H 8. 
(Note that it is different from a standard ONIOM input generated by GaussView, which might be  C     **0**   2.02254521   -0.36259181   -0.01365118 M H 8 ).
**External geometry file is not supported for ONIOM.**
Also do not forget to contain the connectivity information in the tail1 and tail2 section.


# Structure Interpolation
Since Jun. 2022, an **experimental** feature has been added to do an interpolation between two structures. In some cases, this could aid you on finding a good guess for the MECP between two structures. By using this feature, set the "\*LST1" and "\*LST2" part in the input file. Put the geometries like this:
```
*lst1
 C                 -2.05168000   -3.89619700   -2.03479500
...
*
*lst2
 C                 -4.08538500    4.08668300   -0.76532000
...
*
```
Once "lst" part is detected, KST48 will ask you how many intermediate structures you are going to use for interpolation. The two geometries will be aligned, and the linear interpolation will be performed to generate several intermediate structures between lst1 and lst2. At present, only Gaussian is supported. Once the intermediate structures are approved, KST48 will invoke Gaussian to calculate the energies of both states for each structure, and give a list of energies as the output. You can use the one with minimal A-B energy difference as an initial guess for crossing point locating.

# Fix-DE Method
Since the **2024**Oct version (kst48_2024Oct.py), a new feature named "Fix-DE" has been added. It does not optimize an MECP; instead, it locates the lowest-energy structure for state A, while its energy difference between state B is constrained to be the value set by *fix_dE* option. This feature is achieve by optimizing the following function, which is a Lagrangian multiplier method:

$E_1 + \lambda (E_1 - E_2 - E_{target})$

It is especially useful for searching for a structure with an IP or EA exactly controlled, or the lowest-energy structure with a given excitation energy. An example is shown in tutorial.

# KST48Freq Module

KST48_Freq is a submodule for vibrational analysis of a crossing point, which is updated on Oct16, 2024.

The vibrational analysis is only reliable when the structure is a stationary point. However, for an MECP calculation, a converged structure is NOT a stationary point on the 3N-6 degrees of freedom. This problem can be solved by projecting the Hessian matrix into the subspace spanned by the 3N-7 degrees of freedom, in which a correctly-located MECP is indeed a stationary point.

The KST48Freq module is designed for this process. It requires two Gaussian fchk files for the frequency calculation of state A and B. During this process, the force constant matrix for each state is read, and an effective Hessian is obtained by:

$\boldsymbol{H}_{eff}=(1-\boldsymbol{P})((1-\lambda)\boldsymbol{H}_1+\lambda \boldsymbol{H}_2)(1-\boldsymbol{P})$

$\lambda=|\boldsymbol{g}_1|/(|\boldsymbol{g}_1-\boldsymbol{g}_2|)$

in which $\boldsymbol{P}$ is the projection matrix that projects out the 6 degrees of freedom for rotation and translation, and the 1 degree of freedom on the direction $\boldsymbol{g}_1-\boldsymbol{g}_2$

**Usage**:

Simply run `python3 kst48_freq.py`. It is an interactive program:

```Freq
Freq fchk file for state A?a.fchk
Freq fchk file for state B?b.fchk
Cosine for the angle formed by f1 and f2 is 0.9999592481491033
Frequencies for the HEff of CP:
[-2.13821824e-05 -2.13821824e-05 -4.41518653e-06  7.24372107e-06
  7.24372107e-06  1.51873184e-05  2.94043963e-02  2.68187257e+03
  3.04107649e+03]
Now it's time to output something.
Give me a freq log file so that I could replace the vibration information
and output a new file to kst48_freq.out.
Note that only the frequencies and modes are replaced,
and all the other things are remained.
It can be read by GoodVibes to obtain free energy correction.
Your freq .log file?5_A_freq.log
```

Then a Gaussian-type file is outputted to kst48_freq.out, which can be visualized by GaussView. Ideally it should contain 3N-7 frequencies and normal modes. In some rare cases, it may contain less, because normal modes with wavenumber less than 5 cm-1 is removed during outputting.

The output file can be recognized by GoodVibes for the calculation of thermal correction.

# Tutorial
Several tutorials are collected in [KST48 Tutorial.pdf](https://github.com/RimoAccelerator/KST48/blob/main/KST48%20Tutorial.pdf) and [Tutorial Fix_DE Method.pdf](https://github.com/RimoAccelerator/KST48/blob/main/Tutorial%20Fix_DE%20Method.pdf). Some typical applications are included:


+ Singlet-triplet MECP for transition metal complex
+ S1-T1 MECP of a TADF molecule
+ Gaussian ONIOM Calculation: S1/T1 MECP of Styrene in Explicit DCM
+ ORCA CASSCF: S0/S1 crossing point of CH2
+ ORCA SF-TDDFT: S-T crossing point for CH2
+ Fix DE Method


# Citation
Any use of this code MUST cite this Github Page (https://github.com/RimoAccelerator/KST48/), and it is encouraged to cite the author's first article using KST48. The citation is listed as the following:

1. Y. Ma, A. A. Hussein, ChemistrySelect 2022, 7, e202202354.

2. Yumiao Ma.  KST48: A Powerful Tool for MECP locating. https://github.com/RimoAccelerator/KST48, accessed on xxxx.xx.xx

The citation will be updated once KST48 is published in a journal or as preprint.

# Bug Report
You are welcomed to report any bug or problems to ymma@bsj-institute.top.

