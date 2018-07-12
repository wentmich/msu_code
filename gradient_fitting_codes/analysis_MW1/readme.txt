
$ root -l FieldSim_derivatives_V03.C
- Simulate multipole fields vs z and theta using COSY for magnets.
- Decompose harmonics at all z
- Evaluate at misalignment values dx,dy to estimate derivative of dipole
  components. These derivatives can be used to estimate misalignment of 
  magnets
  (notes.xlsx -- Notes on output from code.)
- Evaluate on-axis gradients to file. For example
   b_n02_m0.out, a_n02_m0.out are normal and skew quad components
   

$ runcosy91.sh COSY_field_Gen.fox
- COSY script used to evaluate multipole field
  It is called from .C scripts.
cosy91.fox -- copy of reference COSY beam physics code
lnCOSY91.bin -- link used to run COSY from another directory

$ filon_MP.c
- Subroutines used in harmonic analysis and other numerical techniques.

$ clean.sh (script to delete generated output files)