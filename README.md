# extended_guinier

This script implements a version of the "Extended Guinier" method. Briefly, the method addresses the problem that for IDPs the region for which Guinier analysis is valid is usually a very small range of q in which the data is usually very noisy. An ansatz is developed for the next term of the expansion of I(q) such that this information can also be used to extract the radius of gyration (with the caveat that it is limited to disordered proteins in present form).


## Citation

[W. Zheng, R. B. Best, "An extended Guinier analysis for intrinsically disordered proteins", JMB v430 pp2540-2553 (2018).](https://doi.org/10.1016/j.jmb.2018.03.007)


## Usage:

```
chmod a+x extended_guinier.py
./extended_guinier.py Nres datafile outfile
where:
  Nres - number of residues in protein
  datafile - experimental data
  outfile - script output

  datafile is expected to have the first two columns be q, I(q) and comment lines preceded by "#"
  outfile has q**2, log[I(q)], log[I(q)]_fitted

```

