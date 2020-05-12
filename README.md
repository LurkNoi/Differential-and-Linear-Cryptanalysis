# Differential-and-Linear-Cryptanalysis

## Issues

- [x] Still need to manually simplify POS when using MILP to model S-Box
    - possible solution: solve with pyeda.boolalg.espresso

- [x] Equations is too long for large (8-bit) S-Box
    - possible solution: ignore small bias

## TODO

- [x] simplify logical condition with python

- [ ] modeling Feistel Network


## Reference

- Howard M. Heys. 2002. A tutorial on linear and differential cryptanalysis. DOI: https://doi.org/10.1080/0161-110291890885

- Sun S. et al. Automatic Security Evaluation of Block Ciphers with S-bP Structures Against Related-Key Differential Attacks. DOI: https://doi.org/10.1007/978-3-319-12087-4_3

- Abdelkhalek, A. et al. MILP Modeling for (Large) S-boxes to Optimize Probability of Differential Characteristics. DOI: https://doi.org/10.13154/tosc.v2017.i4.99-129

- Baignères T. et al. How Far Can We Go Beyond Linear Cryptanalysis?. DOI: https://doi.org/10.1007/978-3-540-30539-2_31
