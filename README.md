# Krylov
## by Valerio Fassi
### Krylov is a header only library inspired from Eigen.

## How to use

Clone the repository:
```
git git@github.com:bolt3x/Krylov.git
```

Inside the **src** folder there is the **Krylov.hpp** header,
to use it in your code compile with
```
-I./path/to/krylov/header
```

##Test Cases
You can test the library by using the make command and running that **test** file 
that will be located in the test directory.
Various test cases can be generated:

```
make test=cg
```
This will be generate a test case that uses the CG solver
```
make test=bcgstab
```
This will generate a test case that uses the BiCGStab solver
```
make test=cgs
```
This will generate a test case that uses the CGS solver

##The Sparse Approximate Inverse Preconditioner 

After running the **test** file an images folder will be created.
Here there will be 3 png file representing:

**original.png** the pattern of the original matrix
**dynamic.png** the pattern of the sparse approximate inverse with dynamic pattern 
**static.png** the pattern of the sparse approximate inverse with static pattern
#Example
using https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lanpro/nos4.html
**ORIGINAL MATRIX**
![original](https://user-images.githubusercontent.com/103378889/216479768-a5505586-5ec3-4e38-9fbf-e4c1109d1f6c.png)
**SPAI WITH DYNAMIC PATTERN**
![dynamic](https://user-images.githubusercontent.com/103378889/216479780-e0870163-d931-4700-a348-004d7f8ad3fc.png)
**SPAI WITH STATIC PATTERN**
![static](https://user-images.githubusercontent.com/103378889/216479790-1f47a430-a595-457f-9ad3-4d486583750f.png)

## Caveat
The sparse approximate inverse preconditioner with **DYNAMIC** pattern is not optimized
so at the moment the use of the **STATIC** version is suggested.
