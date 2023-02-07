# Some Notes #

-Good C++ programmin. COde reasonably commented.
-Good the idea of using namespaces
-Diagonal preconditioner is just a scaling of the matrix, normally deas very little. One could try with a symmetric
Gauss siedel.
- Code compiles and runs without issues. More daring tests could have been tried!, with larger matrices.
- Some parmeters for the tests (like the method used, preconditioner and the matrix file name) could have been read from a file. a bit more work
but at least we would have easier way of running the tests.



# Minor stuff #
- For a vector you could have specialised a matrix with just one column.
- In `qr_solver.hpp` instead of `Matrix getQ() { return Q; }` prefer `const Matrix &getQ() const { return Q; }`, in this way you avoid the copy and guarantee that the instance does not change. You can still obtain a copy if you call `Matrix my_copy_of_Q = QR_instance.getQ();`
- In `InversePower` prefer `std::functions` to function pointers, the result is the same but code is more clear.
- Could have used `const` and `auto` for some local variables to make the code safer and more general.
