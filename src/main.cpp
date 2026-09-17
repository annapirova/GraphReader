#include "matrix.h"
#include "mmio.h"


int main(int argc, const char* argv[]) {
	spMtx<double> matrix(argv[1]);
	
	std::cout << "Matrix is downloaded:\n";
	std::cout << "rows: " << matrix.m << " columns: " << matrix.n << "\n";
	std::cout << "non-zeroes: " << matrix.nz << "\n";	

    return 0;
}
