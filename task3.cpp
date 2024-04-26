#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

double f(double x, double y) {
    return 0.5*sin(x) + 2*cos(y);
}

void dr(double* A, double* B, int rows, int n, double dx) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < n; j++) {
            if (j > 0 && j < n - 1) {
                B[i * n + j] = (A[i * n + j + 1] - A[i * n + j - 1]) / (2 * dx);
            }
            else {
                B[i * n + j] = 0.0;
            }
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int n = 1000;
    double dx = 0.01;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rows = n / size; 

    double* A1 = new double[rows * n];
    double* B1 = new double[rows * n];

    double* A = new double[n * n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = f(i * dx, j * dx);
        }
    }
 
    double* B = nullptr;
    if (rank == 0) {
        B = new double[n * n];
    }

    double time1 = MPI_Wtime();

    if (rank == 0) {
        MPI_Scatter(A, rows * n, MPI_DOUBLE, A1, rows * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
        delete[] A;
    }
    else {
        MPI_Scatter(nullptr, rows * n, MPI_DOUBLE, A1, rows * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 
    }

    dr(A1, B1, rows, n, dx); 

    MPI_Gather(B1, rows * n, MPI_DOUBLE, B, rows * n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

    double time2 = MPI_Wtime();
    double all_time = time2 - time1;

    if (rank == 0) {
        cout << "Size: " << n << "\nTime: " << all_time << "\n";
    }

    delete[] A1;
    delete[] B1;
    if (rank == 0) {
        delete[] B;
    }

    MPI_Finalize();
    return 0;
}
