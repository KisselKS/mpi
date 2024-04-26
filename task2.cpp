#include <iostream>
#include <vector>
#include <chrono>
#include <mpi.h>

using namespace std;

double sum_calc(const vector<double>& arr, int a, int b) {
    double sum = 0.0;
    for (int i = a; i < b; i++)
        sum = sum + arr[i];
    return sum;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    vector<int> sizes = { 10000000 };
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int i = 0; i < sizes.size(); i++) {
        int n = sizes[i];
        int p_size = n / size; 
        vector<double> arr(n);

        if (rank == 0) { 
            srand(time(nullptr));
            for (int i = 0; i < n; i++) {
                arr[i] = rand();
            }
        }

        MPI_Bcast(arr.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD); 

        int a = rank * p_size; 
        int b;
        if (rank >= size - 1)
            b = n;
        else a + p_size;

        auto time1 = chrono::steady_clock::now();
      
        double p_sum = sum_calc(arr, a, b);

        double result;

        MPI_Reduce(&p_sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

        auto time2 = chrono::steady_clock::now();
        chrono::duration<double> all_time = time2 - time1;

        if (rank == 0) {
            cout << "Size: " << n << "\nTime: " << all_time.count() << "\n";
        }
    }

    MPI_Finalize();

    return 0;
}