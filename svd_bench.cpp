// benchmark program for svd3x3
// uses templates so that SVD can be benchmarked in either single or double precision

#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <cmath>

#include "svd3x3.h"

#define N_SAMPLES 500000
#define N_REPEATS 100

// call SVD function, making sure that the compiler does not optimize it out
template <typename TReal>
void bench(TReal mat[9])
{
    svd3x3::decomp<TReal> svd3;
    TReal U[9], V[9], s[3];
    svd3(mat, U, s, V);

    volatile TReal sink = U[0] + U[4] + V[0] + V[8] + s[0] + s[1]; // prevent optimization
    (void)sink;
}

// random sample generation
template <typename TReal>
void setup(std::vector<std::array<TReal, 9>>& samples)
{
    std::mt19937 gen(42);
    std::uniform_real_distribution<TReal> dist(-100.0, 100.0);
    for (int i=0;i<N_SAMPLES;++i)
        for(int j=0;j<9;++j)
            samples[i][j] = dist(gen);
}

// compute mean and stddev
template <typename TReal>
TReal mean(const std::vector<TReal>& data)
{
    TReal s=0; for(auto x:data) s+=x;
    return s/data.size();
}

template <typename TReal>
TReal stdev(const std::vector<TReal>& data, TReal m)
{
    TReal s=0; for(auto x:data) s+=(x-m)*(x-m);
    return std::sqrt(s/data.size());
}

int main()
{
    using TReal = double; // or float

    std::vector<std::array<TReal, 9>> samples(N_SAMPLES);
    setup(samples);

    std::vector<double> timings;
    timings.reserve(N_REPEATS);

    for(int rep=0;rep<N_REPEATS;++rep){
        auto start = std::chrono::high_resolution_clock::now();

        for(int i=0;i<N_SAMPLES;++i){
            bench(samples[i].data());
        }

        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count();
        double per_op_us = elapsed_ns * 1e-3 / N_SAMPLES; // μs per SVD
        timings.push_back(per_op_us);
    }

    double m = mean(timings);
    double s = stdev(timings,m);
    std::cout << "Mean time per SVD: " << m << " μs\n";
    std::cout << "Stddev per SVD: " << s << " μs\n";
}
