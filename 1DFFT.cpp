#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

// Define sample size (must be a power of two)
constexpr int SAMPLE_SIZE = 1024;

// Function declarations
void bit_reverse_permutation(vector<complex<double>>& x);
void compute_twiddles(vector<complex<double>>& W, int N);
void fft(vector<complex<double>>& X);

// Main function
int main() {
vector<complex<double>> X(SAMPLE_SIZE);
double f = 100.0;   // Frequency of the cosine wave (100 Hz)
double fs = 1000.0; // Sampling frequency (1000 Hz)

// Generate cosine wave samples
for (int n = 0; n < SAMPLE_SIZE; ++n) {
X[n] = cos(2 * M_PI * f * n / fs);
}

// Perform FFT
fft(X);

// Output the FFT results with formatting
cout << fixed << setprecision(4);
for (size_t i = 0; i < X.size(); ++i) {
cout << "X[" << i << "] = " << X[i] << endl;
}

return 0;
}

// FFT implementation
void fft(vector<complex<double>>& X) {
int N = X.size();
bit_reverse_permutation(X);

vector<complex<double>> W(N / 2);
compute_twiddles(W, N);

int stages = static_cast<int>(log2(N));
for (int s = 1; s <= stages; ++s) {
int m = 1 << s;
int half_m = m >> 1;

for (int k = 0; k < N; k += m) {
for (int j = 0; j < half_m; ++j) {
complex<double> t = W[(N * j) / m] * X[k + j + half_m];
complex<double> u = X[k + j];
X[k + j] = u + t;
X[k + j + half_m] = u - t;
}
}
}
}

// Bit-reversal permutation
void bit_reverse_permutation(vector<complex<double>>& x) {
int N = x.size();
unsigned int j = 0;
for (unsigned int i = 0; i < N; ++i) {
if (i < j)
swap(x[i], x[j]);

unsigned int m = N >> 1;
while (m >= 1 && j >= m) {
j -= m;
m >>= 1;
}
j += m;
}
}

// Compute twiddle factors
void compute_twiddles(vector<complex<double>>& W, int N) {
for (int k = 0; k < N / 2; ++k) {
double theta = -2.0 * M_PI * k / N;
W[k] = polar(1.0, theta);
}
}
