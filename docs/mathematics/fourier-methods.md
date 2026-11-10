# Fourier Methods

[Back to Index](../index.md)

The `Fourier` class in `Numerics.Mathematics` provides fast Fourier transforms (FFT) and the operations built on them: real-data transforms, linear convolution, correlation and the sample autocorrelation function.

## Conventions

All transforms store a complex sequence of length $n$ as a real array of length $2n$ with interleaved real and imaginary parts. The forward transform follows the Numerical Recipes sign convention:

```math
H_k = \sum_{j=0}^{n-1} h_j \, e^{+2\pi i jk/n}
```

This is the complex conjugate of the `numpy.fft.fft` convention for real input; for any input, the forward transform equals $n$ times `numpy.fft.ifft`. The inverse transform returns $n$ times the inverse DFT, so divide by $n$ to recover the original sequence.

## Methods

| Method | Lengths | Description |
|--------|---------|-------------|
| `FFT(data, inverse)` | power of 2 | In-place radix-2 complex FFT |
| `BluesteinFFT(data, inverse)` | any $n \ge 1$ | Complex DFT of any length; passes power-of-2 lengths to `FFT` |
| `RealFFT(data, inverse)` | power of 2 | FFT of real data, returning the packed positive-frequency half |
| `Convolve(data1, data2)` | any | Full linear convolution of two real sequences |
| `Correlation(data1, data2)` | power of 2 | Circular correlation of two real sequences (not normalized) |
| `Autocorrelation(series, lagMax)` | any | Sample autocorrelation function at lags $0..\text{lagMax}$ |

## Transforms of Any Length

`FFT` and `RealFFT` require a power-of-2 length. `BluesteinFFT` computes the same transform for any length using Bluestein's algorithm [[1]](#1): the identity $jk = \tfrac{1}{2}\left(j^2 + k^2 - (k-j)^2\right)$ turns the DFT into a convolution with the chirp $e^{\pm i\pi m^2/n}$, which is evaluated with power-of-2 FFTs of length at least $2n - 1$. The cost is $O(n \log n)$, about four to six times that of a power-of-2 FFT of similar length.

```cs
using Numerics.Mathematics;

// Three complex points 1, 2, 3 (imaginary parts zero)
var data = new double[] { 1, 0, 2, 0, 3, 0 };
Fourier.BluesteinFFT(data);
Console.WriteLine($"H0 = {data[0]:F4} + {data[1]:F4}i");   // H0 = 6.0000 + 0.0000i
Console.WriteLine($"H1 = {data[2]:F4} + {data[3]:F4}i");   // H1 = -1.5000 + -0.8660i

// Inverse transform, then divide by n to recover the input
Fourier.BluesteinFFT(data, inverse: true);
Console.WriteLine($"h1 = {data[2] / 3:F4}");                // h1 = 2.0000
```

## Linear Convolution

`Convolve` returns the full linear convolution $c_k = \sum_j a_j b_{k-j}$ of two real sequences, of length $n_1 + n_2 - 1$, the same as `numpy.convolve(a, b, mode="full")`. Both sequences are zero-padded to a power of 2 of at least $n_1 + n_2 - 1$, so there is no wraparound, and are transformed with `RealFFT`, which costs about half of a complex FFT. This is the operation needed to add independent discrete random variables: the probability mass function of a sum is the convolution of the mass functions.

```cs
// Sum of two fair dice: convolve the two mass functions
var die = new double[] { 1, 1, 1, 1, 1, 1 };
for (int i = 0; i < die.Length; i++) die[i] /= 6.0;
double[] sum = Fourier.Convolve(die, die);       // totals 2..12
Console.WriteLine($"P(total = 7) = {sum[5]:F6}");   // P(total = 7) = 0.166667
```

The round-off error is absolute, about machine epsilon times $\log_2 N$ times the largest output, rather than relative to each term. Entries far below the peak, such as deep tails of a convolved mass function, can carry round-off of order $10^{-16}$ times the peak, and exact zeros can come back as tiny negative values; clamp them when a non-negative result is required.

For a pair of real sequences of equal length $n$, `Convolve` was 1.8 to 4 times faster than the equivalent complex-FFT convolution (the same padding, two complex `FFT` calls and one inverse) for $n$ from 1,000 to 262,144 on .NET 10.

## Autocorrelation

`Autocorrelation` returns the standard sample autocorrelation

```math
r_k = \frac{\sum_{t=0}^{n-1-k} (x_t - \bar{x})(x_{t+k} - \bar{x})}{\sum_{t=0}^{n-1} (x_t - \bar{x})^2}
```

as a two-column array of lags and values. It is computed by FFT on a copy zero-padded to a power of 2 of at least $2n$, so the circular correlation equals the linear one at every lag. The default maximum lag is $\min(10 \log_{10} n, n - 1)$.

```cs
var series = new double[] { 1, 3, 2, 5, 4, 6, 5, 8 };
double[,] acf = Fourier.Autocorrelation(series, 3)!;
for (int k = 0; k <= 3; k++)
    Console.WriteLine($"r({acf[k, 0]}) = {acf[k, 1]:F4}");
// r(0) = 1.0000, r(1) = 0.2447, r(2) = 0.4120, r(3) = -0.1813
```

---

## References

<a id="1">[1]</a> Bluestein, L. I. (1970). A linear filtering approach to the computation of discrete Fourier transform. *IEEE Transactions on Audio and Electroacoustics*, 18(4), 451-455. https://doi.org/10.1109/TAU.1970.1162132

---

[Back to Index](../index.md)
