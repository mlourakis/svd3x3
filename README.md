# Singular Value Decomposition for 3x3 matrices
The [singular value decomposition](https://en.wikipedia.org/wiki/Singular_value_decomposition) (SVD) of 3x3 matrices arises often in applications in computer vision, robotics, graphics, etc.

This is a C++ implementation of SVD customized for 3x3 matrices.
The SVD factors a 3x3 matrix $`M`$ as $`M = U S V^t`$ where $`U`$ and $`V`$ are orthogonal 3x3 matrices and $`S`$ is diagonal. 

The decomposition is based on Nick Higham's [polar decomposition](https://github.com/martinbis11/polar-decomposition-3x3/tree/master) and David Eberly's [symmetric eigen solver](https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf) for 3x3 matrices.

The code is templated so that the decomposition can be computed either in single or double precision. To the best of my knowledge, this is the only fast 3x3 SVD in C++ with this feature; approaches such as that of [McAdams et al.](https://minds.wisconsin.edu/handle/1793/60736) provide single precision floating point only.

The implementation is a header-only library, providing its own matrix operations and having no external dependencies. Matrices are expected to be in the C++ native row-major representation.<br>
To use it in any given project, copy all .h files from the ``include`` directory to a directory where your compiler can find them and ``#include "svd3x3.h"`` in your code. Check ``svd_demo.cpp`` for more details.

## Build
Create a ``build`` directory in the root of the cloned repository and run ``cmake``:

``mkdir build; cd build; cmake ..``

or, for a *release* build,

``cmake .. -DCMAKE_BUILD_TYPE=Release``

Finally build everything:

``make``

To run the demo, once in the ``build`` directory,

``./svd_demo``

## Performance
A benchmark program is included that allows you to measure the performance on your own machine. Run it from ``build`` with

``./svd_bench``

On my 3.2 GHz Intel Core i9-12900K, the benchmark showed that a 3×3 double precision SVD took 0.3727 microseconds (μs) on average.
This is only 1.5× slower than Eric Jang's [3×3 float SVD](https://github.com/ericjang/svd3), which implements the method of [McAdams et al.](https://minds.wisconsin.edu/handle/1793/60736) Float 3×3 SVD itself took 0.2966 μs, roughly 1.24× slower than Eric’s implementation.

## Math background
if the [polar decomposition](https://en.wikipedia.org/wiki/Polar_decomposition) of $`M`$ is $`M=Q H`$ and $`H=V S V^t`$, then an SVD is $`(Q V) S V^t \equiv U S V^t`$

## Cite as
If you use this code in your published work, please cite this repo:<br><br>
<pre>
@misc{lourakis2024svd3x3,
  title={Singular Value Decomposition for 3x3 Matrices},
  author={Manolis Lourakis},
  howpublished={[web page] \url{http://github.com/mlourakis/svd3x3}},
  year={2024},
  note={[Accessed on 4 Dec. 2024]},
}
</pre>
