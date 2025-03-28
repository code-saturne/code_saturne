<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
-->

\page cs_dg_common_constructs Common construct types

[TOC]

In this section, commonly-used construct types whose use may require specific
explainations or recommendations are described.

Indexed arrays {#sec_prg_indexed_arrays}
==============

In many instance data such as mesh connectivity requires managing a variable
number of entries per data element. This is for example the case of
*faces → vertices* connectivity. The average number of vertices
per face is usually quite low, but the maximum number may be significantly
higher, so using an array with regular stride would be very inefficient
for some data sets.

A common solution to this problem is to use indexed arrays, in which an array
containing data is supplemented by a second array containing the start indexes
of entries in the data array for each element.

In C and C++ code, the natural indexing is zero-based, but one-based
indexing may also be used for some arrays as a vestige of prior Fortran code,
or for arrays using global numbers.

The recommendations are the following:

- Local index arrays should be zero-based.
- Global index arrays should be one-based. This should only concern
  indexes read from or written to file.
- When containing cell, face, or vertex connectivity information, data
  arrays may be either zero or one-based: zero based arrays are less
  error-prone so they should be preferred, but where element ids may be
  signed (so as to convey orientation information), one-based arrays are
  necessary. In a given structure, consistency is recommended, so if
  *cells → faces* connectivity requires one-based face numbers,
  an associated *faces → vertices* connectivity may also use
  one-based vertex numbers, even though vertices have no orientation.

Let us consider an array `array_data` indexed by a zero-based
`array_index` array. The values of `array_data` associated with
element <em>i<sub>e</sub></em>, are the values ranging from indexes
<em>i<sub>start</sub> = i<sub>e</sub></em> included to
<em>i<sub>end</sub> = i<sub>e</sub>+1</em> excluded (past-the-end index).

The number of values associated with <em>i<sub>e</sub></em> is determined by:
`array_index[i_e+1] - array_index[i_e]`, whether the index
is zero-based or one-based.

For an indexed array of _n_ elements, the size the index array should thus
be equal to _n+1_ (and not _n_ as would be the case for regular 1-d or
strided arrays), and the total size of `array_data` is equal to
`array_index[n]` for a zero-based index, or
`array_index[n] - array_index[0]` in general.

Similar popular data structures
-------------------------------

Readers familiar with *Compressed Sparse Row* or similar matrix or
graph representations may already have noted the similarity with
the indexed arrays described here. In the case of CSR matrix structures,
2 data arrays are often associated with 1 row index: one array definining
the column indices, and a second one defining the associated values.

This is in reality no different than using an indexed array as described here
to define a *faces → vertices* connectivity, and also associating
data (for example coordinates) to vertices.

In code\_saturne, matrix non-diagonal terms usually correspond to cell faces,
and the CSR matrix representation is very similar to that of a
*cells → faces* connectivity, except for the fact that a
standard CSR representation uses only "unsigned" column ids, whereas
face numbers may be signed in the matching mesh representation so as
to convey orientation (an alternative solution would be to use
a separate array for orientation, in which case the similarity to CSR
would be complete).

Indexed Array Example
---------------------

We illustrate the use of an indexed array to define a *faces → vertices*
connectivity for a simple surface mesh:

\anchor fig_example_indexed_array
![Example indexed array](prog_indexed_array.svg)

In this example we use 0-based (*0* to *n-1*) numbering, as in most
uses of indexed arrays in code\_saturne.

Recommended structure
---------------------

For mesh connectivities (adjacencies) that may be used across
several functions, the \ref cs_adjacency_t structure
allows storing an manipulation indexed arrays.

Parallel dispatch {#sec_prg_parallel_dispatch}
=================

In C++ code using the `.cpp` extension, the `cs_dispatch_context` class may
be used to run code defined as C++ lambda functions in parallel. The
corresponding code can be written in a portable manner, but executed using
different runtimes, such as OpenMP or CUDA.

Examples are provided in the source tree, in
- `tests/cs_dispatch_test.cpp`
- `src/alge/cs_benchmark.cpp`

Note that this class should not be called from code with the `.cxx` extension,
as when using CUDA, that code will ve compiled using a different compiler,
and using the class in both cases may confuse those compiler's avoidance of
duplicating C++ template code instanciation.

Base mechanism
--------------

The parallel dispatch is based on using templated `parallel_for` functions
including a loop (on a CPU) or kernel call (on a GPU), in which
a functor (usually generated automatically as a _lambda function_)
is called.

We explain the working of the parallel dispatch, on a few examples,
initiallt eliding the top level template aspects for clarity.

Simple parallel dispatch example
--------------------------------

Using the code from `cs_dispatch.h`, the following loop:

```{.cpp}
#pragma omp parallel for
for (size_t idx = 0; idx < n; idx++) {
  if (is_disabled[idx])
    continue;

  y[i] += a * x[i];
}
```

can be replaced by:

```{.cpp}
cs_dispatch_context ctx;
ctx.parallel_for(n, [=] (cs_lnum_t idx) {
  if (is_disabled[idx])
    return;

  y[i] += a * x[i];
});
```
Note here that the `continue` statement (skip to next loop element)
is replaced by `return` (exit called function), as the dispatch
mechanisms calls an intermediate function.

### Functors and lambda function
The previous code block is equivalent to:

```{.cpp}
ctx.parallel_for(n, a, is_disabled, x, y);
```

With the following implementation of the parallel_for method (on a CPU with OpenMP parallelism)

```{.cpp}
void
parallel_for(size_t         n,
             double         a,
             const bool    *is_disabled,
             const double  *x,
             double        *y)
{
  auto_generated func(is_disabled, a, x, y);

  #pragma omp parallel for
  for (size_t idx = 0; idx < n; idx++) {
    func(idx);
  }
}
```

and the following functor definition:

```{.cpp}
class auto_generated {

private:

  double            a_;
  const bool       *is_disabled_;
  const cs_real_t  *x_;
  cs_real_t        *y_;

public:

  // Constructor
  auto_generated(double         a,
                 const bool    *is_disabled,
                 const double  *x,
                 double        *y) :
    a_(a), is_disabled_(is_disabled), x_(x), y_(y) {}

  // Operator
  void operator() (size_t  idx)
  {
    if (is_disabled[idx])
      return;

    y[i] += a * x[i];
  }

};
```

A functor is simply a class implementing the `()` operator.

In the `ctx.parallel_for(n, [=] (cs_lnum_t idx) {...});`
syntax above, the `[]` syntax is a capture clause indicating the following code defines a lambda function. Such functions can access variables outside their enclosing scope, as specified in the capture clause. `[&]` means all variables are captured by reference, while `[=]` means all variables are captured by copy. The capture clause can be more complex, specifying different capture clauses for specific variables, but for the purposes of the code\_saturne parallel dispatch mechanism, we always use the `[=]` (capture all variables by copy) syntax, as it is the only one adapted to running functions (kernels) on devices using a different memory space.

Using this lambda capture mechanism, a functor similar to the one described above is generated automatically based on the enclosed code.

### GPU code generation

The practical interest here is that using functors allows separating the loop parallelism logic from that of the called functor. Combined with template mechanisms, the same code can be used to generate parallel code using different back-ends or scheduling mechanisms.
For example, to run on a CUDA-enabled device, the equivalent to the function above can be defined as:

```{.cpp}
cs_dispatch_context ctx;
ctx.parallel_for(n, [=] __host__ __device__ (cs_lnum_t idx) {
  if (is_disabled[idx])
    return;

  y[i] = a * x[i];
});
```

This function is almost identical to the original one, except for the `__host__ __device__` specifiers. To ensure the syntax matches exactly, we use a CS_F_HOST_DEVICE macro which expands to the above specifiers when compiled for CUDA, and is empty for simple CPU code.

So the loop is finally written as :

```{.cpp}
cs_dispatch_context ctx;
ctx.parallel_for(n, [=] CS_F_HOST_DEVICE (cs_lnum_t idx) {
  if (is_disabled[idx])
    return;

  y[i] = a * x[i];
});
```

### CUDA implementation

On a GPU, using CUDA, the parallel_for method is implemented as follows:

```{.cpp}
void
parallel_for(size_t         n,
             double         a,
             const bool    *is_disabled,
             const double  *x,
             double        *y){
  auto_generated func(i_face_cells, a, x, y, i_sum_type);

  cs_cuda_kernel_parallel_for<<<l_grid_size, block_size_, 0, stream_>>>
    (m->n_i_faces, f);
}
```

With the following CUDA kernel:

```{.cpp}
__global__ void
cs_cuda_kernel_parallel_for(cs_lnum_t      n,
                            auto_genrated  f) {
  // grid_size-stride loop
  for (cs_lnum_t id = blockIdx.x * blockDim.x + threadIdx.x; id < n;
       id += blockDim.x * gridDim.x) {
    f(id, args...);
  }
}
```

For readers not familiar with CUDA, note that  `blockId`, `blockDim`, and `threadIdx` are built-in kernel variables describing where a given thread executes.

### Asynchronous execution

The dispatch mechanism is asynchronous, so when the call to `parallel_for` returns, we are only assured that the computation is scheduled. The associated wait ensures we wait for the call to actually complete. On a CPU with OpenMP, the end of an OpenMP section involves an implicit barrier, so  wait is a no-op (does nothing). When  using a CUDA back-end, where kernel launches are actually implicit, it  contains a call to `cudaStreamSynchronize`,

Face-based loops with parallel assembly (scatter sum)
-----------------------------------------------------

For many operations such as balance computations, code\_saturne uses*scatter* based algorithms, where values computed at faces are contributed (_assembled_) to cell-based arrays.
When run on a relatively small number of threads (such as on a CPU), a specific face numbering is used to avoid race conditions (https://doi.org/10.1002/cpe.2852), leading to code similar to the following example.

```{.cpp}
for (int ig = 0; ig < n_groups; ig++) {
  #pragma omp parallel for
  for (int it = 0; it < m->numbering->n_threads; it++) {
    for (cs_lnum_t face_id = m->numbering->index[ig][it][0];
         face_id < m->numbering->index[ig][it][1];
         face_id++) {
      cs_lnum_t i = i_face_cells[face_id][0];
      cs_lnum_t j = i_face_cells[face_id][1];

      y[i] += a[k][0]*x[j];
      y[j] += a[k][1]*x[i];
    }
  }
}
```

When the number of threads becomes too massive, an alternative option is to simply use atomic sums, but can incur a significant performance penalty:

```{.cpp}
#pragma omp parallel for
for (cs_lnum_t face_id = 0; face_id < m->n_i_faces; face_id++) {
  cs_lnum_t i = i_face_cells[face_id][0];
  cs_lnum_t j = i_face_cells[face_id][1];

  #pragma omp atomic
  y[i] += a[k][0]*x[j];
  #pragma omp atomic
  y[j] += a[k][1]*x[i];
}
```

Using the dispatch mechanism, the external loops can be hidden in the implementation, leading to the following code:

```{.cpp}
cs_dispatch_context ctx;
cs_dispatch_sum_type i_sum_type = ctx.get_parallel_for_i_faces_sum_type(m);

ctx.parallel_for_i_faces(m, [=] CS_F_HOST_DEVICE (cs_lnum_t face_id) {
  cs_lnum_t i = i_face_cells[face_id][0];
  cs_lnum_t j = i_face_cells[face_id][1];

  cs_dispatch_sum(&y[i], a[k][0]*x[j], i_sum_type);
  cs_dispatch_sum(&y[j], a[k][1]*x[i], i_sum_type);
});
```

### Implementation
The lambda capture mechanism is the same as described in the previous examples (each different sets of captured arguments leading to the generation of a different functor).
On a CPU with OpenMP-based parallelism, the implementation expands to:

```{.cpp}
void
parallel_for_i_faces(cs_mesh_t            *m,
                     const cs_lnum_2_t    *i_face_cells,
                     const cs_real_2_t    *a,
                     const cs_real_t      *x,
                     cs_real_t            *y,
                     cs_dispatch_sum_type  i_sum_type)
{
  auto_generated func(i_face_cells, a, x, y, i_sum_type);

  for (int ig = 0; ig < n_groups; ig++) {
    #pragma parallel for
    for (int it = 0; it < m->numbering->n_threads; it++) {
      for (cs_lnum_t face_id = m->numbering->index[ig][it][0];
           face_id < m->numbering->index[ig][it][1];
           face_id++) {
        func(face_id);
      }
    }
  }
}
```

On a GPU using CUDA, the function body is replaced by:

```{.cpp}
{
  auto_generated func(i_face_cells, a, x, y, i_sum_type);

  cs_cuda_kernel_parallel_for<<<l_grid_size, block_size_, 0, stream_>>>
    (m->n_i_faces, f);
}
```

Using the same `cs_cuda_kernel_parallel_for` kernel template as for the simple example.

### Actual mechanism

For the dispatch mechanism to actually be usable with many different loops, it must be based on templates. For clarity, in the examples above, template parameters were expanded to actual types, so as to focus on the dispatch and loop mechanisms.

The `cs_dispatch` class also uses a CRTP (curiously recurring template pattern) to optionally allow forcing the compilation of a given section of code only on the CPU (using `cs_host_context`) or GPU (using `cs_device_context`). This can be useful when a multiple variants of a given algorithm are required for performance reasons, and generating only the variant which is expected to perform best on a given architecture is desired.
