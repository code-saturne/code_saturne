<!--
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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

\page cs_dg_node_level_parallelism Node-level parallelism and accelerators

[TOC]

Host and device parallelism {#sec_prg_host_device}
===========================

Most current accelerated computing nodes are based on a host-device architecture,
where the **host** is based on a general-purpose processor, to which one or
more **devices** are adjoined to offload some operations.
On classical, non-accelerated nodes, only the host is present.

The most common devices are GPGPUs (General Purpose Graphical Processing Unit),
but other types of accelerators such as NEC's Vector Engine or FPGA
(Field Programmable Gate Array) could be considered.

Though accelerators can provide for greater performance and energy efficiency,
Leveraging their power cannot be done without adaptations in the programming
model.

Memory models
-------------

Accelerators usually have dedicated memory, with high bandwidth, which is separate
from the host memory. Copying between the host and device memory occurs latency,
and must be done as infrequently as possible. Mainstream accelerator programming
models may provide both separate host and device memory accesses, with
explicit exchanges, and "unified shared memory", allowing access of memory both from
the host and devices in an almost transparent manner, so as to provide
programmability and maintainability. The underlying implementation is usually
based on a [paging](https://en.wikipedia.org/wiki/Memory_paging) mechanism,
which may incur additional latency if not sufficiently well understood
and used (so the associated programming models may provide functions for
prefetching or providing hints to as how the memory is actually used).

Available memory on devices is often more limited than that on the host,
so allocating everything explicitly on the device could cause us to run out
of available memory. Using unified shared memory can avoid this, as memory paging
may provide the "illusion" of having more memory on the device, though this
can seriously degrade performance when this mechanism kicks in.

Ideally, we could use unified shared memory in all cases, and this might be done
in the future, it seems safer for the present to provide control to the developer
over which type of memory is used. So in code_saturne, when allocating memory
which might be needed on and accelerator, the \ref CS_MALLOC_HD function should
be used, specifying the allocation type with a \ref cs_alloc_mode_t argument.
In a manner similar to the older \ref BFT_MALLOC macro, this provides
some instrumentation, and is construed as a portability layer between several
programming models.

Programming models
------------------

### A few possible programming models

As classical C and current C++ or Fortran standards cannot express all the
possible parallelism, programming for accelerators may be done using either:

- Directive-based approaches, such as:
  * [OpenMP](https://en.wikipedia.org/wiki/OpenMP) 5.x,
  * [OpenACC](https://en.wikipedia.org/wiki/OpenACC),
  * Specific compiler `#pragma`s, such as for NEC's Vector Engine.

- Dedicated mainstream language extensions, such as:
  * [CUDA](https://docs.nvidia.com/cuda/cuda-c-programming-guide/#introduction),
  * [HIP](https://developer.amd.com/resources/rocm-learning-center/fundamentals-of-hip-programming/),
  * [DPC++](https://www.intel.com/content/www/us/en/develop/documentation/oneapi-programming-guide/top/oneapi-programming-model/data-parallel-c-dpc.html)

- Dedicated libraries, such as [Kokkos](https://kokkos.org).

- Languages designed specifically for HPC, such as [Chapel](https://chapel-lang.org).

- Domain-specific languages (DSL), which are usually based on a form of preprocessing
  to generate complex code in a mainstream language from simpler patterns tailored
  to a application domain.

Note that the mainstream language extensions listed above, as well as Kokkos,
and many DSLs are all based on C++.

None of these approches is currently as ubiquitous or portable as the C and
Fortran basis with host-based OpenMP directives on which most of code_saturne
is built:

- OpenMP would be expected to be the most portable solution here, but handling of
  accelerators was quite incomplete up until OpenMP 5.2 (and should be further
  improved in OpenMP 6.0), and most compilers installed on our machines do not
  yet support this version of the standard.

- OpenACC currently has more mature accelerator offload support than openMP, but
  seems supported by fewer vendors, so it could be considered mainly as a
  short-term approach.

- Kokkos is well established, but would add a critical and very intrusive
  dependency on all architectures, so is avoided for now.

### Selected programming models

Given these constraints, the current strategy regarding accelerator support is
the following:

- Use CUDA for the most common hot-spots (sparse linear-system solutions and
  gradient reconstruction), so as to benefit from previous work on the code,
  albeit in a non-portable manner.

- Use OpenMP offload when available in other parts of the code
  (requiring a recent compiler to actually use accelerated devices, but allowed
  (and at worst ignored) also on base code).

- Possibly migrate some code to C++, or DPC++ if this standard gains traction,
  to be able to benefit from convergence or at least similarity of C++-based
  approaches. This could be the preferred long-term solution.

- Use the aforementioned memory management wrappers from \ref cs_base_accel.cxx
  so as to be able to combine approaches and manage arrays used in accelerated
  algorithms from non-accelerated functions in other sections of the code.

Host-level parallelism
----------------------

Host-level parallelism is currently based on OpenMP constructs. Parallel loops
using threads are used where possible. Though this is not used yet, OpenMP
tasks could also be used to benefit from additional parallelization
opportunities.

Vectorization may also be used locally to enhance performance, whether handled
automatically by the compiler or handled explicitely through directives.
In practice, as code_saturns is mostly memory-bound, the benefits of
vectorization are limited, so improving the vectorization of various
algorithms is not a priority.

Device-level parallelism
------------------------

Various devices may be considerd, but the main targets are currently
GPGPUs.

As mentioned above, exploiting parallelism can be based on CUDA, DPC++,
or OpenMP directives.

Note that parallelism on GPGPU's is usually based on massive multi-threading,
where operations on an array may usually be divided into a series of
chunks (blocs), where each block is scheduled to run on available processors
(ideally in unspecified order), and computation of a given block is
itself multi-threaded.

Computational kernels launched on a device from the host are usually at least
in part asynchronous with the host (at least with CUDA and OneAPI/DPC++),
so that parallelism between the host and device may be exploited when the
algorithm allows this.
