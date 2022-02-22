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

These arrays are mainly used in the C parts of the code_saturne source, though
the interior and boundary *faces → vertices* connectivity is also
visible in the Fortran code. Remember that in Fortran code, arrays
are always one-based (i.e. the first element of an array has index 1),
while in C code, the natural indexing is zero-based, but one-based
indexing may also be used for arrays visible from Fortran code, or for arrays
using global numbers. In code_saturne, zero-based indexes are often used with
one-based data, for example when defining element connectivities,
where element ids are usually one-based (both as a convention
to simplify mapping to Fortran, and in the case of *cells → faces*
connectivities, so as to use the sign to determine face orientation).
For C code, when there are no mapping constraints due to Fortran,
the recommendations are the following:

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
element *i<sub>e</sub>*, are the values ranging from indexes
*i<sub>start</sub> = i<sub>e</sub>* included to
*i<sub>end</sub> = i<sub>e</sub>+1* excluded (past-the-end index).

The number of values associated with *i<sub>e</sub>* is determined by:
`array_index[i_e+1] - array_index[i_e]`, whether the index
is zero-based or one-based.

For an indexed array of *n* elements, the size the index array should thus
be equal to *n+1* (and not *n* as would be the case for regular 1-d or
strided arrays), and the total size of `array_data` is equal to
`array_index[n]` for a zero-based index, or
`array_index[n] - array_index[0]` in general.

similar popular data structures
-------------------------------

Readers familiar with *Compressed Sparse Row* or similar matrix or
graph representations may already have noted the similarity with
the indexed arrays described here. In the case of CSR matrix structures,
2 data arrays are often associated with 1 row index: one array definining
the column indices, and a second one defining the associated values.

This is in reality no different than using an indexed array as described here
to define a *faces → vertices* connectivity, and also associating
data (for example coordinates) to vertices.

In code_saturne, matrix non-diagonal terms usually correspond to cell faces,
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
uses of indexed arrays in code_saturne.

Recommended structure
---------------------

For mesh connectivities (adjacencies) that may be used across
several functions, the \ref cs_adjacency_t structure
allows storing an manipulation indexed arrays.


