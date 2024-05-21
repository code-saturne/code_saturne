#!/usr/bin/env python3

#------------------------------------------------------------------------------
# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 2009-2024 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.
#-------------------------------------------------------------------------------

import math
import matplotlib.pyplot as plt
import numpy as np
import os.path
import struct

from argparse import ArgumentParser
from scipy.sparse import coo_matrix

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def init_read_binary_file(file_path, num_bytes):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    try:
        f = open(file_path, 'rb')
    #
    except FileNotFoundError:
        print(f'File "{file_path}" does not exist')
    #
    except PermissionError:
        print(f'File "{file_path}" Permission denied')
    #
    except IOError:
        print(f'File "{file_path}" I/O error detected')
    else:
        content = f.read(num_bytes)

    return content

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def continue_read_binary_file(file_path, start_byte, num_bytes):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    # Sanity checks on the file have been done with the initial call

    with open(file_path, 'rb') as f:
        f.seek(start_byte)
        content = f.read(num_bytes)

    return content

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_binary_matrix(file_path):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    data = init_read_binary_file(file_path, 3)

    # Decode the metadata

    intsize, doublesize, binmode = struct.unpack('B B c', data)
    intsize = int(intsize)
    doublesize = int(doublesize)
    binmode = binmode.decode('utf-8')

    # Build integer decode format

    if binmode == 'l':
        binfmt = '<'
    else:
        binfmt = ''

    if intsize == 4:
        intfmt = 'L'
    elif intsize == 8:
        intfmt = 'Q'
    else:
        sys.exit(f"Size of int:{intsize}. One expects 4 or 8.")

    byte_shift = 3
    data = continue_read_binary_file(file_path, byte_shift, intsize)
    n_rows = struct.unpack(f'{binfmt}{intfmt}', data)[0]

    byte_shift += intsize
    data = continue_read_binary_file(file_path, byte_shift, intsize)
    nnz = struct.unpack(f'{binfmt}{intfmt}', data)[0]

    # Output the main information

    print(f'int size: {intsize}')
    print(f'double size: {doublesize}')
    print(f'binary mode: {binmode}')
    print(f"int decode format: '{binfmt}{intfmt}'")
    print(f'n_rows: {n_rows}')
    print(f'nnz: {nnz}')

    byte_shift += intsize
    data = continue_read_binary_file(file_path, byte_shift, intsize*nnz)
    row_nums = struct.unpack(f'{binfmt}{nnz}{intfmt}', data)

    byte_shift += nnz*intsize
    data = continue_read_binary_file(file_path, byte_shift, intsize*nnz)
    col_nums = struct.unpack(f'{binfmt}{nnz}{intfmt}', data)

    byte_shift += nnz*intsize
    data = continue_read_binary_file(file_path, byte_shift, doublesize*nnz)
    values = struct.unpack(f'{nnz}d', data)

    return n_rows, nnz, row_nums, col_nums, values

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def read_binary_vector(file_path):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    data = init_read_binary_file(file_path, 3)

    # Decode the metadata

    intsize, doublesize, binmode = struct.unpack('B B c', data)
    intsize = int(intsize)
    doublesize = int(doublesize)
    binmode = binmode.decode('utf-8')

    # Build integer decode format

    if binmode == 'l':
        binfmt = '<'
    else:
        binfmt = ''

    if intsize == 4:
        intfmt = 'L'
    elif intsize == 8:
        intfmt = 'Q'
    else:
        sys.exit(f"Size of int:{intsize}. One expects 4 or 8.")

    byte_shift = 3
    data = continue_read_binary_file(file_path, byte_shift, intsize)
    n_rows = struct.unpack(f'{binfmt}{intfmt}', data)[0]

    # Output the main information

    print(f'int size: {intsize}')
    print(f'double size: {doublesize}')
    print(f'binary mode: {binmode}')
    print(f"int decode format: '{binfmt}{intfmt}'")
    print(f'n_rows: {n_rows}')

    byte_shift += intsize
    data = continue_read_binary_file(file_path, byte_shift, doublesize*n_rows)
    values = struct.unpack(f'{n_rows}d', data)

    return n_rows, values

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def post_matrix(matrix):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    plt.figure(figsize=(8,6))

    if args.scatter:
        scatter = plt.scatter(matrix.col, -matrix.row, c=matrix.val,
                              cmap='plasma', marker='o', s=1)
        plt.colorbar(scatter)
    else:
        plt.spy(matrix, markersize=0.5)

    plt.show()

    # eigenvalues = np.linalg.eigvals(matrix.toarray())
    # print(eigenvalues)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def post_vector(n_rows, values, tag="vector"):
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    plt.figure(figsize=(8,6))
    plt.grid(True, linestyle=':', lw=0.8)
    plt.title(f"{tag} values")
    plt.ylabel("Occurence")
    plt.xlabel("Values")

    plt.hist(values, bins=min(200, max(10, math.floor(0.001*n_rows))))
    plt.show()

    plt.figure(figsize=(8,6))
    plt.grid(True, linestyle=':', lw=0.8)
    plt.title(f"{tag} values")
    plt.ylabel("Values")
    plt.xlabel("Id")

    plt.plot(values, marker='o', ms=2, mfc='w', linestyle='None')
    plt.show()

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Main part
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if __name__ == '__main__':

    parser = ArgumentParser(description=""" Analyze binary files
    produced by code_saturne using functions from cs_matrix_util.c""")

    parser.add_argument("-m", "--mat", metavar="MAT_FILE_NAME")
    parser.add_argument("-r", "--rhs", metavar="RHS_FILE_NAME")
    parser.add_argument("-s", "--sol", metavar="SOL_FILE_NAME")
    parser.add_argument("--scatter", action='store_true',
                        help="""Use the scatter function rather than the spy
                        function for plotting a matrix""")
    parser.add_argument("--post", action='store_true',
                        help="Activate predefined post-processings")

    args = parser.parse_args()

    if args.mat:
        if os.path.isfile(args.mat):

            n_rows, nnz, rows, cols, values = read_binary_matrix(args.mat)

            print(f"Min values: {min(values)}")
            print(f"Max values: {max(values)}")

            # Create a sparse matrix in COO format from read data

            matrix = coo_matrix((values, (rows, cols)))

            # Post-processing

            if args.post:
                post_matrix(matrix)

        else:
            print(f"File '{args.mat}' is not an existing file.")

    if args.rhs:
        if os.path.isfile(args.rhs):

            n_rows, rhs_values = read_binary_vector(args.rhs)

            print(f"Min. RHS values: {min(rhs_values)}")
            print(f"Max. RHS values: {max(rhs_values)}")
            print(f"Mean RHS values: {np.mean(rhs_values)}")
            print(f"Var. RHS values: {np.var(rhs_values)}")

            # Post-processing

            if args.post:
                post_vector(n_rows, rhs_values, tag='RHS')

        else:
            print(f"File '{args.rhs}' is not an existing file.")

    if args.sol:
        if os.path.isfile(args.sol):

            n_rows, sol_values = read_binary_vector(args.sol)

            print(f"Min. solution values: {min(sol_values)}")
            print(f"Max. solution values: {max(sol_values)}")
            print(f"Mean solution values: {np.mean(sol_values)}")
            print(f"Var. solution values: {np.var(sol_values)}")

            # Post-processing

            if args.post:
                post_vector(n_rows, sol_values, tag='Solution')

        else:
            print(f"File '{args.sol}' is not an existing file.")
