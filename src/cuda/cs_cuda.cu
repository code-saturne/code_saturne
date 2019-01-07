/*============================================================================
 * CUDA offloading support
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) IBM Corp. 2017, 2018

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
*/

/*----------------------------------------------------------------------------
 * Standard C/C++ library headers
 *----------------------------------------------------------------------------*/
#include <assert.h>
#include <stdio.h>

#include <algorithm>
#include <list>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/
#include "cs_cuda.h"
#include "cs_matrix.h"

#include "bft_mem.h"

/*----------------------------------------------------------------------------
 *  Have Eclipse CDT parser deal with CUDA specific syntax and data types.
 *----------------------------------------------------------------------------*/
#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__

namespace {
dim3 threadIdx;
dim3 blockIdx;
dim3 blockDim;
dim3 gridDim;
int warpSize;

template <typename T>
__device__ T __shfl_down_sync(unsigned mask, T var, unsigned int delta,
                              int width = warpSize);

void __syncthreads(void);

double atomicAdd(double *address, double val);
unsigned long long int atomicCAS(unsigned long long int *address,
                                 unsigned long long int compare,
                                 unsigned long long int val);
} // namespace

#define LAUNCH_KERNEL(...) any_kernels<__VA_ARGS__>()
#else
#define LAUNCH_KERNEL(...)                                                     \
  case getKindsID<__VA_ARGS__>():                                              \
    any_kernels<__VA_ARGS__><<<GridDim, BlockDim>>>();                         \
    break
#endif

/*----------------------------------------------------------------------------
 *  Implementation if CUDA is defined
 *----------------------------------------------------------------------------*/
#ifdef HAVE_CUDA_OFFLOAD
#include <cuda.h>
#include <cuda_runtime.h>

#define CHECK(x)                                                               \
  if (cudaError_t err = (x)) {                                                 \
    printf("CUDA error in '%s', at line '%d': %s\n", __FILE__, __LINE__,       \
           cudaGetErrorString(err));                                           \
    assert(false);                                                             \
  }

namespace {
// The number of devices available in a node. This is NOT the devices visible
// for the current rank.
int NumberNodeDevices = -1;
// The device ID associated with this rank
int DeviceID = -1;

// Total number of ranks in current node.
int NumberRanks = 0;

// Device properties.
cudaDeviceProp DeviceProperties;

// Function that flushes the argument to the GPU and launch the relevant
// kernels.
void flush(void);

// Class that manages the memory mapped to the GPU.
class MemoryManager {
public:
  // The types of the map.
  enum MapType { maptype_alloc, maptype_to, maptype_from, maptype_release };

  const char *getMapTypeName(MapType Type) const {
    switch (Type) {
    case maptype_alloc:
      return "Alloc";
    case maptype_to:
      return "To";
    case maptype_from:
      return "From";
    case maptype_release:
      return "Release";
    }
    return nullptr;
  }

private:
  // Information for 1 map
  class SlotInfo {
    // The pointer to the slot.
    char *Begin;
    // The pointer to the end of the slot (non-inclusive).
    char *End;
    // The pointer in the device.
    char *DevBegin;
    // The number of times the slot was referenced.
    size_t RefCount;

  public:
    char *getBegin(void) const { return Begin; }
    char *getEnd(void) const { return End; }
    size_t getSize(void) const { return End - Begin; }
    char *getDevBegin(void) const { return DevBegin; }
    size_t getCount(void) const { return RefCount; }
    void incCount(void) { ++RefCount; }
    void decCount(void) {
      assert(RefCount > 0 && "This slot ref count is already zero!");
      --RefCount;
    }

    SlotInfo() : Begin(nullptr), End(nullptr), DevBegin(nullptr), RefCount(1) {}
    SlotInfo(char *Begin, char *End, char *DevBegin)
        : Begin(Begin), End(End), DevBegin(DevBegin), RefCount(1) {}
    SlotInfo(char *Ptr, char *DevPtr, size_t Size)
        : Begin(Ptr), End(Ptr + Size), DevBegin(DevPtr), RefCount(1) {}
  };

  // List with all current maps.
  std::list<SlotInfo> SlotsInfo;

  void dumpSlots(void) {
#ifdef CS_CUDA_DEBUG
    printf("Mapped slots:\n");
    for (const auto &S : SlotsInfo)
      printf("-> (refcount=%d) %p-%p[:\n", (int)S.getCount(), S.getBegin(),
             S.getEnd());
#endif
  }

  // Variables to communicate reduction results.
  cs_real_t *ReductionHost = nullptr;
  cs_real_t *ReductionDev = nullptr;

  // Variable to control dual buffering of reduction results.
  cs_lnum_t ReductionResultVersion = 0;

  // Memory chunk for the current rank to use.
  void *GPUStorage = nullptr;

  //
  // Lists that keep the used and available memory slots.
  //
  class GPUStorageSlot {
    // Utility to align addresses to 256 bytes.
    static uintptr_t AlignPtr(uintptr_t Addr) {
      return (Addr + 256 - 1) / 256 * 256;
    }
    uintptr_t Ptr;
    size_t Size;

  public:
    uintptr_t getPtr() const { return Ptr; }
    uintptr_t getAlignedPtr() const { return AlignPtr(Ptr); }
    size_t getSize() const { return Size; }
    size_t getAlignedSize() const {
      size_t ASize = Size - (getAlignedPtr() - getPtr());
      return (ASize > getSize()) ? 0 : ASize;
    }
    uintptr_t getBegin() const { return Ptr; }
    uintptr_t getEnd() const { return Ptr + Size; }

    GPUStorageSlot(void *Ptr, size_t Size)
        : Ptr(reinterpret_cast<uintptr_t>(Ptr)), Size(Size) {}
    GPUStorageSlot(uintptr_t Ptr, size_t Size) : Ptr(Ptr), Size(Size) {}
  };

  std::list<GPUStorageSlot> GPUStorageUsedSlots; // Sorted by address
  std::list<GPUStorageSlot> GPUStorageFreeSlots; // Sorted by address

  // Malloc and free GPU storage.
  void *GPUStorageMalloc(size_t Size) {
    flush();

    // Candidate iterator.
    auto FE = GPUStorageFreeSlots.end();
    auto FC = FE;

    // Look for the smallest slot where we can fit the storage.
    for (auto FI = GPUStorageFreeSlots.begin(); FI != FE; ++FI) {
      // We got a slot! Get the one with the minimum size.
      if (FI->getAlignedSize() >= Size)
        FC = (FC == FE) ? FI : ((FC->getSize() > FI->getSize()) ? FI : FC);
    }

    // No slots, so we have no memory available.
    if (FC == FE) {
      assert(false && "Run out of GPU memory space.");
      return nullptr;
    }

    // Create a slot in the used list. Find where to insert
    // and emplace it there.
    auto UI = GPUStorageUsedSlots.begin();
    for (auto UE = GPUStorageUsedSlots.end(); UI != UE; ++UI)
      if (UI->getPtr() > FC->getPtr())
        break;

    uintptr_t AlignPtr = FC->getAlignedPtr();
    GPUStorageUsedSlots.emplace(UI, FC->getPtr(),
                                Size + AlignPtr - FC->getPtr());

    // Adjust the slot from the free-slots list.
    if (FC->getAlignedSize() == Size)
      GPUStorageFreeSlots.erase(FC);
    else
      *FC = GPUStorageSlot(AlignPtr + Size, FC->getAlignedSize() - Size);

    return reinterpret_cast<void *>(AlignPtr);
  }

  void GPUStorageFree(void *InPtr) {
    flush();
    uintptr_t Ptr = reinterpret_cast<uintptr_t>(InPtr);

    // Look for the candidate slot the pointer refers too.
    auto UC = GPUStorageUsedSlots.begin();
    auto UE = GPUStorageUsedSlots.end();
    for (; UC != UE; ++UC)
      if (UC->getAlignedPtr() == Ptr)
        break;

    // If we can't find the pointer, we have an invalid pointer to dealloc.
    if (UC == UE) {
      assert(false && "Can't find GPU storage to deallocate.");
      return;
    }

    // Look in the free slots list where the free space should be inserted.
    auto FI = GPUStorageFreeSlots.begin();
    auto FE = GPUStorageFreeSlots.end();
    for (; FI != FE; ++FI)
      if (FI->getPtr() > UC->getPtr())
        break;

    // Check if the slot can be combined with previous free slot.
    auto PrevI = (FI != GPUStorageFreeSlots.begin()) ? std::prev(FI) : FE;
    if (PrevI != FE && UC->getPtr() == PrevI->getEnd()) {
      // Combine with previous slot.
      *UC = GPUStorageSlot(PrevI->getPtr(), PrevI->getSize() + UC->getSize());
      GPUStorageFreeSlots.erase(PrevI);
    }

    // Check if the slot can be combined with following
    if (FI != FE && FI->getPtr() == UC->getEnd()) {
      *UC = GPUStorageSlot(UC->getPtr(), UC->getSize() + FI->getSize());
      auto RemoveI = FI;
      ++FI;
      GPUStorageFreeSlots.erase(RemoveI);
    }

    // Move contents of used list into free list.
    GPUStorageFreeSlots.splice(GPUStorageFreeSlots.begin(), GPUStorageUsedSlots,
                               UC);
    return;
  }

public:
  void initialize(void) {
    CHECK(cudaMallocHost((void **)&ReductionHost, 3 * sizeof(cs_real_t)));
    CHECK(cudaMalloc((void **)&ReductionDev, 6 * sizeof(cs_real_t)));
    CHECK(cudaMemsetAsync(ReductionDev, 0, 6 * sizeof(cs_real_t)))

    // Allocate a chunk of memory to avoid having to repeat malloc/free
    // operations in the GPU. We assume that ranks are evenly distributed
    // by GPUs. We take 80% of the memory available and split it by ranks.
    assert(NumberRanks && NumberNodeDevices);
    size_t RanksUsingSameGPU =
        (NumberRanks + NumberNodeDevices - 1) / NumberNodeDevices;
    size_t MemoryForAllRanks = DeviceProperties.totalGlobalMem * 80ul / 100ul;
    size_t MemoryForCurrentRank = MemoryForAllRanks / RanksUsingSameGPU;

    printf(" --> Allocating %lu/%lu bytes in DeviceID: %d\n",
           MemoryForCurrentRank, DeviceProperties.totalGlobalMem, DeviceID);

    CHECK(cudaMalloc(&GPUStorage, MemoryForCurrentRank));

    assert(GPUStorage);

    // We start we a very empty slot.
    GPUStorageFreeSlots.emplace_back(GPUStorage, MemoryForCurrentRank);
  }
  void finalize(void) {
    if (ReductionHost)
      CHECK(cudaFreeHost(ReductionHost));
    if (ReductionDev)
      CHECK(cudaFree(ReductionDev));
    if (GPUStorage)
      CHECK(cudaFree(GPUStorage));

    printf("Have %d used slots and %d free slots\n", GPUStorageUsedSlots.size(),
           GPUStorageFreeSlots.size());
    printf("Used Slots:\n");
    for (auto &S : GPUStorageUsedSlots)
      printf("---> Used [0x%016lx-0x%016lx[\n", S.getPtr(), S.getEnd());
    printf("Free Slots:\n");
    for (auto &S : GPUStorageFreeSlots)
      printf("---> Free [0x%016lx-0x%016lx[\n", S.getPtr(), S.getEnd());
  }

  // Get reduction results.
  cs_real_t *getReductionDevPtr(void) const {
    assert(ReductionDev && "Reduction results storage is not defined!");
    return ReductionDev;
  }
  cs_real_t *getReductionResults(void) {
    // Copy the data from the GPU from the correct version of the buffer.
    cs_real_t *DevPtr = ReductionDev + (ReductionResultVersion * 3);

    flush();
    CHECK(cudaMemcpy(ReductionHost, DevPtr, 3 * sizeof(cs_real_t),
                     cudaMemcpyDeviceToHost));

    // Toggle version
    ReductionResultVersion ^= 0x01;

    return ReductionHost;
  }
  cs_lnum_t getReductionVersion(void) const { return ReductionResultVersion; }

  // Map/unmap the specified pointer to the device. Return the device pointer or
  // null if the pointer is not valid. All movements are asynchronous unless the
  // client specifies the synchronous flag.
  void map(const void *InputPointer, size_t Size, MapType Type,
           bool Synchronous = false) {
    // static int ID = 0;

    // If we don't have devices we don't need to do anything.
    if (!NumberNodeDevices)
      return;

    char *Pointer =
        const_cast<char *>(reinterpret_cast<const char *>(InputPointer));

    assert(Pointer && "Invalid slot information");

    char *Begin = Pointer;
    char *End = Pointer + Size;

    //
    // Detect the provided pointer in the list of maps.
    //
    auto II = SlotsInfo.begin();
    auto IE = SlotsInfo.end();

    // Zero size slot.
    if (Begin == End)
      for (; II != IE && II->getBegin() != Begin; ++II)
        ;
    // Non-zero size slot.
    else
      for (; II != IE; ++II) {
        // There are no overlaps.
        if (Begin >= II->getEnd() || End <= II->getBegin())
          continue;

        // The slot is contained in existing slot.
        if (II->getBegin() <= Begin && II->getEnd() >= End)
          break;

        // Partial overlap of slot.
        printf("Slot [%p-%p[ partially overlaps with [%p-%p[\n", Begin, End,
               II->getBegin(), II->getEnd());
        assert(false);
      }

    //
    // Allocate/deallocate and move data as needed.
    //

    //
    // Slot exists.
    //
    if (II != IE) {
      if (Type == maptype_alloc || Type == maptype_to) {
        II->incCount();
        dumpSlots();
        return;
      }

      assert((Type == maptype_from || Type == maptype_release) &&
             "Unexpected map type.");

      // Is this the last reference? If, so copy the data back and destroy the
      // slot.
      II->decCount();
      if (!II->getCount()) {
        if (Type == maptype_from) {
          flush();
          if (Synchronous) {
            CHECK(cudaMemcpy(Begin, II->getDevBegin(), II->getSize(),
                             cudaMemcpyDeviceToHost));
          } else {
            CHECK(cudaMemcpyAsync(Begin, II->getDevBegin(), II->getSize(),
                                  cudaMemcpyDeviceToHost));
          }
        }

        GPUStorageFree(II->getDevBegin());
        SlotsInfo.erase(II);
      }
      dumpSlots();
      return;
    }

    //
    // Slot does not exist.
    //

    // Slot does not exist so there is no point in releasing it.
    if (Type == maptype_release) {
      dumpSlots();
      return;
    }

    assert((Type == maptype_alloc || Type == maptype_to) &&
           "Can't move from nonexistent slot.");
    char *DevicePtr = reinterpret_cast<char *>(GPUStorageMalloc(Size));

    assert(DevicePtr && "Malloc returned invalid device pointer.");

    if (Type == maptype_to) {
      if (Synchronous) {
        flush();
        CHECK(cudaMemcpy(DevicePtr, Begin, Size, cudaMemcpyHostToDevice));
      } else {
        CHECK(cudaMemcpyAsync(DevicePtr, Begin, Size, cudaMemcpyHostToDevice));
      }
    }

    SlotsInfo.emplace_front(Begin, End, DevicePtr);

    dumpSlots();
    return;
  }

  // Return the device pointer for the specified host pointer.
  void *getDevicePtr(const void *InputHostPointer,
                     bool MustExist = true) const {
    char *HostPointer =
        const_cast<char *>(reinterpret_cast<const char *>(InputHostPointer));

    assert(HostPointer && "Invalid host pointer!");

    for (auto &S : SlotsInfo) {
      if (HostPointer >= S.getBegin() && HostPointer < S.getEnd()) {
        char *DevPtr = S.getDevBegin();
        assert(DevPtr && "Invalid device pointer!");
        // Return the device pointer. We need to use the same offset the host
        // pointer is using.
        return DevPtr + (HostPointer - S.getBegin());
      }
    }

    assert(!MustExist && "No device pointer found for specified host pointer.");
    return nullptr;
  }

  template <typename T>
  T *getDevicePtr(T *InputHostPointer, bool MustExist = true) const {
    const T *HostPointer = InputHostPointer;
    return reinterpret_cast<T *>(
        getDevicePtr(reinterpret_cast<const void *>(HostPointer), MustExist));
  }
};

MemoryManager MM;

// Number of threads that cooperate in the computation of one vector element.
// Using too little will reduce the number of coalesced accessed while accessing
// the indexes. Using too much will result in wasting resources, and some
// threads won't have anything to do. This number should be around the average
// number of non zero elements in each row of the input matrix. Must divide 32.
#define CUDA_GPU_GROUPSIZE (8)
// Minimum number of elements (rows) a group takes care of. This will influence
// how many blocks will run overall.
#define CUDA_GPU_MIN_ROWS_PER_GROUP (8)
// Number of threads per block. We should have enough threads for the block to
// be busy, but we should not use too much, otherwise it will decrease the
// likelihood of getting more blocks being scheduled for the same SM. Must be a
// multiple of 32.
#define CUDA_GPU_NUMTHD (128)
#define CUDA_GPU_GROUPS (CUDA_GPU_NUMTHD / CUDA_GPU_GROUPSIZE)
#define CUDA_GPU_WARP_SIZE 32

unsigned setBlockAndGridDimensions(cs_lnum_t n_rows, dim3 &blockSize,
                                   dim3 &gridSize) {
  unsigned NumGroups = CUDA_GPU_GROUPS;
  unsigned MinRowsPerGroup = CUDA_GPU_MIN_ROWS_PER_GROUP;

  // Target number of rows per block.
  unsigned RowsPerBlock = MinRowsPerGroup * NumGroups;

  // Required number of blocks - ceil(n_rows/RowsPerBlock).
  gridSize = {(n_rows + RowsPerBlock - 1) / RowsPerBlock, 1, 1};

  // If it does not meet the requirements, we have to adjust the number of
  // rows per group.
  if (gridSize.x > DeviceProperties.maxGridSize[0]) {
    uint64_t AdjustedRowsPerGroup = (((uint64_t)MinRowsPerGroup * gridSize.x) +
                                     DeviceProperties.maxGridSize[0] - 1) /
                                    DeviceProperties.maxGridSize[0];

    RowsPerBlock = AdjustedRowsPerGroup * NumGroups;
    gridSize.x = (n_rows + RowsPerBlock - 1) / RowsPerBlock;
  }

  assert(gridSize.x <= DeviceProperties.maxGridSize[0] &&
         "Error in adjusting number of rows per block.");
  blockSize = {CUDA_GPU_NUMTHD, 1, 1};
  return RowsPerBlock;
}

// Matrix formats used in matrix vector multiplications.
enum KernelKinds {
  // Matrix-vector:
  MV_CSR,
  MV_MSR,
  MV_CSR_No_Diag,
  MV_MSR_No_Diag,
  MV_Seidel,
  MV_Seidel_With_Red,
  // Dot product:
  DP_xx,
  DP_xx_xy,
  DP_xy_yz,
  // Vector operations:
  VO_vc_equal_zero,
  VO_vc_equal_va,
  VO_vc_sub_equal_va,
  VO_vc_mul_equal_va,
  VO_vc_equal_va_mul_vb,
  VO_vc_add_equal_s_mul_vb,
  VO_2x_vc_add_equal_s_mul_vb,
  VO_vc_equal_va_add_s_mul_vb,

  // Invalid kind.
  InvalidKernelKind,
};

const char *getKernelKindName(KernelKinds Kind) {
  switch (Kind) {
  // Matrix-vector:
  case MV_CSR:
    return "MV_CSR";
  case MV_MSR:
    return "MV_MSR";
  case MV_CSR_No_Diag:
    return "MV_CSR_No_Diag";
  case MV_MSR_No_Diag:
    return "MV_MSR_No_Diag";
  case MV_Seidel:
    return "MV_Seidel";
  case MV_Seidel_With_Red:
    return "MV_Seidel_With_Red";
  // Dot product:
  case DP_xx:
    return "DP_xx";
  case DP_xx_xy:
    return "DP_xx_xy";
  case DP_xy_yz:
    return "DP_xy_yz";
  // Vector operations:
  case VO_vc_equal_zero:
    return "VO_vc_equal_zero";
  case VO_vc_equal_va:
    return "VO_vc_equal_va";
  case VO_vc_sub_equal_va:
    return "VO_vc_sub_equal_va";
  case VO_vc_mul_equal_va:
    return "VO_vc_mul_equal_va";
  case VO_vc_equal_va_mul_vb:
    return "VO_vc_equal_va_mul_vb";
  case VO_vc_add_equal_s_mul_vb:
    return "VO_vc_add_equal_s_mul_vb";
  case VO_2x_vc_add_equal_s_mul_vb:
    return "VO_2x_vc_add_equal_s_mul_vb";
  case VO_vc_equal_va_add_s_mul_vb:
    return "VO_vc_equal_va_add_s_mul_vb";
  // Invalid kind.
  case InvalidKernelKind:
    return "InvalidKernelKind";
  }
  return nullptr;
}

// Kernel to perform matrix-vector multiplication.
template <KernelKinds Kind>
__device__ void matrix_vector_multiplication(
    const cs_lnum_t *restrict dev_row_index,
    const cs_lnum_t *restrict dev_col_id, const cs_real_t *restrict dev_val,
    const cs_real_t *restrict dev_d_val, const cs_real_t *restrict dev_x,
    cs_real_t *restrict dev_y, cs_lnum_t n_rows, cs_lnum_t n_cols,
    cs_lnum_t n_rows_per_block) {

  const unsigned bdimx = CUDA_GPU_GROUPSIZE;
  const unsigned bdimy = CUDA_GPU_GROUPS;
  const unsigned tidx = threadIdx.x % CUDA_GPU_GROUPSIZE;
  const unsigned tidy = threadIdx.x / CUDA_GPU_GROUPSIZE;

  // We expect a 2D block, where y is the row a group of threads cooperate on,
  // and x is the ID of each thread working on the same resulting vector
  // element.
  const cs_lnum_t StartRow = n_rows_per_block * blockIdx.x + tidy;
  const cs_lnum_t EndRowExp = StartRow + n_rows_per_block;
  const cs_lnum_t EndRow = (n_rows < EndRowExp) ? n_rows : EndRowExp;

  for (cs_lnum_t ii = StartRow; ii < EndRow; ii += bdimy) {
    unsigned AM = __activemask();
    const cs_lnum_t r0 = dev_row_index[ii];
    const cs_lnum_t r1 = dev_row_index[ii + 1];

    const cs_lnum_t *restrict col_id = dev_col_id + r0;
    const cs_real_t *restrict m_row = dev_val + r0;
    const cs_lnum_t n_cols = r1 - r0;
    cs_real_t sii = 0.0;

    for (cs_lnum_t jj = tidx; jj < n_cols; jj += bdimx)
      if (Kind == MV_MSR || Kind == MV_MSR_No_Diag || Kind == MV_CSR ||
          col_id[jj] != ii)
        sii += (m_row[jj] * dev_x[col_id[jj]]);

    for (cs_lnum_t kk = 1; kk < bdimx; kk *= 2)
      sii += __shfl_down_sync(AM, sii, kk, bdimx);

    if (!tidx) {

      if (Kind == MV_MSR)
        sii += dev_d_val[ii] * dev_x[ii];

      dev_y[ii] = sii;
    }
  }
  return;
}

#if (__CUDA_ARCH__ < 600)
// Atomic double add for older GPUs.
__device__ double oldAtomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

__device__ void _fw_and_bw_lu_gs(const cs_real_t mat[], int db_size,
                                 cs_real_t x[], const cs_real_t b[]) {

  /* forward */
  for (int ii = 0; ii < db_size; ii++) {
    x[ii] = b[ii];
    for (int jj = 0; jj < ii; jj++)
      x[ii] -= x[jj] * mat[ii * db_size + jj];
  }

  /* backward */
  for (int ii = db_size - 1; ii >= 0; ii--) {
    for (int jj = db_size - 1; jj > ii; jj--)
      x[ii] -= x[jj] * mat[ii * db_size + jj];
    x[ii] /= mat[ii * (db_size + 1)];
  }
}

#if (DB_SIZE_MAX > CUDA_GPU_GROUPSIZE)
#   error MSR block kernels need the block size to be compatible with number of CUDA group threads.
#endif
#undef DB_SIZE_MAX
template <KernelKinds Kind, unsigned DiagBlockSize>
__device__ void matrix_vector_seidel_block(
    const cs_lnum_t *restrict dev_row_index,
    const cs_lnum_t *restrict dev_col_id, const cs_real_t *restrict dev_val,
    const cs_real_t *restrict dev_ad_inv, const cs_real_t *restrict dev_ad,
    const cs_real_t *restrict dev_rhs, cs_real_t *restrict dev_vx,
    cs_lnum_t dev_red_version, cs_real_t *restrict dev_res2x,
    cs_lnum_t diag_block_size, cs_lnum_t dev_db0, cs_lnum_t dev_db1,
    cs_lnum_t dev_db2, cs_lnum_t dev_db3, cs_lnum_t n_rows, cs_lnum_t n_cols_in,
    cs_lnum_t n_rows_per_block) {

  cs_real_t *restrict res = dev_res2x + (dev_red_version * 3);
  cs_real_t *restrict res_to_clear = dev_res2x + ((dev_red_version ^ 0x01) * 3);

  const unsigned bdimx = CUDA_GPU_GROUPSIZE;
  const unsigned bdimy = CUDA_GPU_GROUPS;
  const unsigned tidx = threadIdx.x % CUDA_GPU_GROUPSIZE;
  const unsigned tidy = threadIdx.x / CUDA_GPU_GROUPSIZE;

  // Shared array containing the ad_inv coefficients.
  __shared__ cs_real_t
      shared_ad_inv_storage[(DiagBlockSize * DiagBlockSize) * CUDA_GPU_GROUPS];
  cs_real_t *shared_ad_inv =
      shared_ad_inv_storage + tidy * (DiagBlockSize * DiagBlockSize);

  // Each thread will have its local vx for the matrix operations.
  cs_real_t local_vx0[DiagBlockSize];
  cs_real_t local_vx[DiagBlockSize];

  // Arrays to enable coalesced loading of arrays from memory.
  __shared__ cs_real_t shared_vxm1_storage[DiagBlockSize * CUDA_GPU_GROUPS];
  __shared__ cs_real_t shared_ad_storage[DiagBlockSize * CUDA_GPU_GROUPS];
  cs_real_t *shared_vxm1 = shared_vxm1_storage + tidy * DiagBlockSize;
  cs_real_t *shared_ad = shared_ad_storage + tidy * DiagBlockSize;

  // Shared memory to help in block reduction.
  __shared__ cs_real_t shared_res[CUDA_GPU_NUMTHD / CUDA_GPU_WARP_SIZE];

  // We expect a 2D block, where y is the row a group of threads cooperate on,
  // and x is the ID of each thread working on the same resulting vector
  // element.
  const cs_lnum_t StartRow = n_rows_per_block * blockIdx.x + tidy;
  const cs_lnum_t EndRowExp = StartRow + n_rows_per_block;
  const cs_lnum_t EndRow = (n_rows < EndRowExp) ? n_rows : EndRowExp;

  // Zero the buffer for next reduction.
  if (Kind == MV_Seidel_With_Red && !tidx && StartRow < 3)
    res_to_clear[StartRow] = 0.0;

  cs_real_t local_res = 0.0;

  for (cs_lnum_t iit = StartRow; iit < EndRow; iit += bdimy) {
    unsigned AM = __activemask();

    cs_lnum_t ii = iit;
    // Version with reduction uses a reverse scan.
    if (Kind == MV_Seidel_With_Red)
      ii = n_rows - iit - 1;

    const cs_lnum_t r0 = dev_row_index[ii];
    const cs_lnum_t r1 = dev_row_index[ii + 1];

    const cs_lnum_t *restrict col_id = dev_col_id + r0;
    const cs_real_t *restrict m_row = dev_val + r0;
    const cs_lnum_t n_cols = r1 - r0;

    //  Initialize the local vx0 and vxm1 to zero.
    for (auto &v : local_vx0)
      v = 0;

    // Get RHS into vx0.
    for (cs_lnum_t kk = tidx; kk < diag_block_size; kk += bdimx) {
      local_vx0[kk] = dev_rhs[ii * dev_db1 + kk];

      if (Kind == MV_Seidel_With_Red) {
        shared_vxm1[kk] = dev_vx[ii * dev_db1 + kk];
        shared_ad[kk] = dev_ad[ii * dev_db1 + kk];
      }
    }

    // Load ad_inv coefficients.
    for (cs_lnum_t kk = tidx; kk < dev_db3; kk += bdimx)
      shared_ad_inv[kk] = *(dev_ad_inv + dev_db3 * ii + kk);

    // Perform matrix vector operations. dev_vx access is not coalesced.
    for (cs_lnum_t jj = tidx; jj < n_cols; jj += bdimx) {
      // Coalesced accesses for the column and row arrays.
      const auto row = m_row[jj];
      const auto col = col_id[jj];
      for (cs_lnum_t kk = 0; kk < diag_block_size; kk++)
        local_vx0[kk] -= (row * dev_vx[col * dev_db1 + kk]);
    }

    // Make a reduction of local vx0.
    for (cs_lnum_t jj = 1; jj < bdimx; jj *= 2)
      for (cs_lnum_t kk = 0; kk < diag_block_size; kk += 1)
        local_vx0[kk] += __shfl_down_sync(AM, local_vx0[kk], jj, bdimx);

    // Master thread of each group does the forward/backward compute and store
    // the result to memory.
    if (!tidx) {

      _fw_and_bw_lu_gs(shared_ad_inv, dev_db0, local_vx, local_vx0);

      for (cs_lnum_t kk = 0; kk < diag_block_size; kk += 1) {
        if (Kind == MV_Seidel_With_Red) {
          register cs_real_t r =
              shared_ad[kk] * (local_vx[kk] - shared_vxm1[kk]);
          local_res += (r * r);
        }

        dev_vx[ii * dev_db1 + kk] = local_vx[kk];
      }
    }
  }

  // Reduce the result across groups, warps and blocks.
  if (Kind == MV_Seidel_With_Red) {
    // Make a reduction across groups - we have a result for each group.
    for (cs_lnum_t kk = CUDA_GPU_GROUPSIZE; kk < CUDA_GPU_WARP_SIZE; kk *= 2)
      local_res += __shfl_down_sync(0xffffffff, local_res, kk);

    const unsigned numWarps = blockDim.x / CUDA_GPU_WARP_SIZE;
    const unsigned laneID = threadIdx.x % CUDA_GPU_WARP_SIZE;
    const unsigned warpID = threadIdx.x / CUDA_GPU_WARP_SIZE;

    // First lane of each warp record its contribution to shared memory.
    if (!laneID)
      shared_res[warpID] = local_res;

    __syncthreads();

    // We only need the first warp from now on.
    if (warpID)
      return;

    // Each thread in the warp picks one element from shared memory.
    local_res = 0.0;
    if (laneID < numWarps)
      local_res = shared_res[laneID];

    // Make a reduction across the warp, but now we only need to cover numWarps.
    for (cs_lnum_t kk = 1; kk < numWarps; kk *= 2)
      local_res += __shfl_down_sync(0xffffffff, local_res, kk);

    // Write results atomically to memory.
    if (!laneID) {
#if (__CUDA_ARCH__ < 600)
      oldAtomicAdd(res + 0, local_res);
#else
      atomicAdd(res + 0, local_res);
#endif
    }
  }

  return;
}

// Kernel to perform dot products.
template <KernelKinds Kind>
__device__ void
dot_product(cs_lnum_t version, cs_lnum_t n_rows, const cs_real_t *restrict x,
            const cs_real_t *restrict y, const cs_real_t *restrict z,
            cs_real_t *restrict res2x, cs_lnum_t n_rows_per_block) {

  cs_real_t *restrict res = res2x + (version * 3);
  cs_real_t *restrict res_to_clear = res2x + ((version ^ 0x01) * 3);

  // We use 1 shared-memory element per warp.
  __shared__ cs_real_t sh_xx[CUDA_GPU_NUMTHD / CUDA_GPU_WARP_SIZE];
  __shared__ cs_real_t sh_xy[CUDA_GPU_NUMTHD / CUDA_GPU_WARP_SIZE];
  __shared__ cs_real_t sh_yz[CUDA_GPU_NUMTHD / CUDA_GPU_WARP_SIZE];

  // We expect a 1D block.
  const cs_lnum_t StartRow = n_rows_per_block * blockIdx.x + threadIdx.x;
  const cs_lnum_t EndRowExp = StartRow + n_rows_per_block;
  const cs_lnum_t EndRow = (n_rows < EndRowExp) ? n_rows : EndRowExp;

  // Zero the buffer for next reduction.
  if (StartRow < 3)
    res_to_clear[threadIdx.x] = 0.0;

  cs_real_t xx = 0.0, xy = 0.0, yz = 0.0;

  for (cs_lnum_t ii = StartRow; ii < EndRow; ii += blockDim.x) {
    cs_real_t xi, yi, zi;

    xi = x[ii];
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      yi = y[ii];
    if (Kind == DP_xy_yz)
      zi = z[ii];

    if (Kind == DP_xx || Kind == DP_xx_xy)
      xx += xi * xi;
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      xy += xi * yi;
    if (Kind == DP_xy_yz)
      yz += yi * zi;
  }

  // Make a reduction across the warp
  for (cs_lnum_t kk = 1; kk < CUDA_GPU_WARP_SIZE; kk *= 2) {
    if (Kind == DP_xx || Kind == DP_xx_xy)
      xx += __shfl_down_sync(0xffffffff, xx, kk);
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      xy += __shfl_down_sync(0xffffffff, xy, kk);
    if (Kind == DP_xy_yz)
      yz += __shfl_down_sync(0xffffffff, yz, kk);
  }

  const unsigned numWarps = CUDA_GPU_NUMTHD / CUDA_GPU_WARP_SIZE;
  const unsigned laneID = threadIdx.x % CUDA_GPU_WARP_SIZE;
  const unsigned warpID = threadIdx.x / CUDA_GPU_WARP_SIZE;

  // First lane of each warp record its contribution to shared memory.
  if (!laneID) {
    if (Kind == DP_xx || Kind == DP_xx_xy)
      sh_xx[warpID] = xx;
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      sh_xy[warpID] = xy;
    if (Kind == DP_xy_yz)
      sh_yz[warpID] = yz;
  }

  __syncthreads();

  // We only need the first warp from now on.
  if (warpID)
    return;

  // Each thread in the warp picks one element from shared memory.
  xx = 0.0, xy = 0.0, yz = 0.0;
  if (laneID < numWarps) {
    if (Kind == DP_xx || Kind == DP_xx_xy)
      xx = sh_xx[laneID];
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      xy = sh_xy[laneID];
    if (Kind == DP_xy_yz)
      yz = sh_yz[laneID];
  }

  // Make a reduction across the warp, but now we only need to cover numWarps.
  for (cs_lnum_t kk = 1; kk < numWarps; kk *= 2) {
    if (Kind == DP_xx || Kind == DP_xx_xy)
      xx += __shfl_down_sync(0xffffffff, xx, kk);
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      xy += __shfl_down_sync(0xffffffff, xy, kk);
    if (Kind == DP_xy_yz)
      yz += __shfl_down_sync(0xffffffff, yz, kk);
  }

  // Write results atomically to memory.
  if (!laneID) {
#if (__CUDA_ARCH__ < 600)
    if (Kind == DP_xx || Kind == DP_xx_xy)
      oldAtomicAdd(res + 0, xx);
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      oldAtomicAdd(res + 1, xy);
    if (Kind == DP_xy_yz)
      oldAtomicAdd(res + 2, yz);
#else
    if (Kind == DP_xx || Kind == DP_xx_xy)
      atomicAdd(res + 0, xx);
    if (Kind == DP_xx_xy || Kind == DP_xy_yz)
      atomicAdd(res + 1, xy);
    if (Kind == DP_xy_yz)
      atomicAdd(res + 2, yz);
#endif
  }

  return;
}

// Kernel to perform vector operations.
template <KernelKinds Kind>
__device__ void
vector_operation(cs_lnum_t n_rows, cs_real_t s, cs_real_t *restrict vc1,
                 cs_real_t *restrict vc2, const cs_real_t *restrict va,
                 const cs_real_t *restrict vb1, const cs_real_t *restrict vb2,
                 cs_lnum_t n_rows_per_block) {

  // We expect a 1D block.
  const cs_lnum_t StartRow = n_rows_per_block * blockIdx.x + threadIdx.x;
  const cs_lnum_t EndRowExp = StartRow + n_rows_per_block;
  const cs_lnum_t EndRow = (n_rows < EndRowExp) ? n_rows : EndRowExp;

  for (cs_lnum_t ii = StartRow; ii < EndRow; ii += blockDim.x) {
    switch (Kind) {
    case VO_vc_equal_zero:
      vc1[ii] = 0.0;
      break;
    case VO_vc_equal_va:
      vc1[ii] = va[ii];
      break;
    case VO_vc_sub_equal_va:
      vc1[ii] -= va[ii];
      break;
    case VO_vc_mul_equal_va:
      vc1[ii] *= va[ii];
      break;
    case VO_vc_equal_va_mul_vb:
      vc1[ii] = va[ii] * vb1[ii];
      break;
    case VO_vc_add_equal_s_mul_vb:
      vc1[ii] += s * vb1[ii];
      break;
    case VO_2x_vc_add_equal_s_mul_vb:
      vc1[ii] += s * vb1[ii];
      vc2[ii] += s * vb2[ii];
      break;
    case VO_vc_equal_va_add_s_mul_vb:
      vc1[ii] = va[ii] + (s * vb1[ii]);
      break;
    }
  }

  return;
}

#define CS_CUDA_GPU_MAXIMUM_ARGUMENTS_PER_KERNEL 16
#define CS_CUDA_GPU_MAXIMUM_NUM_KERNELS 5

// Class that manages the arguments for a kernel.
class KernelArgsBase {
  // Union to track the size we need to fit all the types we use as argument.
  union KernelArgType {
    cs_lnum_t lnum_t;
    cs_lnum_t *lnum_t_ptr;
    cs_real_t real_t;
    cs_real_t *real_t_ptr;
  };

  const KernelKinds KernelKind;
  size_t NumberOfArgs = 0;
  KernelArgType Args[CS_CUDA_GPU_MAXIMUM_ARGUMENTS_PER_KERNEL];

public:
  template <typename ArgT> int setArg(ArgT Arg) {

    assert(NumberOfArgs < CS_CUDA_GPU_MAXIMUM_ARGUMENTS_PER_KERNEL &&
           "Storing too many args.");

    KernelArgType *Storage = &Args[NumberOfArgs++];
    *reinterpret_cast<ArgT *>(Storage) = Arg;
    return 0;
  }

  KernelArgsBase(const KernelKinds KernelKind) : KernelKind(KernelKind) {}
  KernelArgsBase() : KernelKind(InvalidKernelKind) {}

  template <typename ArgT> __host__ __device__ ArgT getArg(int Index) {

    assert(Index < CS_CUDA_GPU_MAXIMUM_ARGUMENTS_PER_KERNEL &&
           "Storing too many args.");

    KernelArgType *Storage = &Args[Index];
    return *reinterpret_cast<ArgT *>(Storage);
  }

  __host__ __device__ KernelKinds getKind(void) const { return KernelKind; }
};

template <KernelKinds Kind> struct KernelArgs : public KernelArgsBase {
  template <typename... T> KernelArgs(T... InputArgs) : KernelArgsBase(Kind) {
    int expand[] = {0, setArg(InputArgs)...};
    (void)expand;
  }
};

// Class to encode arguments for a series of kernels
struct KernelArgsSeries {
  dim3 GridDim;
  dim3 BlockDim;
  unsigned RowsPerBlock = 0;
  int NumKernels = 0;
  KernelArgsBase Args[CS_CUDA_GPU_MAXIMUM_NUM_KERNELS];

  // Add bundle of arguments.
  template <KernelKinds Kind, typename... T>
  void add(cs_lnum_t n_rows, T... args) {

    assert(NumKernels < CS_CUDA_GPU_MAXIMUM_NUM_KERNELS &&
           "Not expecting to deal with such a long series of kernels.");

    // Create launching dimensions.
    if (!NumKernels)
      RowsPerBlock = setBlockAndGridDimensions(n_rows, BlockDim, GridDim);

    // Set the new arguments in place.
    new (&Args[NumKernels++]) KernelArgs<Kind>(args...);
  }

  // Transfer and execute in the GPU.
  void flush(void);

  // Default Ctor
  KernelArgsSeries() {}

private:
  // Utility function that returns a 64-bit integer ID for a list of kinds. This
  // assumes a kernel can't have more than 5 kinds.
  template <KernelKinds FirstKind = InvalidKernelKind,
            KernelKinds... OtherKinds>
  static constexpr uint64_t getKindsID(uint64_t ID = 0) {
    return (FirstKind == InvalidKernelKind)
               ? ID
               : getKindsID<OtherKinds...>((ID << 12) |
                                           ((unsigned)FirstKind & 0xfff));
  }
};

KernelArgsSeries *KernelArgsSeriesHostVersions = nullptr;
KernelArgsSeries *KernelArgsSeriesHost = nullptr;
// We declare the device side as a char to prevent dynamic initialisation.
__constant__ char KernelArgsSeriesGPU[sizeof(KernelArgsSeries)];

template <KernelKinds Kind>
__device__ int any_kernel(KernelArgsBase &Arg, unsigned n_rows_per_block) {
  switch (Kind) {
  // Matrix-vector:
  case MV_CSR_No_Diag:
  case MV_CSR:
  case MV_MSR_No_Diag:
    matrix_vector_multiplication<Kind>(
        /* dev_row_index */ Arg.getArg<cs_lnum_t *>(0),
        /* dev_col_id */ Arg.getArg<cs_lnum_t *>(1),
        /* dev_val */ Arg.getArg<cs_real_t *>(2),
        /* dev_d_val */ nullptr,
        /* dev_x */ Arg.getArg<cs_real_t *>(3),
        /* dev_y */ Arg.getArg<cs_real_t *>(4),
        /* n_rows */ Arg.getArg<cs_lnum_t>(5),
        /* n_cols */ Arg.getArg<cs_lnum_t>(6),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case MV_MSR:
    matrix_vector_multiplication<Kind>(
        /* dev_row_index */ Arg.getArg<cs_lnum_t *>(0),
        /* dev_col_id */ Arg.getArg<cs_lnum_t *>(1),
        /* dev_val */ Arg.getArg<cs_real_t *>(2),
        /* dev_d_val */ Arg.getArg<cs_real_t *>(3),
        /* dev_x */ Arg.getArg<cs_real_t *>(4),
        /* dev_y */ Arg.getArg<cs_real_t *>(5),
        /* n_rows */ Arg.getArg<cs_lnum_t>(6),
        /* n_cols */ Arg.getArg<cs_lnum_t>(7),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case MV_Seidel:
    matrix_vector_seidel_block<Kind, 3>(
        /* dev_row_index */ Arg.getArg<cs_lnum_t *>(0),
        /* dev_col_id */ Arg.getArg<cs_lnum_t *>(1),
        /* dev_val */ Arg.getArg<cs_real_t *>(2),
        /* dev_ad_inv */ Arg.getArg<cs_real_t *>(3),
        /* dev_ad */ nullptr,
        /* dev_rhs */ Arg.getArg<cs_real_t *>(4),
        /* dev_vx */ Arg.getArg<cs_real_t *>(5),
        /* dev_red_version */ 0,
        /* dev_red */ nullptr,
        /* diag_block_size */ Arg.getArg<cs_lnum_t>(6),
        /* dev_db0 */ Arg.getArg<cs_lnum_t>(7),
        /* dev_db1 */ Arg.getArg<cs_lnum_t>(8),
        /* dev_db2 */ Arg.getArg<cs_lnum_t>(9),
        /* dev_db3 */ Arg.getArg<cs_lnum_t>(10),
        /* n_rows */ Arg.getArg<cs_lnum_t>(11),
        /* n_cols */ Arg.getArg<cs_lnum_t>(12),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case MV_Seidel_With_Red:
    matrix_vector_seidel_block<Kind, 3>(
        /* dev_row_index */ Arg.getArg<cs_lnum_t *>(0),
        /* dev_col_id */ Arg.getArg<cs_lnum_t *>(1),
        /* dev_val */ Arg.getArg<cs_real_t *>(2),
        /* dev_ad_inv */ Arg.getArg<cs_real_t *>(3),
        /* dev_ad */ Arg.getArg<cs_real_t *>(4),
        /* dev_rhs */ Arg.getArg<cs_real_t *>(5),
        /* dev_vx */ Arg.getArg<cs_real_t *>(6),
        /* dev_red_version */ Arg.getArg<cs_lnum_t>(7),
        /* dev_red */ Arg.getArg<cs_real_t *>(8),
        /* diag_block_size */ Arg.getArg<cs_lnum_t>(9),
        /* dev_db0 */ Arg.getArg<cs_lnum_t>(10),
        /* dev_db1 */ Arg.getArg<cs_lnum_t>(11),
        /* dev_db2 */ Arg.getArg<cs_lnum_t>(12),
        /* dev_db3 */ Arg.getArg<cs_lnum_t>(13),
        /* n_rows */ Arg.getArg<cs_lnum_t>(14),
        /* n_cols */ Arg.getArg<cs_lnum_t>(15),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  //
  // Dot product:
  //
  case DP_xx:
    dot_product<Kind>(
        /* version */ Arg.getArg<cs_lnum_t>(0),
        /* n_rows */ Arg.getArg<cs_lnum_t>(1),
        /* x */ Arg.getArg<cs_real_t *>(2),
        /* y */ nullptr,
        /* z */ nullptr,
        /* res */ Arg.getArg<cs_real_t *>(3),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case DP_xx_xy:
    dot_product<Kind>(
        /* version */ Arg.getArg<cs_lnum_t>(0),
        /* n_rows */ Arg.getArg<cs_lnum_t>(1),
        /* x */ Arg.getArg<cs_real_t *>(2),
        /* y */ Arg.getArg<cs_real_t *>(3),
        /* z */ nullptr,
        /* res */ Arg.getArg<cs_real_t *>(4),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case DP_xy_yz:
    dot_product<Kind>(
        /* version */ Arg.getArg<cs_lnum_t>(0),
        /* n_rows */ Arg.getArg<cs_lnum_t>(1),
        /* x */ Arg.getArg<cs_real_t *>(2),
        /* y */ Arg.getArg<cs_real_t *>(3),
        /* z */ Arg.getArg<cs_real_t *>(4),
        /* res */ Arg.getArg<cs_real_t *>(5),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  //
  // Vector operations:
  //
  case VO_vc_equal_zero:
    vector_operation<Kind>(
        /* n_rows */ Arg.getArg<cs_lnum_t>(0),
        /* s */ 0.0,
        /* vc1 */ Arg.getArg<cs_real_t *>(1),
        /* vc2 */ nullptr,
        /* va */ nullptr,
        /* vb1 */ nullptr,
        /* vb2 */ nullptr,
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case VO_vc_equal_va:
  case VO_vc_sub_equal_va:
  case VO_vc_mul_equal_va:
    vector_operation<Kind>(
        /* n_rows */ Arg.getArg<cs_lnum_t>(0),
        /* s */ 0.0,
        /* vc1 */ Arg.getArg<cs_real_t *>(1),
        /* vc2 */ nullptr,
        /* va */ Arg.getArg<cs_real_t *>(2),
        /* vb1 */ nullptr,
        /* vb2 */ nullptr,
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case VO_vc_equal_va_mul_vb:
    vector_operation<Kind>(
        /* n_rows */ Arg.getArg<cs_lnum_t>(0),
        /* s */ 0.0,
        /* vc1 */ Arg.getArg<cs_real_t *>(1),
        /* vc2 */ nullptr,
        /* va */ Arg.getArg<cs_real_t *>(2),
        /* vb1 */ Arg.getArg<cs_real_t *>(3),
        /* vb2 */ nullptr,
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case VO_vc_add_equal_s_mul_vb:
    vector_operation<Kind>(
        /* n_rows */ Arg.getArg<cs_lnum_t>(0),
        /* s */ Arg.getArg<cs_real_t>(1),
        /* vc1 */ Arg.getArg<cs_real_t *>(2),
        /* vc2 */ nullptr,
        /* va */ nullptr,
        /* vb1 */ Arg.getArg<cs_real_t *>(3),
        /* vb2 */ nullptr,
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case VO_2x_vc_add_equal_s_mul_vb:
    vector_operation<Kind>(
        /* n_rows */ Arg.getArg<cs_lnum_t>(0),
        /* s */ Arg.getArg<cs_real_t>(1),
        /* vc1 */ Arg.getArg<cs_real_t *>(2),
        /* vc2 */ Arg.getArg<cs_real_t *>(3),
        /* va */ nullptr,
        /* vb1 */ Arg.getArg<cs_real_t *>(4),
        /* vb2 */ Arg.getArg<cs_real_t *>(5),
        /* n_rows_per_block */ n_rows_per_block);
    break;
  case VO_vc_equal_va_add_s_mul_vb:
    vector_operation<Kind>(
        /* n_rows */ Arg.getArg<cs_lnum_t>(0),
        /* s */ Arg.getArg<cs_real_t>(1),
        /* vc1 */ Arg.getArg<cs_real_t *>(2),
        /* vc2 */ nullptr,
        /* va */ Arg.getArg<cs_real_t *>(3),
        /* vb1 */ Arg.getArg<cs_real_t *>(4),
        /* vb2 */ nullptr,
        /* n_rows_per_block */ n_rows_per_block);
    break;
  }
  __syncthreads();
  return 0;
}

template <KernelKinds... Kinds> __global__ void any_kernels(void) {

  auto *KA = reinterpret_cast<KernelArgsSeries *>(&KernelArgsSeriesGPU[0]);
  const unsigned n_rows_per_block = KA->RowsPerBlock;
  unsigned idx = 0;

  int dummy[] = {any_kernel<Kinds>(KA->Args[idx++], n_rows_per_block)...};
  (void)dummy;
}

void KernelArgsSeries::flush(void) {
  // Nothing to flush.
  if (!NumKernels)
    return;

  assert(KernelArgsSeriesHost == this &&
         "Not expecting to flush this version.");

  CHECK(cudaMemcpyToSymbolAsync(KernelArgsSeriesGPU, this,
                                sizeof(KernelArgsSeries), 0,
                                cudaMemcpyHostToDevice));

  uint64_t ID = 0;
  for (unsigned idx = 0; idx < NumKernels; ++idx)
    ID = (ID << 12) | ((unsigned)Args[idx].getKind() & 0xfff);

  switch (ID) {
    LAUNCH_KERNEL(DP_xx);
    LAUNCH_KERNEL(DP_xx_xy);
    LAUNCH_KERNEL(MV_Seidel);
    LAUNCH_KERNEL(MV_Seidel_With_Red);
    LAUNCH_KERNEL(MV_CSR, VO_vc_sub_equal_va, VO_vc_equal_va_mul_vb,
                  VO_vc_equal_va, DP_xx_xy);
    LAUNCH_KERNEL(MV_CSR, DP_xy_yz);
    LAUNCH_KERNEL(VO_2x_vc_add_equal_s_mul_vb, DP_xx);
    LAUNCH_KERNEL(VO_vc_equal_va_mul_vb, DP_xx_xy);
    LAUNCH_KERNEL(VO_vc_equal_va_add_s_mul_vb);
    LAUNCH_KERNEL(VO_2x_vc_add_equal_s_mul_vb, VO_vc_equal_va_mul_vb, DP_xx_xy);
    LAUNCH_KERNEL(MV_CSR);
    LAUNCH_KERNEL(MV_MSR, VO_vc_sub_equal_va, VO_vc_equal_va_mul_vb,
                  VO_vc_equal_va, DP_xx_xy);
    LAUNCH_KERNEL(MV_MSR, DP_xy_yz);
    LAUNCH_KERNEL(MV_MSR, VO_vc_sub_equal_va);
    LAUNCH_KERNEL(MV_MSR);
    LAUNCH_KERNEL(VO_vc_equal_zero);
    LAUNCH_KERNEL(VO_vc_equal_va_add_s_mul_vb, VO_vc_equal_zero);
    LAUNCH_KERNEL(VO_vc_equal_va, DP_xx_xy);
    LAUNCH_KERNEL(VO_2x_vc_add_equal_s_mul_vb);
  default:
    // couldn't find kernel.
    printf("The implementation for following sequence of kernels could not be "
           "found:\n");
    for (unsigned idx = 0; idx < NumKernels; ++idx)
      printf(" --> %s\n", getKernelKindName(Args[idx].getKind()));
    assert(false);
    break;
  }

  // We are using dual buffering here so that we can reset the structure for the
  // next set of kernel.
  if (this == KernelArgsSeriesHostVersions)
    KernelArgsSeriesHost = KernelArgsSeriesHostVersions + 1;
  else
    KernelArgsSeriesHost = KernelArgsSeriesHostVersions;
  KernelArgsSeriesHost->NumKernels = 0;
}

void flush(void) { KernelArgsSeriesHost->flush(); }
} // namespace

extern "C" {
void cs_cuda_map_alloc(const void *Pointer, size_t Size) {
  MM.map(Pointer, Size, MemoryManager::maptype_alloc);
}
void cs_cuda_map_to(const void *Pointer, size_t Size) {
  MM.map(Pointer, Size, MemoryManager::maptype_to);
}
void cs_cuda_map_from(const void *Pointer, size_t Size) {
  MM.map(Pointer, Size, MemoryManager::maptype_from);
}
void cs_cuda_map_from_sync(const void *Pointer, size_t Size) {
  MM.map(Pointer, Size, MemoryManager::maptype_from, /*Synchronous = */ true);
}
void cs_cuda_map_release(const void *Pointer, size_t Size) {
  MM.map(Pointer, Size, MemoryManager::maptype_release);
}

int cs_cuda_mat_vec_p_l_csr(bool exclude_diag,
                            const cs_lnum_t *restrict row_index,
                            const cs_lnum_t *restrict col_id,
                            const cs_real_t *restrict val,
                            const cs_real_t *restrict x, cs_real_t *restrict y,
                            cs_lnum_t n_rows, cs_lnum_t n_cols) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  // Matrix vector multiplication rows from one block may depend on columns
  // computed by other blocks, so we need to flush here.
  flush();

  cs_cuda_map_to(x, n_cols * sizeof(cs_real_t));
  cs_cuda_map_alloc(y, n_rows * sizeof(cs_real_t));

  if (exclude_diag)
    KernelArgsSeriesHost->add<MV_CSR_No_Diag>(
        n_rows, MM.getDevicePtr(row_index), MM.getDevicePtr(col_id),
        MM.getDevicePtr(val), MM.getDevicePtr(x), MM.getDevicePtr(y), n_rows,
        n_cols);
  else
    KernelArgsSeriesHost->add<MV_CSR>(n_rows, MM.getDevicePtr(row_index),
                                      MM.getDevicePtr(col_id),
                                      MM.getDevicePtr(val), MM.getDevicePtr(x),
                                      MM.getDevicePtr(y), n_rows, n_cols);

  cs_cuda_map_release(x, n_cols * sizeof(cs_real_t));
  cs_cuda_map_from_sync(y, n_rows * sizeof(cs_real_t));
  return 1;
}

int cs_cuda_mat_vec_p_l_msr(bool exclude_diag,
                            const cs_lnum_t *restrict row_index,
                            const cs_lnum_t *restrict col_id,
                            const cs_real_t *restrict x_val,
                            const cs_real_t *restrict d_val,
                            const cs_real_t *restrict x, cs_real_t *restrict y,
                            cs_lnum_t n_rows, cs_lnum_t n_cols) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  // Matrix vector multiplication rows from one block may depend on columns
  // computed by other blocks, so we need to flush here.
  flush();

  cs_cuda_map_to(x, n_cols * sizeof(cs_real_t));
  cs_cuda_map_alloc(y, n_rows * sizeof(cs_real_t));

  if (exclude_diag || !d_val)
    KernelArgsSeriesHost->add<MV_MSR_No_Diag>(
        n_rows, MM.getDevicePtr(row_index), MM.getDevicePtr(col_id),
        MM.getDevicePtr(x_val), nullptr, MM.getDevicePtr(x), MM.getDevicePtr(y),
        n_rows, n_cols);
  else
    KernelArgsSeriesHost->add<MV_MSR>(
        n_rows, MM.getDevicePtr(row_index), MM.getDevicePtr(col_id),
        MM.getDevicePtr(x_val), MM.getDevicePtr(d_val), MM.getDevicePtr(x),
        MM.getDevicePtr(y), n_rows, n_cols);

  cs_cuda_map_release(x, n_cols * sizeof(cs_real_t));
  cs_cuda_map_from_sync(y, n_rows * sizeof(cs_real_t));
  return 1;
}

int cs_cuda_seidel_forward(const cs_lnum_t *restrict row_index,
                           const cs_lnum_t *restrict col_id,
                           const cs_real_t *restrict val,
                           const cs_real_t *restrict ad_inv,
                           const cs_real_t *restrict rhs,
                           cs_real_t *restrict vx, cs_lnum_t diag_block_size,
                           const int *db_size, cs_lnum_t n_rows,
                           cs_lnum_t n_cols) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  // Matrix vector multiplication rows from one block may depend on columns
  // computed by other blocks, so we need to flush here.
  flush();

  cs_cuda_map_to(vx, diag_block_size * n_cols * sizeof(cs_real_t));
  cs_cuda_map_to(rhs, diag_block_size * n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<MV_Seidel>(
      n_rows, MM.getDevicePtr(row_index), MM.getDevicePtr(col_id),
      MM.getDevicePtr(val), MM.getDevicePtr(ad_inv), MM.getDevicePtr(rhs),
      MM.getDevicePtr(vx), diag_block_size, db_size[0], db_size[1], db_size[2],
      db_size[3], n_rows, n_cols);

  flush();

  cs_cuda_map_release(rhs, 0);
  cs_cuda_map_from_sync(vx, 0);

  return 1;
}

int cs_cuda_seidel_backward(
    const cs_lnum_t *restrict row_index, const cs_lnum_t *restrict col_id,
    const cs_real_t *restrict val, const cs_real_t *restrict ad_inv,
    const cs_real_t *restrict ad, const cs_real_t *restrict rhs,
    cs_real_t *restrict vx, cs_real_t *restrict red, cs_lnum_t diag_block_size,
    const int *db_size, cs_lnum_t n_rows, cs_lnum_t n_cols) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  // Matrix vector multiplication rows from one block may depend on columns
  // computed by other blocks, so we need to flush here.
  flush();

  cs_cuda_map_to(vx, diag_block_size * n_cols * sizeof(cs_real_t));
  cs_cuda_map_to(rhs, diag_block_size * n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<MV_Seidel_With_Red>(
      n_rows, MM.getDevicePtr(row_index), MM.getDevicePtr(col_id),
      MM.getDevicePtr(val), MM.getDevicePtr(ad_inv), MM.getDevicePtr(ad),
      MM.getDevicePtr(rhs), MM.getDevicePtr(vx), MM.getReductionVersion(),
      MM.getReductionDevPtr(), diag_block_size, db_size[0], db_size[1],
      db_size[2], db_size[3], n_rows, n_cols);

  flush();

  // When we get the reduction results, a synchronization is imposed, so we
  // don't have to do it for vx.
  cs_cuda_map_release(rhs, 0);
  cs_cuda_map_from(vx, 0);

  cs_real_t *res = MM.getReductionResults();
  *red = res[0];

  return 1;
}

int cs_cuda_dot_product_xx(cs_real_t *restrict xx, const cs_real_t *restrict x,
                           cs_lnum_t n_rows) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(x, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<DP_xx>(n_rows, MM.getReductionVersion(), n_rows,
                                   MM.getDevicePtr(x), MM.getReductionDevPtr());

  // flush();

  cs_real_t *res = MM.getReductionResults();
  *xx = res[0];

  cs_cuda_map_release(x, n_rows * sizeof(cs_real_t));
  return 1;
}

int cs_cuda_dot_product_xx_xy(cs_real_t *restrict xx, cs_real_t *restrict xy,
                              const cs_real_t *restrict x,
                              const cs_real_t *restrict y, cs_lnum_t n_rows) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(x, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(y, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<DP_xx_xy>(n_rows, MM.getReductionVersion(), n_rows,
                                      MM.getDevicePtr(x), MM.getDevicePtr(y),
                                      MM.getReductionDevPtr());

  // flush();

  cs_real_t *res = MM.getReductionResults();
  *xx = res[0];
  *xy = res[1];

  cs_cuda_map_release(x, n_rows * sizeof(cs_real_t));
  cs_cuda_map_release(y, n_rows * sizeof(cs_real_t));
  return 1;
}

int cs_cuda_dot_product_xy_yz(cs_real_t *restrict xy, cs_real_t *restrict yz,
                              const cs_real_t *restrict x,
                              const cs_real_t *restrict y,
                              const cs_real_t *restrict z, cs_lnum_t n_rows) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(x, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(y, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(z, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<DP_xy_yz>(
      n_rows, MM.getReductionVersion(), n_rows, MM.getDevicePtr(x),
      MM.getDevicePtr(y), MM.getDevicePtr(z), MM.getReductionDevPtr());

  // flush();

  cs_real_t *res = MM.getReductionResults();
  *xy = res[1];
  *yz = res[2];

  cs_cuda_map_release(x, n_rows * sizeof(cs_real_t));
  cs_cuda_map_release(y, n_rows * sizeof(cs_real_t));
  cs_cuda_map_release(z, n_rows * sizeof(cs_real_t));
  return 1;
}

int cs_cuda_vector_vc_equal_zero_if_exists(cs_lnum_t n_rows, cs_lnum_t n_elems,
                                           cs_real_t *restrict vc) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  // Only have to do something if vc exists in the device.
  cs_real_t *DevPtr = MM.getDevicePtr(vc, /* MustExist = */ false);
  if (!DevPtr)
    return 0;

  KernelArgsSeriesHost->add<VO_vc_equal_zero>(n_rows,
                                              /* n_rows =           */ n_elems,
                                              /* vc1 =              */ DevPtr);

  // flush();
  return 1;
}
int cs_cuda_vector_vc_equal_va(cs_lnum_t n_rows, cs_real_t *restrict vc,
                               const cs_real_t *restrict va) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_alloc(vc, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<VO_vc_equal_va>(
      n_rows,
      /* n_rows =           */ n_rows,
      /* vc1 =              */ MM.getDevicePtr(vc),
      /* va =               */ MM.getDevicePtr(va));

  // flush();
  cs_cuda_map_release(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from_sync(vc, n_rows * sizeof(cs_real_t));
  return 1;
}
int cs_cuda_vector_vc_sub_equal_va(cs_lnum_t n_rows, cs_real_t *restrict vc,
                                   const cs_real_t *restrict va) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(vc, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(va, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<VO_vc_sub_equal_va>(
      n_rows,
      /* n_rows =           */ n_rows,
      /* vc1 =              */ MM.getDevicePtr(vc),
      /* va =               */ MM.getDevicePtr(va));

  // flush();
  cs_cuda_map_release(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from_sync(vc, n_rows * sizeof(cs_real_t));
  return 1;
}
int cs_cuda_vector_vc_mul_equal_va(cs_lnum_t n_rows, cs_real_t *restrict vc,
                                   const cs_real_t *restrict va) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(vc, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(va, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<VO_vc_mul_equal_va>(
      n_rows,
      /* n_rows =           */ n_rows,
      /* vc1 =              */ MM.getDevicePtr(vc),
      /* va =               */ MM.getDevicePtr(va));

  // flush();
  cs_cuda_map_release(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from_sync(vc, n_rows * sizeof(cs_real_t));
  return 1;
}
int cs_cuda_vector_vc_equal_va_mul_vb(cs_lnum_t n_rows, cs_real_t *restrict vc,
                                      const cs_real_t *restrict va,
                                      const cs_real_t *restrict vb) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(vb, n_rows * sizeof(cs_real_t));
  cs_cuda_map_alloc(vc, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<VO_vc_equal_va_mul_vb>(
      n_rows,
      /* n_rows =           */ n_rows,
      /* vc1 =              */ MM.getDevicePtr(vc),
      /* va =               */ MM.getDevicePtr(va),
      /* vb1 =              */ MM.getDevicePtr(vb));

  // flush();
  cs_cuda_map_release(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_release(vb, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from_sync(vc, n_rows * sizeof(cs_real_t));
  return 1;
}
int cs_cuda_vector_vc_add_equal_s_mul_vb(cs_lnum_t n_rows, cs_real_t s,
                                         cs_real_t *restrict vc1,
                                         const cs_real_t *restrict vb1,
                                         cs_real_t *restrict vc2,
                                         const cs_real_t *restrict vb2) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(vc1, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(vb1, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(vc2, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(vb2, n_rows * sizeof(cs_real_t));

  if (vc2)
    KernelArgsSeriesHost->add<VO_2x_vc_add_equal_s_mul_vb>(
        n_rows,
        /* n_rows =           */ n_rows,
        /* s =                */ s,
        /* vc1 =              */ MM.getDevicePtr(vc1),
        /* vc2 =              */ MM.getDevicePtr(vc2),
        /* vb1 =              */ MM.getDevicePtr(vb1),
        /* vb2 =              */ MM.getDevicePtr(vb2));
  else
    KernelArgsSeriesHost->add<VO_vc_add_equal_s_mul_vb>(
        n_rows,
        /* n_rows =           */ n_rows,
        /* s =                */ s,
        /* vc1 =              */ MM.getDevicePtr(vc1),
        /* vb1 =              */ MM.getDevicePtr(vb1));

  // flush();
  cs_cuda_map_release(vb1, n_rows * sizeof(cs_real_t));
  cs_cuda_map_release(vb2, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from(vc1, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from_sync(vc2, n_rows * sizeof(cs_real_t));
  return 1;
}
int cs_cuda_vector_vc_equal_va_add_s_mul_vb(cs_lnum_t n_rows, cs_real_t s,
                                            cs_real_t *restrict vc,
                                            const cs_real_t *restrict va,
                                            const cs_real_t *restrict vb) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_cuda_map_to(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_to(vb, n_rows * sizeof(cs_real_t));
  cs_cuda_map_alloc(vc, n_rows * sizeof(cs_real_t));

  KernelArgsSeriesHost->add<VO_vc_equal_va_add_s_mul_vb>(
      n_rows,
      /* n_rows =           */ n_rows,
      /* s =                */ s,
      /* vc1 =              */ MM.getDevicePtr(vc),
      /* va =               */ MM.getDevicePtr(va),
      /* vb1 =              */ MM.getDevicePtr(vb));

  // flush();
  cs_cuda_map_release(va, n_rows * sizeof(cs_real_t));
  cs_cuda_map_release(vb, n_rows * sizeof(cs_real_t));
  cs_cuda_map_from_sync(vc, n_rows * sizeof(cs_real_t));
  return 1;
}

int cs_cuda_move_to_device_if_exists(cs_lnum_t n_rows, cs_lnum_t n_elems,
                                     cs_real_t *restrict vc) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_real_t *DevPtr = MM.getDevicePtr(vc, /* MustExist = */ false);
  if (!DevPtr)
    return 0;

  flush();
  CHECK(cudaMemcpyAsync(DevPtr, vc, n_elems * sizeof(cs_real_t),
                        cudaMemcpyHostToDevice));
  return 1;
}

int cs_cuda_move_from_device_if_exists(cs_lnum_t n_rows, cs_lnum_t n_elems,
                                       cs_real_t *restrict vc) {

  if (!(NumberNodeDevices && n_rows > CS_CUDA_GPU_THRESHOLD))
    return 0;

  cs_real_t *DevPtr = MM.getDevicePtr(vc, /* MustExist = */ false);
  if (!DevPtr)
    return 0;

  flush();
  CHECK(cudaMemcpy(vc, DevPtr, n_elems * sizeof(cs_real_t),
                   cudaMemcpyDeviceToHost));
  return 1;
}

void *cs_cuda_host_alloc(size_t size) {
  if (NumberNodeDevices) {
    void *ptr;
    CHECK(cudaMallocHost(&ptr, size));
    assert(ptr && "CUDA host malloc returned a nullptr!");
    return ptr;
  }
  return nullptr;
}
int cs_cuda_host_free(void *ptr) {
  if (NumberNodeDevices) {
    CHECK(cudaFreeHost(ptr));
    return 1;
  }
  return 0;
}

void cs_cuda_initialize(void) {

  int NumberVisibleDevices = 0;
  CHECK(cudaGetDeviceCount(&NumberVisibleDevices));

  if (NumberVisibleDevices) {
    // We assume that the way devices are distributed by socket is taken care of
    // by the scheduler. Therefore we should at most one device. The scheduler
    // must also make sure to use a GPU in the same socket.
    assert(NumberVisibleDevices == 1 &&
           "Expecting to have only one visible device by this rank.");

    DeviceID = 0;

    // Assume OpenMPI is being used.
    // TODO: Find a portable way to do it.
    const char *LocalSizeStr = getenv("OMPI_COMM_WORLD_LOCAL_SIZE");
    assert(LocalSizeStr && "Cannot find OMPI_COMM_WORLD_LOCAL_SIZE in "
                           "environment - expecting OpenMPI is in use.");

    NumberRanks = atoi(LocalSizeStr);
    assert(NumberRanks > 0 && "Invalid number of ranks.");

    // Hardcode number of devices in the node to 4 - MPS sets number of devices
    // to 1 which messes up with out memory allocations.
    const char *NumberNodeDevicesString =
        getenv("CS_NUMBER_OF_GPUS_IN_THE_SYSTEM");
    if (NumberNodeDevices) {
      NumberNodeDevices = atoi(NumberNodeDevicesString);
      assert(NumberNodeDevices > 0 && "No CUDA devices found.");
    } else
      NumberNodeDevices = 4;

    CHECK(cudaSetDevice(DeviceID));
  } else
    return;

  assert(DeviceID >= 0 && "Invalid Device ID.");

  //
  // Sanity check to see if the GPU is working correctly. We copy a random
  // number in and out of the GPU to see if it is working properly.
  //
  int hRin = 0, hRout = 0;
  int *dR = nullptr;

  srand(time(NULL));
  while (hRin == hRout)
    hRin = rand();

  CHECK(cudaMalloc((void **)&dR, sizeof(int)));
  CHECK(cudaMemcpy(dR, &hRin, sizeof(int), cudaMemcpyHostToDevice));
  CHECK(cudaMemcpy(&hRout, dR, sizeof(int), cudaMemcpyDeviceToHost));
  CHECK(cudaFree(dR));

  if (hRin == hRout) {
    printf(" --> Device %d seems to be working properly\n", DeviceID);

    // Get   properties of the target device.
    CHECK(cudaGetDeviceProperties(&DeviceProperties, DeviceID));

    // Initialize persistent storage.
    MM.initialize();

    // Initialize kernel argument buffers. We have two versions so that we can
    // do dual buffering, i.e. move one to the device while resetting the other.
    CHECK(cudaMallocHost((void **)&KernelArgsSeriesHostVersions,
                         2 * sizeof(KernelArgsSeries)));
    new (KernelArgsSeriesHostVersions) KernelArgsSeries[2];
    KernelArgsSeriesHost = KernelArgsSeriesHostVersions;
    assert(KernelArgsSeriesHost && KernelArgsSeriesHostVersions &&
           "Invalid pointers.");

  } else {
    printf(" --> Device %d is NOT working properly\n", DeviceID);
    NumberNodeDevices = 0;
    assert(false);
  }
  return;
}

void cs_cuda_finalize(void) {
  // We only need to clear state if we have any devices.
  if (NumberNodeDevices) {
    MM.finalize();
    if (KernelArgsSeriesHostVersions)
      CHECK(cudaFreeHost(KernelArgsSeriesHostVersions));
  }
}

void cs_cuda_attempt_host_alloc(cs_real_t **ptr, int elems) {
  *ptr = nullptr;
  if (NumberNodeDevices) {
    CHECK(cudaMallocHost(ptr, elems * sizeof(cs_real_t)));
    assert(*ptr && "CUDA host malloc returned a nullptr!");
  } else {
    BFT_MALLOC(*ptr, elems, cs_real_t);
    assert(*ptr && "Malloc returned invalid pointer.");
  }
}

void cs_cuda_attempt_host_free(cs_real_t *ptr) {
  if (NumberNodeDevices) {
    CHECK(cudaFreeHost(ptr));
  } else
    BFT_FREE(ptr);
}

} // extern "C"

#else

extern "C" {
// CUDA malloc/free defaults to a regular malloc/free if CUDA is not configured.
void cs_cuda_attempt_host_alloc(cs_real_t **ptr, int elems) {
  *ptr = nullptr;
  BFT_MALLOC(*ptr, elems, cs_real_t);
  assert(*ptr && "Malloc returned invalid pointer.");
}
void cs_cuda_attempt_host_free(cs_real_t *ptr) { BFT_FREE(ptr); }
} // extern "C"

#endif // HAVE_CUDA_OFFLOAD
