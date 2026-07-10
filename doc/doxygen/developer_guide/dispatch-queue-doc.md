# `cs_dispatch_queue.h`

`cs_dispatch_queue.h` provides 3 objects that implement a 
[SYCL-like device task management](
https://registry.khronos.org/SYCL/specs/sycl-2020/html/sycl-2020.html) system 
based on 3 objects:

- `cs_dispatch_queue`: The main object to interact with, allows running tasks 
  on the host and the device with dependency management.

- `cs_task`/`cs_host_task`: Allow interacting with tasks launched by 
  `cs_dispatch_queue` by waiting for end of a task or extracting an event 
  to synchronize with.

- `cs_event`: Represents an event triggered at the end of a task for other tasks 
  to synchronize with. It is used as an input of `cs_dispatch_queue` methods to
  express dependencies between tasks.

## `cs_dispatch_queue`

`cs_dispatch_queue` is meant to be used in a way similar to 
[`sycl::queue`](https://github.khronos.org/SYCL_Reference/iface/queue.html).
It can be used to declare parallel tasks on a device that depend on 
(ie. wait for) each other. It holds a context that can be accessed as a member 
named `initializer_context` which is used to initialize each task created by the
`cs_dispatch_queue`.

- **In device mode**, each task is spawned from the `initializer_context` and 
  its own [CUDA stream](
  https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__STREAM.html) for
  asynchronous, parallel execution.

- **In host mode**, tasks are run synchronously.

Tasks will be waited upon at destruction. This behavior ensures resources held 
for a task are not destroyed before its termination.

`cs_device_queue` wraps all of `cs_dispatch_queue`'s algorithms to run them
as separated tasks represented by `cs_task`. An additional `single_task` method
can be used to declare tasks run on the host and represented by `cs_host_task`.

## `cs_task` and `cs_host_task`

- `cs_task` holds resources relative to the execution of a
  task started with `cs_dispatch_queue` (ie. a CUDA stream in device mode). 

- `cs_host_task` is a specialization of `cs_task` that holds additional data
  necessary for the execution of a host task.



## `cs_event`

`cs_event` represent a time point in a task to synchronize with.

- **In device mode**, they are implemented on top of `cudaEvent_t`.

- **In host mode**, they are only used for time measurement as the execution
  is synchronous.

Both `cs_task` and `cs_host_task` provide an operator to generate a `cs_event` 
for synchronization or time measurement.

## Usage example

The system can be used as SYCL (without memory handle based dependency 
management) on top of CUDA:

```cpp
inline void
_cs_sycl_like_example()
{
  // Inspired by:
  // https://enccs.github.io/sycl-workshop/task-graphs-synchronization/#how-to-specify-dependencies

  cs_dispatch_queue Q;
  Q.initializer_context.set_use_gpu(true);

  constexpr std::size_t N = 16 * 1024 * 1024;
  int                  *a, *b;

  CS_MALLOC_HD(a, N, int, CS_ALLOC_HOST_DEVICE_SHARED);
  CS_MALLOC_HD(b, N, int, CS_ALLOC_DEVICE);

  {
    // Start task A asynchronously
    auto task_a =
      Q.parallel_for(N, [=] CS_F_HOST_DEVICE(std::size_t id) { a[id] = 1; });

    // Start task B asynchronously
    auto task_b =
      Q.parallel_for(N, [=] CS_F_HOST_DEVICE(std::size_t id) { b[id] = 2; });

    // Start task C asynchronously, which depends on tasks A and B
    auto task_c =
      Q.parallel_for(N,
                     { task_a, task_b },
                     [=] CS_F_HOST_DEVICE(std::size_t id) { a[id] += b[id]; });

    // Start task D asynchronously, which depends on task C.
    // NB: a single_task runs on the host, not the device.
    int  value  = 3;
    auto task_d = Q.single_task(
      { task_c },
      [=](int *data) -> void {
        for (int i = 1; i < N; i++) {
          data[0] += data[i];
        }

        data[0] /= value;
      },
      a);

    // NB: Tasks are waited upon at destruction so this is not necessary:
    // task_d.wait();
  }

  CS_FREE_HD(a);
  CS_FREE_HD(b);
}
```
