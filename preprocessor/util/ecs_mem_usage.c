/*============================================================================
 * Base memory usage information (System and Library dependent)
 *============================================================================*/

/*
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
*/

/*----------------------------------------------------------------------------*/

#include "ecs_def.h"

/* OS type */

#if defined(__linux__) || defined(__linux) || defined(linux)
#define ECS_OS_Linux

#endif

/*
 * Standard C library headers
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined (ECS_OS_Linux) && defined(HAVE_SYS_STAT_H) \
 && defined(HAVE_SYS_TYPES_H) && defined(HAVE_UNISTD_H) \

#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#elif defined(HAVE_GETRUSAGE)
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#endif

#if defined(HAVE_UNISTD_H) && defined(HAVE_SBRK)
#if defined (ECS_OS_Linux)
#define __USE_MISC 1
#endif
#include <unistd.h>
#endif

#if defined(HAVE_STDDEF_H)
#include <stddef.h>
#endif

/*
 * Optional library and ECS headers
 */

#include "ecs_mem_usage.h"

/*-----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-------------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static int  _ecs_mem_usage_global_initialized = 0;

static size_t _ecs_mem_usage_global_max_pr = 0;

#if defined(USE_SBRK)
static void  *_ecs_mem_usage_global_init_sbrk = NULL;
#endif

#if defined (ECS_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                             && defined(HAVE_SYS_TYPES_H)
static int  _ecs_mem_usage_proc_file_init = 0;
#endif

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

#if defined (ECS_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                             && defined(HAVE_SYS_TYPES_H)

/*!
 * \brief Initialize current process memory use count depending on system.
 */

static void
_ecs_mem_usage_pr_size_init(void)
{
  char  buf[512]; /* should be large enough for "/proc/%lu/status"
                     then beginning of file content */
  int fd;
  size_t  r_size, i;
  _Bool   status_has_peak = false;
  const pid_t  pid = getpid();

  /*
    Under Linux with procfs, one line of the pseudo-file "/proc/pid/status"
    (where pid is the process number) is of the following form:
    VmSize:     xxxx kB
    This line may be the 12th to 13th for a 2.6.x kernel.
    On more recent 2.6.x kernels, another line (the 12th) is of the form:
    VmPeak:     xxxx kB
  */

  if (_ecs_mem_usage_proc_file_init != 0)
    return;

  sprintf(buf, "/proc/%lu/status", (unsigned long) pid);

  fd = open(buf, O_RDONLY);

  if (fd != -1) {

    r_size = read(fd, buf, 512);

    if (r_size > 32) { /* Leave a margin for "VmPeak" or "VmSize:" line */
      r_size -= 32;
      for (i = 0; i < r_size; i++) {
        if (buf[i] == 'V' && strncmp(buf+i, "VmPeak:", 7) == 0) {
          status_has_peak = true;
          break;
        }
      }
      for (i = 0; i < r_size; i++) {
        if (buf[i] == 'V' && strncmp(buf+i, "VmSize:", 7) == 0)
          break;
      }
      /* If VmSize was found, proc file may be used */
      if (i < r_size) {
        if (status_has_peak == true)
          _ecs_mem_usage_proc_file_init = 1;
      }
    }

    (void)close(fd);
  }

  /* If initialization failed for some reason (proc file unavailable or does
     or does not contain the required fields), mark method as unusable */
  if (_ecs_mem_usage_proc_file_init == 0)
    _ecs_mem_usage_proc_file_init = -1;
}

/*!
 * \brief Finalize current process memory use count depending on system.
 */

static void
_ecs_mem_usage_pr_size_end(void)
{
  if (_ecs_mem_usage_proc_file_init != 1)
    return;
}

#else  /* defined (ECS_OS_Linux) && ... */

#define _ecs_mem_usage_pr_size_init()
#define _ecs_mem_usage_pr_size_end()

#endif /* defined (ECS_OS_Linux) && ... */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Initialize memory usage count depending on system.
 *
 * This functions checks if it has already been called, so
 * it is safe to call more than once (though it is not
 * thread-safe). Only the first call is effective.
 */

void
ecs_mem_usage_init(void)
{
  if (_ecs_mem_usage_global_initialized != 0)
    return;

#if defined(USE_SBRK)

  /*
    We use sbrk() to know the size of the heap. This is not of any use
    to guess at allocated memory when some part of the memory may
    be allocated with mmap(), such as with glibc on Linux.
  */

  _ecs_mem_usage_global_init_sbrk = (void *) sbrk(0);

#endif /* (USE_SBRK) */

  _ecs_mem_usage_global_initialized = 1;
}

/*!
 * \brief End memory usage count depending on system.
 */

void
ecs_mem_usage_end(void)
{
  _ecs_mem_usage_pr_size_end();
}

/*!
 * \brief Indicates if ecs_mem_usage_...() functions are initialized.
 *
 * \returns 1 if ecs_mem_usage_init has been called, 0 otherwise.
 */

int
ecs_mem_usage_initialized(void)
{
  return _ecs_mem_usage_global_initialized;
}

/*!
 * \brief Return current process memory use (in kB) depending on system.
 *
 * If the information is not available (depending on availability of
 * non-portable function calls), 0 is returned.
 */

#if defined (ECS_OS_Linux) && defined(HAVE_SYS_STAT_H) \
                           && defined(HAVE_SYS_TYPES_H)

size_t
ecs_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;

  /*
    Under Linux with procfs, one line of the pseudo-file "/proc/pid/status"
    (where pid is the process number) is of the following form:
    VmSize:     xxxx kB
    With more recent kernels, we also have a line of the form:
    VmPeak:     xxxx kB
  */

  {
    if (_ecs_mem_usage_proc_file_init == 0)
      _ecs_mem_usage_pr_size_init();

    if (_ecs_mem_usage_proc_file_init == 1) {

      char  buf[81]; /* should be large enough for "/proc/%lu/status" */
      const pid_t  pid = getpid();

      FILE *fp;
      unsigned long val;
      char *s;

      sprintf(buf, "/proc/%lu/status", (unsigned long) pid);
      fp = fopen(buf, "r");

      if (fp != NULL) {

        int fields_read = 0;

        while (fields_read < 2) {
          s = fgets(buf, 80, fp);
          if (s == NULL)
            break;
          if (strncmp(s, "VmSize:", 7) == 0) {
            sscanf (s + 7, "%lu", &val);
            sys_mem_usage = (size_t) val;
            fields_read += 1;
          }
          else if (strncmp(s, "VmPeak:", 7) == 0) {
            sscanf (s + 7, "%lu", &val);
            if ((size_t) val > _ecs_mem_usage_global_max_pr)
              _ecs_mem_usage_global_max_pr = (size_t) val;
            fields_read += 1;
          }
        }

        fclose(fp);

      }
    }

    _ecs_mem_usage_pr_size_end();
  }

  if (sys_mem_usage > _ecs_mem_usage_global_max_pr)
    _ecs_mem_usage_global_max_pr = sys_mem_usage;

  return sys_mem_usage;
}

#elif defined(USE_SBRK)

size_t
ecs_mem_usage_pr_size(void)
{
  size_t alloc_size = 0;

  if (_ecs_mem_usage_global_initialized) {
    void    *end_addr;

    end_addr = (void *) sbrk(0);

#if defined(HAVE_PTRDIFF_T)
    alloc_size = (size_t)(  (ptrdiff_t)end_addr
                          - (ptrdiff_t)_ecs_mem_usage_global_init_sbrk) / 1024;
#else
    alloc_size = (end_addr - _ecs_mem_usage_global_init_sbrk) / 1024;
#endif

  }

  if (alloc_size > _ecs_mem_usage_global_max_pr)
    _ecs_mem_usage_global_max_pr = alloc_size;

  return alloc_size;
}

#elif defined(HAVE_GETRUSAGE)

size_t
ecs_mem_usage_pr_size(void)
{
  size_t sys_mem_usage = 0;
  struct rusage usage;

  getrusage(RUSAGE_SELF, &usage);

  sys_mem_usage = usage.ru_maxrss / 1024;

  return sys_mem_usage;
}

#else /* Default case */

size_t
ecs_mem_usage_pr_size(void)
{
  return 0;
}

#endif /* ECS_OS_Linux, ... */

/*
 * \brief Return maximum process memory use (in kB) depending on OS.
 *
 * The returned value is the maximum returned by ecs_mem_usage_pr_size()
 * during the program's lifetime. With memory allocations which return
 * memory to the system (such as the GNU glibc on Linux systems),
 * this value will be correct only if allocation is tracked. This should
 * be the case if malloc hooks are used with the glibc allocation
 * functions (ECS library's default configuration/installation option),
 * but may give results lower than the true maximum in other cases.
 */

size_t
ecs_mem_usage_max_pr_size(void)
{
  (void) ecs_mem_usage_pr_size();

  return _ecs_mem_usage_global_max_pr;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
