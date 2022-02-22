/*============================================================================
 * Base memory allocation wrappers with optional tracing
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

/*
 * Standard C library headers
 */

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * Optional library and ECS headers
 */

#include "ecs_def.h"
#include "ecs_mem.h"
#include "ecs_mem_usage.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*-----------------------------------------------------------------------------
 * Local type definitions
 *----------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*
 * Structure defining an allocated memory block (for memory tracing)
 */

struct _ecs_mem_block_t {

  void    *p_bloc;  /* Allocated memory block start adress */
  size_t   size;    /* Allocated memory block length */

};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*-----------------------------------------------------------------------------
 * Local macro documentation
 *----------------------------------------------------------------------------*/

/*! \fn ECS_MALLOC(_ptr, _ni, _type)
 * \brief Allocate memory for _ni elements of type _type.
 *
 * This macro calls ecs_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr  pointer to allocated memory.
 * \param [in]  _ni   number of elements.
 * \param [in]  _type element type.
 */

/*! \fn ECS_REALLOC(_ptr, _ni, _type)
 * \brief Reallocate memory for _ni elements of type _type.
 *
 * This macro calls ecs_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 * \param [in]      _ni   number of elements.
 * \param [in]      _type element type.
 */

/*! \fn ECS_FREE(_ptr)
 * \brief Free allocated memory.
 *
 * This macro calls ecs_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 */

/*-----------------------------------------------------------------------------
 * Local macro definitions
 *----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *----------------------------------------------------------------------------*/

static int  _ecs_mem_global_initialized = 0;

static FILE *_ecs_mem_global_file = NULL;

static struct _ecs_mem_block_t  *_ecs_mem_global_block_array = NULL;

static unsigned long  _ecs_mem_global_block_nbr = 0 ;
static unsigned long  _ecs_mem_global_block_max = 512 ;

static size_t  _ecs_mem_global_alloc_cur = 0;
static size_t  _ecs_mem_global_alloc_max = 0;

static size_t  _ecs_mem_global_n_allocs = 0;
static size_t  _ecs_mem_global_n_reallocs = 0;
static size_t  _ecs_mem_global_n_frees = 0;

/*-----------------------------------------------------------------------------
 * Local function definitions
 *----------------------------------------------------------------------------*/

/*
 * Given a character string representing a file name, returns
 * pointer to that part of the string corresponding to the base name.
 *
 * parameters:
 *   file_name:      <-- full name of source file.
 *
 * return:
 *   pointer to part of file name corresponding to base name.
 */

static const char *
_ecs_mem_basename(const char  *file_name)
{
  int i;

  if (file_name == NULL)
    return NULL;

  for (i = strlen(file_name) - 1;
       i > 0 && file_name[i] != ECS_PATH_SEP;
       i--);

  if (file_name[i] == ECS_PATH_SEP)
    i++;

  return (file_name + i);
}

/*
 * Determines values associated with an array representing a
 * long integer
 *
 * parameters:
 *   counter: <-- counter to update.
 *   value:   --> counter values in output unit [out].
 *   unit:    --> counter value unit : ' ', 'k', 'm', 'g', 't', or 'p'
 *                for bytes, Kilobytes, Megabytes, Gigabytes, Terabytes,
 *                or Petabytes [out].
 */

static void
_ecs_mem_size_val(const size_t    counter,
                  unsigned long   value[2],
                  char           *unit)
{
  int i;
  size_t _counter[2];
  const char units[] = {' ', 'k', 'm', 'g', 't', 'p', 'e'};

  for (i = 0, _counter[0] = counter, _counter[1] = 0;
       _counter[0] >= 1024 && i < 6;
       i++) {
    _counter[1] = _counter[0] % 1024;
    _counter[0] /= 1024;
  }

  value[0] = _counter[0];
  value[1] = _counter[1];
  *unit = units[i];
}

/*
 * Memory usage summary.
 */

static void _ecs_mem_summary(FILE  *f)
{
  char unit;
  unsigned long value[2];
  size_t mem_usage;

  if (f == NULL)
    return;

  fprintf(f, "\n\n");
  fprintf(f,
          "Memory allocation summary\n"
          "-------------------------\n\n");

  /* Available memory usage information */

  _ecs_mem_size_val(_ecs_mem_global_alloc_cur, value, &unit);
  fprintf(f,
          "Theoretical current allocated memory:   %8lu.%lu %cB\n",
          value[0], value[1], unit);

  _ecs_mem_size_val(_ecs_mem_global_alloc_max, value, &unit);
  fprintf(f,
          "Theoretical maximum allocated memory:   %8lu.%lu %cB\n",
          value[0], value[1], unit);

  fprintf(f,
          "\n"
          "Number of allocations:   %lu\n"
          "          reallocations: %lu\n"
          "          frees:         %lu\n\n",
          (unsigned long)_ecs_mem_global_n_allocs,
          (unsigned long)_ecs_mem_global_n_reallocs,
          (unsigned long)_ecs_mem_global_n_frees);

  if (ecs_mem_usage_initialized() == 1) {

    /* Maximum measured memory */

    mem_usage = ecs_mem_usage_max_pr_size();
    if (mem_usage > 0) {
      fprintf(f,
              "Maximum program memory measure:  %8lu kB\n",
              (unsigned long)mem_usage);
    }

    /* Current measured memory */

    mem_usage = ecs_mem_usage_pr_size();
    if (mem_usage > 0)
      fprintf(f,
              "Current program memory measure:   %8lu kB\n",
              (unsigned long)mem_usage);
  }

}

/*
 * Return the _ecs_mem_block structure corresponding to a given
 * allocated block.
 *
 * parameters:
 *   p_in: <-- allocated block's start adress.
 *
 * returns:
 *   corresponding _ecs_mem_block structure.
 */

static struct _ecs_mem_block_t *
_ecs_mem_block_info(const void *p_get)
{
  struct _ecs_mem_block_t  *pinfo = NULL;
  unsigned long idx;

  if (_ecs_mem_global_block_array != NULL) {

    for (idx = _ecs_mem_global_block_nbr - 1;
         idx > 0 && (_ecs_mem_global_block_array + idx)->p_bloc != p_get;
         idx--);

    if ((_ecs_mem_global_block_array + idx)->p_bloc != p_get) {
      _ecs_mem_summary(stderr);
      ecs_error(__FILE__, __LINE__, 0,
                _("Adress [%10p] does not correspond to "
                  "the beginning of an allocated block."),
                p_get);
    }
    else {
      pinfo = _ecs_mem_global_block_array + idx;
      assert(p_get == pinfo->p_bloc);
    }

  }

  return pinfo;
}

/*
 * Return the size of a given allocated block.
 *
 * parameters:
 *   p_in: <-- allocated block's start adress.
 *
 * returns:
 *   block size.
 */

static size_t
_ecs_mem_block_size(const void  *p_in)
{
  struct _ecs_mem_block_t *pinfo = _ecs_mem_block_info(p_in);

  if (pinfo != NULL)
    return pinfo->size;
  else
    return 0;
}

/*
 * Fill a _ecs_mem_block_t structure for an allocated pointer.
 */

static void
_ecs_mem_block_malloc(void          *p_new,
                      const size_t   size_new)
{
  struct _ecs_mem_block_t *pinfo;

  assert(size_new != 0);

  if (_ecs_mem_global_block_array == NULL)
    return;

  if (_ecs_mem_global_block_nbr >= _ecs_mem_global_block_max) {

    _ecs_mem_global_block_max *= 2;
    _ecs_mem_global_block_array
      = (struct _ecs_mem_block_t *) realloc(_ecs_mem_global_block_array,
                                            sizeof(struct _ecs_mem_block_t)
                                            * _ecs_mem_global_block_max);

    if (_ecs_mem_global_block_array == NULL) {
      _ecs_mem_summary(stderr);
      ecs_error(__FILE__, __LINE__, errno,
                _("Memory allocation failure"));
      return;
    }

  }

  _ecs_mem_global_block_nbr += 1;

  pinfo = _ecs_mem_global_block_array + _ecs_mem_global_block_nbr - 1;

  /* Start adress and size of allocated block */

  pinfo->p_bloc = p_new;
  pinfo->size   = size_new;
}

/*
 * Update a _ecs_mem_block_t structure for an reallocated pointer.
 */

static void
_ecs_mem_block_realloc(const void    *p_old,
                       void          *p_new,
                       size_t         size_new)
{
  struct _ecs_mem_block_t *pinfo;

  assert(size_new != 0);

  pinfo = _ecs_mem_block_info(p_old);

  if (pinfo != NULL) {
    pinfo->p_bloc = p_new;
    pinfo->size   = size_new;
  }
}

/*
 * Free a _ecs_mem_block_t structure for a freed pointer.
 */

static void
_ecs_mem_block_free(const void *p_free)
{
  struct _ecs_mem_block_t *pinfo, *pmove;
  unsigned long idx;

  if (_ecs_mem_global_block_array == NULL)
    return;

  for (idx = _ecs_mem_global_block_nbr - 1;
       idx > 0 && (_ecs_mem_global_block_array + idx)->p_bloc != p_free;
       idx--);

  if ((_ecs_mem_global_block_array + idx)->p_bloc != p_free) {
    _ecs_mem_summary(stderr);
    ecs_error(__FILE__, __LINE__, 0,
              _("Adress [%10p] does not correspond to "
                "the beginning of an allocated block."),
              p_free);
  }

  else {

    /* We move the contents of the array's final block to the position
       of the freed block, and shorten the array's useful part by one. */

    pinfo = _ecs_mem_global_block_array + idx;
    pmove = _ecs_mem_global_block_array + _ecs_mem_global_block_nbr - 1;
    pinfo->p_bloc = pmove->p_bloc;
    pinfo->size   = pmove->size;

    _ecs_mem_global_block_nbr -= 1;

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Initialize memory handling.
 *
 * This function should be called before any other ecs_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * \param log_file_name name of optional log_file (if NULL, no log).
 */

void
ecs_mem_init(const char *log_file_name)
{
  size_t alloc_size;

  if (_ecs_mem_global_initialized == 1) {
    _ecs_mem_summary(stderr);
    ecs_error(__FILE__, __LINE__, 0,
              _("ecs_mem_init() has already been called"));
  }
  _ecs_mem_global_initialized = 1;

  alloc_size = sizeof(struct _ecs_mem_block_t) * _ecs_mem_global_block_max;

  _ecs_mem_global_block_array
    = malloc(sizeof(struct _ecs_mem_block_t) * _ecs_mem_global_block_max);

  if (_ecs_mem_global_block_array == NULL) {
    _ecs_mem_summary(stderr);
    ecs_error(__FILE__, __LINE__, errno,
              _("Failure to allocate \"%s\" (%lu bytes)"),
              "_ecs_mem_global_block_array", (unsigned long)alloc_size);
    return;
  }

  if (log_file_name != NULL) {

    _ecs_mem_global_file = fopen(log_file_name, "w");

    /*
      If the file could not be opened, we do not abort, as it is not
      absolutely necessary. We silently continue.
      (We could warn the user, but this would require either using
      ecs_printf(), which we prefer to keep independent of the ecs_mem_...()
      functions to avoid evental crossed definitions when user-defined, or
      "warning handling" similar to error handling, with a possibility
      of user-defined warning handlers, as we are not sure if the calling
      code uses stderr (especially in a distributed environment). This
      is probably not worth the bother.
    */

    if (_ecs_mem_global_file == NULL)
      fprintf(stderr,
              _("Failure to open memory log file \"%s\"\n"),
              log_file_name);

  }

  /* Log file header */

  if (_ecs_mem_global_file != NULL) {

    fprintf(_ecs_mem_global_file,
            "       :     FILE NAME              : LINE  :"
            "  POINTER NAME                          : N BYTES   :"
            " (+- N BYTES) : TOTAL BYTES  : [    ADRESS]\n"
            "-------:----------------------------:-------:"
            "----------------------------------------:-----------:"
            "-----------------------------:--------------");

  }

}

/*!
 * \brief End memory handling.
 *
 * This function should be called after all other ecs_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */

void ecs_mem_end(void)
{
  if (_ecs_mem_global_initialized == 0) {
    _ecs_mem_summary(stderr);
    ecs_error(__FILE__, __LINE__, 0,
              _("ecs_mem_end() called before ecs_mem_init()"));
  }
  _ecs_mem_global_initialized = 0;

  if (_ecs_mem_global_file != NULL) {

    unsigned long  non_free = 0;
    struct _ecs_mem_block_t  *pinfo;

    /* Memory usage summary */

    _ecs_mem_summary(_ecs_mem_global_file);

    /* List of non-freed pointers */

    if (_ecs_mem_global_block_array != NULL) {

      fprintf(_ecs_mem_global_file, "List of non freed pointers:\n");

      for (pinfo = _ecs_mem_global_block_array;
           pinfo < _ecs_mem_global_block_array + _ecs_mem_global_block_nbr;
           pinfo++) {

        fprintf(_ecs_mem_global_file,"[%10p]\n", pinfo->p_bloc);
        non_free++;

      }

      fprintf(_ecs_mem_global_file,
              "Number of non freed pointers remaining: %lu\n",
              non_free);

    }

    fclose(_ecs_mem_global_file);
  }

  /* Reset defaults in case of later initialization */

  if (_ecs_mem_global_block_array != NULL) {
    free(_ecs_mem_global_block_array);
    _ecs_mem_global_block_array = NULL;
  }

  _ecs_mem_global_block_nbr   = 0 ;
  _ecs_mem_global_block_max   = 512 ;

  _ecs_mem_global_alloc_cur = 0;
  _ecs_mem_global_alloc_max = 0;

  _ecs_mem_global_n_allocs = 0;
  _ecs_mem_global_n_reallocs = 0;
  _ecs_mem_global_n_frees = 0;

}

/*!
 * \brief Indicates if ecs_mem_...() functions are initialized.
 *
 * \returns 1 if ecs_mem_init has been called, 0 otherwise.
 */

int
ecs_mem_initialized(void)
{
  return _ecs_mem_global_initialized;
}

/*!
 * \brief Allocate memory for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the ecs_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to allocated memory.
 */

void *
ecs_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num)
{
  void       *p_loc;
  size_t      alloc_size = ni * size;

  if (ni == 0)
    return NULL;

  /* Allocate memory and check return */

  p_loc = malloc(alloc_size);

  if (p_loc == NULL) {
    _ecs_mem_summary(stderr);
    ecs_error(file_name, line_num, errno,
              _("Failure to allocate \"%s\" (%lu bytes)"),
              var_name, (unsigned long)alloc_size);
    return NULL;
  }
  else if (_ecs_mem_global_initialized == 0)
    return p_loc;

  /* Memory allocation counting */

  _ecs_mem_global_alloc_cur += alloc_size;

  if (_ecs_mem_global_alloc_max < _ecs_mem_global_alloc_cur)
    _ecs_mem_global_alloc_max = _ecs_mem_global_alloc_cur;

  if (_ecs_mem_global_file != NULL) {
    fprintf(_ecs_mem_global_file, "\n  alloc: %-27s:%6d : %-39s: %9lu",
            _ecs_mem_basename(file_name), line_num,
            var_name, (unsigned long)alloc_size);
    fprintf(_ecs_mem_global_file, " : (+%9lu) : %12lu : [%10p]",
            (unsigned long)alloc_size,
            (unsigned long)_ecs_mem_global_alloc_cur,
            p_loc);
    fflush(_ecs_mem_global_file);
  }

  _ecs_mem_block_malloc(p_loc, alloc_size);

  _ecs_mem_global_n_allocs += 1;

  /* Return pointer to allocated memory */

  return p_loc;
}

/*!
 * \brief Reallocate memory for ni elements of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the ecs_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, ecs_alloc() called).
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to reallocated memory.
 */

void *
ecs_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num)
{
  void      *p_loc;

  long size_diff;
  size_t old_size;
  size_t new_size = ni * size;

  /*
    Behave as ecs_malloc() if the previous pointer is equal to NULL.
    Note that the operation will then appear as a first allocation
    ('alloc') in the _ecs_mem_global_file trace file.
  */

  if (ptr == NULL)
    return ecs_mem_malloc(ni,
                          size,
                          var_name,
                          file_name,
                          line_num);

  /* If the old size equals the new size, nothing needs to be done. */

  old_size = _ecs_mem_block_size(ptr);

  if (new_size == old_size)
    return ptr;

  /*
    We may also simply free memory. Note that in this case, the operation
    appears as 'free' in the _ecs_mem_global_file trace file.
  */

  else if (ni == 0)
    return ecs_mem_free(ptr,
                        var_name,
                        file_name,
                        line_num);

  /* In the final case, we have a true reallocation */

  else {

    size_diff = new_size - old_size;

    p_loc = realloc(ptr, new_size);

    if (p_loc == NULL) {
      _ecs_mem_summary(stderr);
      ecs_error(file_name, line_num, errno,
                _("Failure to reallocate \"%s\" (%lu bytes)"),
                var_name, (unsigned long)new_size);
      return NULL;
    }
    else if (_ecs_mem_global_initialized == 0)
      return p_loc;

    _ecs_mem_global_alloc_cur += size_diff;

    if (size_diff > 0) {
      if (_ecs_mem_global_alloc_max < _ecs_mem_global_alloc_cur)
        _ecs_mem_global_alloc_max = _ecs_mem_global_alloc_cur;
    }

    if (_ecs_mem_global_file != NULL) {
      char sgn = (size_diff > 0) ? '+' : '-';
      fprintf(_ecs_mem_global_file, "\nrealloc: %-27s:%6d : %-39s: %9lu",
              _ecs_mem_basename(file_name), line_num,
              var_name, (unsigned long)new_size);
      fprintf(_ecs_mem_global_file, " : (%c%9lu) : %12lu : [%10p]",
              sgn,
              (unsigned long) ((size_diff > 0) ? size_diff : -size_diff),
              (unsigned long)_ecs_mem_global_alloc_cur,
              p_loc);
      fflush(_ecs_mem_global_file);
    }

    _ecs_mem_block_realloc(ptr, p_loc, new_size);

    _ecs_mem_global_n_reallocs += 1;

    return p_loc;
  }

}

/*!
 * \brief Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the ecs_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, ecs_alloc() called).
 * \param [in] var_name  allocated variable name string
 * \param [in] file_name name of calling source file
 * \param [in] line_num  line number in calling source file
 *
 * \returns NULL pointer.
 */

void *
ecs_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num)
{
  size_t  size_info;

  /* NULL pointer case (non-allocated location) */

  if (ptr == NULL)
    return NULL;

  /* General case (free allocated memory) */

  if (_ecs_mem_global_initialized != 0) {

    size_info = _ecs_mem_block_size(ptr);

    _ecs_mem_global_alloc_cur -= size_info;

    if (_ecs_mem_global_file != NULL) {
      fprintf(_ecs_mem_global_file,"\n   free: %-27s:%6d : %-39s: %9lu",
              _ecs_mem_basename(file_name), line_num,
              var_name, (unsigned long)size_info);
      fprintf(_ecs_mem_global_file, " : (-%9lu) : %12lu : [%10p]",
              (unsigned long)size_info,
              (unsigned long)_ecs_mem_global_alloc_cur,
              ptr);
      fflush(_ecs_mem_global_file);
    }

    _ecs_mem_block_free(ptr);

    _ecs_mem_global_n_frees += 1;
  }

  free(ptr);

  return NULL;
}

/*!
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through ecs_mem_...() (in kB).
 */

size_t
ecs_mem_size_current(void)
{
  return (_ecs_mem_global_alloc_cur / 1024);
}

/*!
 * \brief Return maximum theoretical dynamic memory allocated.
 *
 * \return maximum memory handled through ecs_mem_...() (in kB).
 */

size_t
ecs_mem_size_max(void)
{
  return (_ecs_mem_global_alloc_max / 1024);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
