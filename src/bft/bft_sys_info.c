/*============================================================================
 * Base system information (System and Library dependent)
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004-2006  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*-----------------------------------------------------------------------------*/

#if defined(HAVE_CONFIG_H)
#include "cs_config.h"
#endif

#include "bft_config_defs.h"

/*
 * Standard C library headers
 */

#include <string.h>

#if defined(__linux__)
#include <stdio.h>
#endif

#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)
#  if defined(__uxpv__) && defined(HAVE_SYS_TYPES_H)
#  include <sys/types.h> /* Workaround: missing include on VPP500 */
#  endif
#include <sys/sysinfo.h>
#endif

#if defined HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>
#endif

/*
 * Optional library and BFT headers
 */

#include "bft_sys_info.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-------------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------------
 * Local static strings
 *-----------------------------------------------------------------------------*/

#define BFT_SYS_INFO_STRING_LENGTH 80

static char _bft_sys_info_cpu_string[BFT_SYS_INFO_STRING_LENGTH + 1] = "";

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Return basic available CPU info depending on system.
 *
 * \return Pointer to static string containing CPU info.
 */

#if defined(__linux__)

const char *
bft_sys_info_cpu(void)
{

  FILE *fp;
  char buf[BFT_SYS_INFO_STRING_LENGTH + 1] ; /* Should be large enough for the
                                                /proc/cpuinfo line we use */
  char *s;
  int   i;

  fp = fopen("/proc/cpuinfo", "r");

  if (fp != NULL) {

    s = fgets(buf, BFT_SYS_INFO_STRING_LENGTH, fp);

    while (s != NULL && strncmp(s, "model name", 10) != 0)
      s = fgets(buf, BFT_SYS_INFO_STRING_LENGTH, fp);

    if (s != NULL) {
      for ( ; *s != '\0' && *s != ':' ; s++);
      if (*s == ':')
        s++;
      for ( ; *s != '\0' && *s == ' ' ; s++);
      for (i = strlen(s) - 1;
           i > 0 && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r');
           s[i--] = '\0');
      strcpy(_bft_sys_info_cpu_string, s);
    }

    fclose (fp);

  }

  return _bft_sys_info_cpu_string;
}

#else

const char *
bft_sys_info_cpu(void)
{
#if defined HAVE_SYS_UTSNAME_H

  struct utsname  sys_config;

  if (uname(&sys_config) != -1)
    strncpy(_bft_sys_info_cpu_string, sys_config.machine,
            BFT_SYS_INFO_STRING_LENGTH);

  else
    strcpy(_bft_sys_info_cpu_string, "");

#else /* HAVE_SYS_UTSNAME_H */

  strcpy(_bft_sys_info_cpu_string, "");

#endif /* HAVE_SYS_UTSNAME_H */

  return _bft_sys_info_cpu_string;
}

#endif /* bft_OS*/

/*!
 * \brief Return system memory info depending on system.
 *
 * \return System memory (in kB), or 0 if information unavailable.
 */

#if defined(__linux__)

size_t
bft_sys_info_mem_ram(void)
{
#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)

  struct sysinfo info;

  sysinfo(&info);
  return(info.totalram / 1024);

#else
  return(0);
#endif
}

#else

size_t
bft_sys_info_mem_ram(void)
{
  return(0);
}

#endif /* bft_OS*/

/*!
 * \brief Return swap memory info depending on system.
 *
 * \return Swap memory (in kB), or 0 if information unavailable.
 */

#if defined(__linux__)

size_t
bft_sys_info_mem_swap(void)
{
#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)

  struct sysinfo info;

  sysinfo(&info);
  return(info.totalswap / 1024);

#else
  return(0);
#endif
}

#else

size_t
bft_sys_info_mem_swap(void)
{
  return(0);
}

#endif /* OS type */

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
