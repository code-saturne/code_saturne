/*============================================================================
 * Base system information (System and Library dependent)
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <string.h>
#include <time.h>

#if defined(__linux__)
#include <stdio.h>
#endif

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined HAVE_SYS_UTSNAME_H
#include <sys/utsname.h>
#endif

#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)
#include <sys/sysinfo.h>
#endif

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
#include <pwd.h>
#endif

#if defined(__bgq__) && defined(__xlc__)
#include <spi/include/kernel/location.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_printf.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_system_info.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local type definitions
 *============================================================================*/

/*============================================================================
 *  Global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Remove leading and trailing whitespace from a string.
 *
 * parameters:
 *   s <-> string to be cleaned
 *----------------------------------------------------------------------------*/

static void
_string_clean(char  *s)
{
  assert(s != NULL);

  int l = strlen(s);
  int i = l - 1;
  while (i > -1 && (   s[i] == ' ' || s[i] == '\t'
                    || s[i] == '\n' || s[i] == '\r')) {
    s[i--] = '\0';
  }

  i = 0;
  while (i < l && (   s[i] == ' ' || s[i] == '\t'
                   || s[i] == '\n' || s[i] == '\r')) {
    i++;
  }

  if (i > 0) {
    int j = 0;
    while (i <= l)
      s[j++] = s[i++];
  }
}

/*----------------------------------------------------------------------------
 * Return basic available CPU info depending on system.
 *
 * parameters:
 *   cpu_str     --> string with CPU description
 *   cpu_str_max <-- maximum length of string with CPU description
 *----------------------------------------------------------------------------*/

static void
_sys_info_cpu(char      *cpu_str,
              unsigned   cpu_str_max)
{
  strcpy(cpu_str, "");

#if defined(__linux__) && !defined(__ve__)

  {
    unsigned _cpu_str_max = (cpu_str_max > 0) ? cpu_str_max - 1 : 0;

    FILE *fp;
    char *s;
    int   i;

    fp = fopen("/proc/cpuinfo", "r");

    if (fp != NULL) {

      s = fgets(cpu_str, _cpu_str_max, fp);

      while (s != NULL && strncmp(s, "model name", 10) != 0)
        s = fgets(cpu_str, _cpu_str_max, fp);

      if (s != NULL) {
        for ( ; *s != '\0' && *s != ':' ; s++);
        if (*s == ':')
          s++;
        for ( ; *s != '\0' && *s == ' ' ; s++);
        for (i = strlen(s) - 1;
             i > 0 && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r');
             s[i--] = '\0');
      }

      fclose (fp);

    }
  }

#elif defined(HAVE_SYS_UTSNAME_H)

  {
    struct utsname  sys_config;

    if (uname(&sys_config) != -1)
      strncpy(cpu_str, sys_config.machine, cpu_str_max);
  }

#endif /* HAVE_SYS_UTSNAME_H */
}

/*----------------------------------------------------------------------------
 * Return Linux info based on /etc/issue.
 *
 * Only the information prior to a first escape sequence is returned.
 *
 * parameters:
 *   issue_str     --> string with system description
 *   issue_str_max <-- maximum length of string with system description
 *----------------------------------------------------------------------------*/

static void
_sys_info_issue(char      *issue_str,
                unsigned   issue_str_max)
{
  issue_str[0] = '\0';

#if defined(__linux__)

  {
    unsigned _issue_str_max = (issue_str_max > 0) ? issue_str_max - 1 : 0;

    FILE *fp;
    char *s;

    fp = fopen("/etc/issue", "r");

    if (fp != NULL) {

      issue_str[0] = ' ';
      issue_str[1] = '(';

      s = fgets(issue_str + 2, _issue_str_max - 4, fp);

      if (s != NULL) {
        int l = strlen(s);
        for (int i = 0; i < l; i++)
          if (s[i] == '\\') {
            s[i] = '\0';
            l = i;
          }
        _string_clean(issue_str + 2);
        l = strlen(issue_str);
        if (l > 2) {
          issue_str[l] = ')';
          issue_str[l+1] = '\0';
        }
        else /* If no info was kept, empty string */
          issue_str[0] = '\0';
      }
    }

    fclose (fp);

  }

#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print available system information.
 *
 * \param[in]  comm  associated MPI communicator
 * \param[in]  log   if true, standard logging; otherwise, single output
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

static void
_system_info(MPI_Comm comm,
             bool     log)

#else

static void
_system_info(bool  log)

#endif

{
  time_t          date;
  size_t          ram;
  int             log_id;
  int             n_logs = (log) ? 2 : 1;

  cs_log_t logs[] = {CS_LOG_DEFAULT, CS_LOG_PERFORMANCE};

#if defined(HAVE_UNAME)
  struct utsname  sys_config;
#endif
#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
  struct passwd   *pwd_user = NULL;
#endif

#if !defined(PATH_MAX)
#define PATH_MAX 1024
#endif

  char  str_date[81];
  char  str_cpu[81];
  char  str_issue[81];
  char  str_directory[PATH_MAX] = "";

# if defined(HAVE_MPI)
  int mpi_flag = 0;
#endif

#if defined(__bgq__) && defined(__xlc__)
  Personality_t personality;
  Kernel_GetPersonality(&personality, sizeof(personality));
#endif

  /* Date */

  if (   time(&date) == -1
      || strftime(str_date, 80, "%c", localtime(&date)) == 0)
    strcpy(str_date, "");

  /* Working directory */

#if defined(HAVE_GETCWD)
  if (getcwd(str_directory, 1024) == NULL)
    strcpy(str_directory, "");
#endif

  /* Print local configuration */
  /*---------------------------*/

  if (log) {
    for (log_id = 0; log_id < n_logs; log_id++)
      cs_log_printf(logs[log_id],
                    "\n%s\n", _("Local case configuration:\n"));
  }

  for (log_id = 0; log_id < n_logs; log_id++)
    cs_log_printf(logs[log_id],
                  "  %s%s\n", _("Date:                "), str_date);

  /* System and machine */

  _sys_info_issue(str_issue, 81);

#if defined(HAVE_UNAME)

  if (uname(&sys_config) != -1) {
    for (log_id = 0; log_id < n_logs; log_id++) {
      cs_log_printf(logs[log_id],
                    "  %s%s %s%s\n", _("System:              "),
                    sys_config.sysname, sys_config.release,
                    str_issue);
      cs_log_printf(logs[log_id],
                    "  %s%s\n", _("Machine:             "),
                    sys_config.nodename);
    }
  }
#endif

  _sys_info_cpu(str_cpu, 81);

  for (log_id = 0; log_id < n_logs; log_id++)
    cs_log_printf(logs[log_id],
                  "  %s%s\n", _("Processor:           "), str_cpu);

  /* Available memory */

#if defined(__bgq__) && defined(__xlc__)
  ram = personality.DDR_Config.DDRSizeMB;
#elif defined(__linux__) \
   && defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)
  {
    struct sysinfo info;
    sysinfo(&info);
    ram = info.totalram / (size_t)(1024*1024);
  }
#else
  ram = 0; /* TODO: complete for Windows systems */
#endif

  if (ram > 0) {
    for (log_id = 0; log_id < n_logs; log_id++)
      cs_log_printf(logs[log_id],
                    "  %s%llu %s\n", _("Memory:              "),
                    (unsigned long long)ram, _("MB"));
  }

  /* User info */

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)

  /* Functions not available on IBM Blue Gene or Cray XT,
     but a stub may exist, so we make sure we ignore it */
#if   defined(__bg__) || defined(_CRAYC) \
   || defined(__CRAYXT) || defined(__CRAYXE) || defined(__CRAYXC)
  pwd_user = NULL;
#else
  pwd_user = getpwuid(geteuid());
#endif

  if (pwd_user != NULL) {

    size_t l_info = 0;

    cs_log_printf(CS_LOG_DEFAULT,
                  "  %s%s", _("User:                "), pwd_user->pw_name);

    if (pwd_user->pw_gecos != NULL) {
      for (l_info = 0;
           (   pwd_user->pw_gecos[l_info] != '\0'
            && pwd_user->pw_gecos[l_info] != ',');
           l_info++);
      if (pwd_user->pw_gecos[l_info] == ',')
        pwd_user->pw_gecos[l_info] = '\0';
      cs_log_printf(CS_LOG_DEFAULT,
                    " (%s)", pwd_user->pw_gecos);
    }

    cs_log_printf(CS_LOG_DEFAULT, "\n");
  }

#endif /* defined(HAVE_GETPWUID) && defined(HAVE_GETEUID) */

  /* Directory info */

  for (log_id = 0; log_id < n_logs; log_id++)
    cs_log_printf(logs[log_id],
                  "  %s%s\n", _("Directory:           "), str_directory);

  /* MPI Info */

#if defined(HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag != 0) {

    int n_ranks = 1, n_world_ranks = 1;

    MPI_Comm_size(comm, &n_ranks);
    MPI_Comm_size(MPI_COMM_WORLD, &n_world_ranks);

#   if defined(__bgq__) && defined(__xlc__)

    {
      int a_torus, b_torus, c_torus, d_torus, e_torus;
      int n_flags = personality.Network_Config.NetFlags;
      int n_hw_threads = Kernel_ProcessorCount();

      for (log_id = 0; log_id < n_logs; log_id++) {
        cs_log_printf(logs[log_id],
                      "  %s%d\n", _("MPI ranks:           "), n_ranks);
        if (n_world_ranks > n_ranks)
          cs_log_printf(logs[log_id],
                        "  %s%d\n", _("MPI_COMM_WORLD size: "),
                        n_world_ranks);
        cs_log_printf(logs[log_id],
                      "  %s%d\n", _("Hardware threads:    "), n_hw_threads);
      }

      if (n_flags & ND_ENABLE_TORUS_DIM_A) a_torus = 1; else a_torus = 0;
      if (n_flags & ND_ENABLE_TORUS_DIM_B) b_torus = 1; else b_torus = 0;
      if (n_flags & ND_ENABLE_TORUS_DIM_C) c_torus = 1; else c_torus = 0;
      if (n_flags & ND_ENABLE_TORUS_DIM_D) d_torus = 1; else d_torus = 0;
      if (n_flags & ND_ENABLE_TORUS_DIM_E) e_torus = 1; else e_torus = 0;

      for (log_id = 0; log_id < n_logs; log_id++) {
        cs_log_printf(logs[log_id],
                      "  %s<%d,%d,%d,%d,%d>\n", _("Block shape:         "),
                      personality.Network_Config.Anodes,
                      personality.Network_Config.Bnodes,
                      personality.Network_Config.Cnodes,
                      personality.Network_Config.Dnodes,
                      personality.Network_Config.Enodes);
        cs_log_printf(logs[log_id],
                      "  %s<%d,%d,%d,%d,%d>\n", _("Torus links enabled: "),
                      a_torus, b_torus, c_torus, d_torus, e_torus);
      }
    }

#   else

    {
      int appnum = -1;

#     if defined(MPI_VERSION) && (MPI_VERSION >= 2)
      void *attp = NULL;
      int flag = 0;
      MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_APPNUM, &attp, &flag);
      if (flag != 0)
        appnum = *(int *)attp;
#     endif

      for (log_id = 0; log_id < n_logs; log_id++) {
        if (appnum > -1 && log_id == 0)
          cs_log_printf(logs[log_id],
                        "  %s%d (%s %d)\n",
                        _("MPI ranks:           "), n_ranks,
                        _("appnum attribute:"), appnum);
        else
          cs_log_printf(logs[log_id],
                        "  %s%d\n", _("MPI ranks:           "), n_ranks);

        if (n_world_ranks > n_ranks)
          cs_log_printf(logs[log_id],
                        "  %s%d\n", _("MPI_COMM_WORLD size: "),
                        n_world_ranks);
      }
    }

#   endif /* defined(HAVE_MPI) */

  }

#endif /* defined(HAVE_MPI) */

#if defined(HAVE_OPENMP)
  {
    int t_id = omp_get_thread_num();
    if (t_id == 0) {
      for (log_id = 0; log_id < n_logs; log_id++) {
        cs_log_printf(logs[log_id],
                      "  %s%d\n", _("OpenMP threads:      "),
                      omp_get_max_threads());
        if (omp_get_dynamic())
          cs_log_printf(logs[log_id],
                        "  %s\n", _("Dynamic scheduling allowed"));
        cs_log_printf(logs[log_id],
                      "  %s%d\n", _("Processors/node:     "),
                      omp_get_num_procs());
      }
    }
  }
#endif
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print available system information.
 *
 * \param[in]  comm  associated MPI communicator
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
cs_system_info(MPI_Comm comm)

#else

void
cs_system_info(void)

#endif

{
#if defined(HAVE_MPI)
  _system_info(comm, true);
#else
  _system_info(true);
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Print available system information, without additional logging
 *
 * \param[in]  comm  associated MPI communicator
 */
/*----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)

void
cs_system_info_no_log(MPI_Comm comm)

#else

void
cs_system_info_no_log(void)

#endif

{
#if defined(HAVE_MPI)
  _system_info(comm, false);
#else
  _system_info(false);
#endif
}

/*-----------------------------------------------------------------------------*/

END_C_DECLS
