//============================================================================
//
//     This file is part of the Code_Saturne CFD tool.
//
//     Copyright (C) 2006-2011 EDF S.A., France
//
//     contact: saturne-support@edf.fr
//
//     The Code_Saturne CFD tool is free software; you can redistribute it
//     and/or modify it under the terms of the GNU General Public License
//     as published by the Free Software Foundation; either version 2 of
//     the License, or (at your option) any later version.
//
//     The Code_Saturne CFD tool is distributed in the hope that it will be
//     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
//     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     You should have received a copy of the GNU General Public License
//     along with the Code_Saturne Kernel; if not, write to the
//     Free Software Foundation, Inc.,
//     51 Franklin St, Fifth Floor,
//     Boston, MA  02110-1301  USA
//
//============================================================================

//============================================================================
// Spawn a child process and establish a connection.
//============================================================================

#include "cs_config.h"

// System headers

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/wait.h>

#if defined(HAVE_POSIX_SPAWN)
#include <spawn.h>
#endif

// Local headers

#include "cfd_proxy_defs.h"
#include "cfd_proxy_comm.h"
#include "cfd_proxy_forward.h"

#include "cfd_proxy_child.h"

//----------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

//============================================================================
//                      Local structure definitions
//============================================================================

struct _cfd_proxy_child_t {

  pid_t                  pid;     // Child pid
  cfd_proxy_comm_t      *comm;    // Associated communicator

};

//============================================================================
//                      Private Function Definitions
//============================================================================

//----------------------------------------------------------------------------
// Wait for child process and print status info
//----------------------------------------------------------------------------

static int
_wait_child(cfd_proxy_child_t  *child,
            bool                verbose)
{
  int status = 0;

  if (verbose) {
    cfd_proxy_printf(_("Waiting for child process (%lu) to finish ..."),
                     (unsigned long)(child->pid));
    cfd_proxy_printf_flush();
  }

  child->pid = waitpid(child->pid, &status, 0);

  if (status != 0) {

    if (verbose)
      cfd_proxy_printf(_(" [error]\n"));

#if defined(WIFEXITED)
    if (WIFEXITED(status))
      cfd_proxy_printf(_("Child exited with status %d\n"),
                       WEXITSTATUS(status));
    else if (WIFSIGNALED(status))
      cfd_proxy_printf(_("Child killed by signal %d\n"), WTERMSIG(status));
    else if (WIFSTOPPED(status))
      cfd_proxy_printf(_("Child stopped by signal %d\n"), WSTOPSIG(status));
#else
    cfd_proxy_printf(_("Child exited or killed with status %d\n"), status);
#endif
    cfd_proxy_printf_flush();
  }
  else if (verbose)
    cfd_proxy_printf(_(" [ok]\n"));

  cfd_proxy_printf_flush();

  return status;
}

//============================================================================
//                      Public Function Definitions
//============================================================================

//----------------------------------------------------------------------------
// Spawn a child process and establish a connection
//
// returns:
//   child process handle in case of success, NULL in case of error.
//----------------------------------------------------------------------------

cfd_proxy_child_t *
cfd_proxy_child_start(const char                   *path,
                      char                   *const argv[restrict],
                      char                   *const envp[restrict],
                      cfd_proxy_comm_type_t         comm_type,
                      int                           comm_verbosity)
{
  int i;
  int retval = 0;

  cfd_proxy_child_t *child = NULL;
  char ** _argv = NULL;
  char *socketopt = NULL;
  char *keyopt = NULL;
  const char socketoptbase[] = "--proxy-socket=";
  const char keyoptbase[] = "--proxy-key=";

  // Check path

  {
    struct stat s;

    if (stat(path, &s) != 0) {
      cfd_proxy_error(__FILE__, __LINE__, errno,
                      _("Impossible to run file: %s"), path);
      return NULL;
    }

    else if (S_ISDIR(s.st_mode)) {
      cfd_proxy_error(__FILE__, __LINE__, 0,
                      _("File is a directory: %s"), path);
      return NULL;
    }

    else if (!((S_IXUSR | S_IXGRP | S_IXOTH) & s.st_mode)) {
      cfd_proxy_error(__FILE__, __LINE__, 0,
                      _("File is not executable: %s"), path);
      return NULL;
    }
  }

  // Initialize communicator

  CFDP_MALLOC(child, 1, cfd_proxy_child_t);

  child->pid = 0;
  child->comm = cfd_proxy_comm_initialize(comm_type, comm_verbosity);

  if (child->comm == NULL) {
    CFDP_FREE(child);
    return NULL;
  }

  // Add communication info to argv

  CFDP_MALLOC(socketopt,
              (  strlen(socketoptbase)
               + strlen(cfd_proxy_comm_get_name(child->comm)) + 1),
             char);
  sprintf(socketopt, "%s%s",
          socketoptbase,
          cfd_proxy_comm_get_name(child->comm));

  CFDP_MALLOC(keyopt, strlen(keyoptbase) + sizeof(int)*8 + 1, char);
  sprintf(keyopt, "%s%d", keyoptbase, cfd_proxy_comm_get_key(child->comm));

  for (i = 0; argv[i] != NULL; i++);

  CFDP_MALLOC(_argv, i + 3, char *);

  for (i = 0; argv[i] != NULL; i++)
    _argv[i] = argv[i];

  _argv[i++] = socketopt;
  _argv[i++] = keyopt;
  _argv[i++] = NULL;

#if defined(HAVE_POSIX_SPAWN)

  posix_spawn_file_actions_t file_actions;
  posix_spawnattr_t attrp;

  posix_spawn_file_actions_init(&file_actions);
  posix_spawnattr_init(&attrp);

  retval = posix_spawn(&(child->pid), path, &file_actions, &attrp,
                       _argv, envp);

  if (retval != 0)
    cfd_proxy_error(__FILE__, __LINE__, retval,
                    _("Error spawning child process %s."), path);

  posix_spawnattr_destroy(&attrp);
  posix_spawn_file_actions_destroy(&file_actions);

#elif defined(HAVE_FORK_EXECVE)

  child->pid = fork();

  if (child->pid < 0)
    cfd_proxy_error(__FILE__, __LINE__, errno,
                    _("Error spawning child process %s."), path);

  if (child->pid == 0) {

    fclose(stderr);

    retval = execve(path, (char *const *)_argv, envp);

    if (retval == -1) {
      cfd_proxy_error(__FILE__, __LINE__, errno,
                      _("Error spawning child process %s."), path);
      exit(EXIT_FAILURE);
    }

  }

#endif

  // Free additional command-line arguments

  CFDP_FREE(_argv);
  CFDP_FREE(keyopt);
  CFDP_FREE(socketopt);

  // Wait for connection

  if (child->pid != 0) {
    retval = cfd_proxy_comm_connect(child->comm,
                                    "CFD_Proxy_comm_socket");
    if (retval != 0) {
      _wait_child(child, false);
      child->comm = cfd_proxy_comm_finalize(child->comm);
      CFDP_FREE(child);
    }
  }

  return child;
}

//----------------------------------------------------------------------------
// End connection with a child process and free associated structure
//----------------------------------------------------------------------------

int
cfd_proxy_child_stop(cfd_proxy_child_t **child)
{
  cfd_proxy_child_t *_child = *child;

  int retval = 0;

  if (_child->comm != NULL)
    _child->comm = cfd_proxy_comm_finalize(_child->comm);

  // Wait for child to finish if not already finished

  retval = _wait_child(_child, true);

  if (_child->pid == -1) {
    cfd_proxy_printf(_("\n"));
    cfd_proxy_error(__FILE__, __LINE__, errno,
                    _("Error waiting for child process."));
    retval = 1;
  }

  // User info

  CFDP_FREE(_child);

  *child = NULL;

  return retval;
}

//----------------------------------------------------------------------------
// Forward all calls from the client and their responses.
//----------------------------------------------------------------------------

void
cfd_proxy_child_forward_all(cfd_proxy_child_t  *child)
{
  if (child->comm != NULL)
    cfd_proxy_forward_all(child->comm);
}

//----------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif /* __cplusplus */

