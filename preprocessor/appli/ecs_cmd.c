/*============================================================================
 * Command line option parser and tracking.
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

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

#include "cs_config.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <ctype.h>       /* isdigit()          */
#include <math.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>      /* atoi()             */
#include <string.h>      /* strlen()           */
#include <time.h>        /* time(), strftime() */

#if defined(HAVE_MKDIR)
#include <sys/stat.h>
#include <sys/types.h>
#endif

#if defined(HAVE_DUP2)
#include <unistd.h>
#endif

#if defined(HAVE_STAT)
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

#if defined(HAVE_GETCWD)
#include <unistd.h>
#endif

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
#include <pwd.h>
#endif

#if defined(HAVE_UNAME) || defined(HAVE_SYS_UTSNAME_H)
#include <sys/utsname.h>
#endif

#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)
#  if defined(__uxpv__) && defined(HAVE_SYS_TYPES_H)
#  include <sys/types.h> /* Workaround: missing include on VPP500 */
#  endif
#include <sys/sysinfo.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "ecs_file.h"
#include "ecs_mem.h"

#include "ecs_pre.h"

#if defined(HAVE_CGNS)
#include <cgnslib.h>
#endif

#if defined(HAVE_MED)
#include "ecs_med.h"
#endif

/*----------------------------------------------------------------------------
 * Headers for the current file
 *----------------------------------------------------------------------------*/

#include "ecs_cmd.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Static global variables
 *============================================================================*/

static char _sys_info_cpu_string[81] = "";

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*!
 * \brief Return basic available CPU info depending on system.
 *
 * \return Pointer to static string containing CPU info.
 */

#if defined(__linux__) || defined(__linux) || defined(linux)

static const char *
_sys_info_cpu(void)
{
  FILE *fp;
  char buf[81] ; /* Should be large enough for the
                    /proc/cpuinfo line we use */
  int   i;
  char *s = NULL;

  _sys_info_cpu_string[0] = '\0';

  fp = fopen("/proc/cpuinfo", "r");

  if (fp != NULL) {

    s = fgets(buf, 80, fp);

    while (s != NULL && strncmp(s, "model name", 10) != 0)
      s = fgets(buf, 80, fp);

    fclose (fp);
  }

  if (s != NULL) {
    for ( ; *s != '\0' && *s != ':' ; s++);
    if (*s == ':')
      s++;
    for ( ; *s != '\0' && *s == ' ' ; s++);
    for (i = strlen(s) - 1;
         i > 0 && (s[i] == ' ' || s[i] == '\n' || s[i] == '\r');
         s[i--] = '\0');
    strcpy(_sys_info_cpu_string, s);
  }

  return _sys_info_cpu_string;
}

#else

static const char *
_sys_info_cpu(void)
{
  _sys_info_cpu_string[0] = '\0';

#if defined HAVE_SYS_UTSNAME_H

  struct utsname  sys_config;

  if (uname(&sys_config) != -1)
    strncpy(_sys_info_cpu_string, sys_config.machine, 80);

  else
    strcpy(_sys_info_cpu_string, "");

#else /* HAVE_SYS_UTSNAME_H */

  strcpy(_sys_info_cpu_string, "");

#endif /* HAVE_SYS_UTSNAME_H */

  return _sys_info_cpu_string;
}

#endif /* ECS_ARCH */

/*----------------------------------------------------------------------------
 * Print version number
 *----------------------------------------------------------------------------*/

static void
_print_version(void)
{

  char str_date [ECS_STR_SIZE];

  int               ind_mois;
  int               nb_extensions = 0;
  struct tm         time_cnv;

  const char nom_mois[12][4]
    = {"Jan", "Feb", "Mar", "Apr", "May", "Jun",
       "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  strcpy(str_date, ecs_glob_build_date);

  /* Date de compilation */

  for (ind_mois = 0; ind_mois < 12; ind_mois++) {
    if (strncmp(str_date, nom_mois[ind_mois], 3) == 0) {
      time_cnv.tm_mon = ind_mois;
      break;
    }
  }

  sscanf(str_date + 3, "%d", &(time_cnv.tm_mday));
  sscanf(str_date + 6, "%d", &(time_cnv.tm_year));

  time_cnv.tm_year -= 1900;

  strcpy(str_date, __TIME__);

  sscanf(str_date, "%d", &(time_cnv.tm_hour));
  sscanf(str_date + 3, "%d", &(time_cnv.tm_min));
  sscanf(str_date + 6, "%d", &(time_cnv.tm_sec));

  time_cnv.tm_isdst = -1;

  /* Recomputed and internationalized build date */

  mktime(&time_cnv);
  strftime(str_date, ECS_STR_SIZE - 1, "%c", &time_cnv);

#if defined (PACKAGE_VERSION)
  printf(_("\n  %s version %s   (built %s)\n\n"),
         PACKAGE_NAME, PACKAGE_VERSION, str_date);
#else
  printf(_("\n  (built %s)\n\n"), str_date);
#endif

#if defined(HAVE_CCM)
  printf(_("  STAR-CCM+ file format support\n"));
  nb_extensions++;
#endif

#if defined(HAVE_CGNS)

# if defined(CGNS_VERSION)
  int cgns_ver_maj  = CGNS_VERSION/1000;
  int cgns_ver_min  = (CGNS_VERSION % 1000) / 100;
  int cgns_ver_rel  = (CGNS_VERSION % 100) / 10;
  printf(_("  CGNS %d.%d.%d file format support\n"),
         cgns_ver_maj, cgns_ver_min, cgns_ver_rel);
# elif defined(NofValidAreaTypes)
  printf(_("  CGNS %d.%d file format support\n"), 2, 2);
# else
  printf(_("  CGNS %d.%d file format support\n"), 2, 1);
# endif

  nb_extensions++;

#endif

#if defined(HAVE_MED)
  ecs_med__version_shlib();
  if (ecs_glob_med_ver_rel < 0)
    printf(_("  MED %d.%d (HDF5 %d.%d.%d) file format support\n"),
           ecs_glob_med_ver_maj, ecs_glob_med_ver_min,
           ecs_glob_hdf5_ver_maj, ecs_glob_hdf5_ver_min,
           ecs_glob_hdf5_ver_rel);
  else
    printf(_("  MED %d.%d.%d (HDF5 %d.%d.%d) file format support\n"),
           ecs_glob_med_ver_maj, ecs_glob_med_ver_min,
           ecs_glob_med_ver_rel, ecs_glob_hdf5_ver_maj,
           ecs_glob_hdf5_ver_min, ecs_glob_hdf5_ver_rel);
  nb_extensions++;
#endif

  if (ecs_file_version_zlib() != NULL) {
    printf(_("  Reading of compressed files ('.gz') with Zlib %s\n"),
           ecs_file_version_zlib());
    if (strcmp(ecs_file_version_zlib(),
               ecs_file_version_build_zlib()) != 0)
    printf(_("    (compiled with Zlib %s)\n"),
           ecs_file_version_build_zlib());
    nb_extensions++;
  }
  if (nb_extensions > 0)
    printf("\n");

}

/*----------------------------------------------------------------------------
 * Print command line, title, and version
 *----------------------------------------------------------------------------*/

static void
_print_preamble(int    argc,
                char  *argv[])
{
  int  iarg, ltot;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ltot = 0;

  for (iarg = 0; iarg < argc; iarg++) {

    ltot += strlen(argv[iarg]) + 1;

    if (ltot > 80) {
      printf("\n");
      ltot = strlen(argv[iarg]) + 1;
    }

    printf("%s ", argv[iarg]);

  }

  printf(_("\n\n"
           "  .------------------------------.\n"
           "  |                              |\n"
           "  |   Code_Saturne Preprocessor  |\n"
           "  |                              |\n"
           "  `------------------------------'\n"));

  _print_version();
}

/*----------------------------------------------------------------------------
 * Utilitarian function for printing
 *----------------------------------------------------------------------------*/

static void
_fct_prt(const char  *opt,
         const char  *arg,
         const char  *texte)
{
  size_t l = strlen(opt);
  size_t pad = (l > 12) ? 16+12-l :16;

  printf("  %-12s ", opt);
  ecs_print_padded_str(arg, pad);
  printf(" %s\n", texte);
}

/*----------------------------------------------------------------------------
 * Print usage
 *----------------------------------------------------------------------------*/

static void
_print_help(void)
{
  char opt_str[81];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\nUsage:  %s [<options>] <file>\n"),
         "cs_preprocess");

  /* General options */
  /*-----------------*/

  printf(_("\n\nGeneral options:\n\n"));

  _fct_prt("-h", "", _(": this help message"));
  _fct_prt("--help", "", _(": same"));

  printf("\n");

  _fct_prt("--log", _("[file]"),
           _(": redirect terminal output to a file"));

  sprintf(opt_str, _("  (default file: \"%s\")"),
          "preprocessor.log");
  _fct_prt("", "", opt_str);

  printf("\n");

  _fct_prt("--no-write", "",
           _(": do not write preprocessor output"));

  printf("\n");

  _fct_prt("-o", _("<file>"),
           _(": output file name"));
  sprintf(opt_str, _("  (default file: \"%s\")"),
          "mesh_input");
  _fct_prt("", "", opt_str);

  _fct_prt("--out", _("<file>"), _(": same"));

  printf("\n");

  _fct_prt("--reorient", "",
           _(": if necessary, correct orientation of"));
  _fct_prt("", "", _("  cells and faces"));

  printf("\n");

  _fct_prt("--version", "",
           _(": print version number"));

  /* Mesh selection options */
  /*------------------------*/

  printf(_("\n\nMesh selection options:\n\n"));

  _fct_prt("--format", _("<keyword>"),
           _(": selection of mesh file format"));

  _fct_prt("--num", "<n> [...]",
           _(": selection of mesh numbers in file"));
  _fct_prt("", "",
           _("  (if the format allows it)"));

  _fct_prt("--grp-cel", _("<keyword>"),
           _(": add groups of cells"));
  _fct_prt("", "", _("   * based on sections: keyword \"section\""));
  _fct_prt("", "", _("   * based on zones:    keyword \"zone\""));
  _fct_prt("", "", _("  (based on format features/conventions)"));

  _fct_prt("--grp-fac", _("<keyword>"),
           _(": add groups of faces"));
  _fct_prt("", "", _("   * based on sections: keyword \"section""\""));
  _fct_prt("", "", _("   * based on zones:    keyword \"zone""\""));
  _fct_prt("", "", _("  (based on format features/conventions)"));

  printf(_("\n\nAvailable mesh formats:\n"));
  printf(_("                                  extension:    keyword:\n"));

  ecs_pre__aff_formats();

  /* Post-processing options */
  /*-------------------------*/

  printf(_("\n\nPostprocessing options:\n\n"));

  _fct_prt("--case",
           _("<name>"),
           _(": case name (without this option,"));

  sprintf(opt_str,
          _("  the default name is: \"%s\""),
          "preprocess");
  _fct_prt("", "", opt_str);

  printf("\n");

  _fct_prt("--post-error", _("[format]"),
           _(": select output format of error meshes"));

  _fct_prt("--post-volume", _("[format]"),
           _(": activate output of volume meshes"));

  printf(_("\n\nAvailable output formats:\n"));
  printf(_("                                  keyword:\n"));
#if defined(HAVE_CGNS)
  printf("   CGNS                           cgns\n");
#endif
  printf("   EnSight Gold                   ensight\n");
#if defined(HAVE_MED)
  printf("   MED                            med\n");
#endif

  /* Environment variables */
  /*-----------------------*/

  printf("\n\n%s:\n\n",
         _("Environment variables"));

  printf(_("  CS_PREPROCESS_MIN_EDGE_LEN    : "
           "length under which an edge is considered\n"
           "                                  "
           "degenerate (default: 1.e-15)\n\n"));

  printf(_("  CS_PREPROCESS_MEM_LOG         : "
           "name of memory operations trace file\n\n"));

  printf(_("  CS_PREPROCESS_IGNORE_IDEAS_COO_SYS : "
           "ignore I-deas coordinate systems\n\n"));
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche un message indiquant
 *   qu'une meme option est mentionnee au moins 2 fois
 *----------------------------------------------------------------------------*/

static void
ecs_loc_cmd__aff_opt_en_double(int          argc,
                               char        *argv[],
                               const char  *opt)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  _print_preamble(argc, argv);

  _print_help();

  ecs_error(__FILE__, __LINE__, 0,
            _("Error in command line specification.\n\n"
              "Option \"%s\" is set at least twice."), opt);
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche un message indiquant
 *   qu'il manque un argument a une option de la ligne de commande
 *----------------------------------------------------------------------------*/

static void
ecs_loc_cmd__aff_manque_arg(int          argc,
                            char        *argv[],
                            const char  *opt)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  _print_preamble(argc, argv);

  _print_help();

  ecs_error(__FILE__, __LINE__, 0,
            _("Error in command line specification.\n\n"
              "Option \"%s\" requires an argument."), opt);
}

/*----------------------------------------------------------------------------
 *  Fonction qui lit les sous-options d'un post traitement
 *----------------------------------------------------------------------------*/

static void
_read_post_opt(int    argc,
               char  *argv[],
               char   post_type[8],
               int   *argpos)
{
  int iarg = *argpos + 1;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  strcpy(post_type, "ens"); /* default */

  if (iarg < argc) {

    if (!strcmp ("ensight", argv[iarg]))
      *argpos += 1;

    if (!strcmp ("cgns", argv[iarg])) {
#if defined(HAVE_CGNS)
      strcpy(post_type, "cgns");
      *argpos += 1;
#else
      ecs_error(__FILE__, __LINE__, 0,
                _("CGNS output format not available in this build."));
#endif
    }

    if (!strcmp ("med", argv[iarg])) {
#if defined(HAVE_MED)
      strcpy(post_type, "med");
      *argpos += 1;
#else
      ecs_error(__FILE__, __LINE__, 0,
                _("MED output format not available in this build."));
#endif
    }

  }
}

/*----------------------------------------------------------------------------
 *  Initialize command-line options
 *----------------------------------------------------------------------------*/

static ecs_cmd_t *
_cmd_initialize(void)
{
  size_t   i, lng;
  ecs_cmd_t  *cmd;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ECS_MALLOC(cmd, 1, ecs_cmd_t);

  /* Initialisations */
  /*=================*/

  cmd->fic_maillage                  = NULL;

  lng = strlen("preprocess") + 1;
  ECS_MALLOC(cmd->nom_cas, lng, char);
  strcpy(cmd->nom_cas, "preprocess");

  cmd->nom_out                       = NULL;

  cmd->nbr_dump = 0;

  cmd->n_num_maillage = 0;
  cmd->num_maillage = NULL;
  cmd->fmt_maillage = ECS_PRE_FORMAT_NUL;
  cmd->grp_cel_section = false;
  cmd->grp_cel_zone    = false;
  cmd->grp_fac_section = false;
  cmd->grp_fac_zone    = false;

  for (i = 0; i < 8; i++) {
    cmd->post_err[i] = '\0';
    cmd->post_vol[i] = '\0';
  }

  strcpy(cmd->post_err, "ens"); /* default */

  cmd->correct_orient = false;

  return cmd;
}

/*----------------------------------------------------------------------------
 *  Fonction qui affiche la configuration du cas de lancement
 *----------------------------------------------------------------------------*/

static void
ecs_loc_cmd__aff_config(ecs_cmd_t  *cmd)
{
  time_t           date;

#if !defined(PATH_MAX)
#define PATH_MAX 1024
#endif

  char             str_date     [ECS_STR_SIZE]  = "";
  char             str_system   [ECS_STR_SIZE]  = "";
  char             str_machine  [ECS_STR_SIZE]  = "";
  char             str_ram      [ECS_STR_SIZE]  = "";
  char             str_user     [ECS_STR_SIZE]  = "";
  char             str_directory[PATH_MAX] = "";

  size_t           ram = 0;
  int              l_user;

#if defined(HAVE_UNAME)
  struct utsname   sys_config;
#endif

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
  struct passwd   *pwd_user = NULL;
#endif

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /* Determination de la configuration du cas de lancement */
  /*-------------------------------------------------------*/

  /* Date */

  if (time(&date) == -1 ||
      strftime(str_date, ECS_STR_SIZE - 1, "%c", localtime(&date)) == 0)
    strcpy(str_date, "");

  /* Systeme et machine */

#if defined(HAVE_UNAME)
  if (uname(&sys_config) != -1) {
    strcpy(str_system, sys_config.sysname);
    strcat(str_system, " ");
    strcat(str_system, sys_config.release);
    strcpy(str_machine, sys_config.nodename);
  }
#endif

  /* Nom de login */

#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)

  int l_info;

  pwd_user = NULL;

  pwd_user = getpwuid(geteuid());

  if (pwd_user != NULL) {

    str_user[ECS_STR_SIZE - 1] = '\0';
    strncpy(str_user, pwd_user->pw_name, ECS_STR_SIZE - 1);

    if (pwd_user->pw_gecos != NULL) {

      l_user = strlen(str_user);
      for (l_info = 0;
           (   pwd_user->pw_gecos[l_info] != '\0'
            && pwd_user->pw_gecos[l_info] != ',');
           l_info++);

      if (l_user + l_info + 3 < ECS_STR_SIZE) {
        strcat(str_user, " (");
        strncpy(str_user + l_user + 2, pwd_user->pw_gecos, l_info);
        str_user[l_user + 2 + l_info]     = ')';
        str_user[l_user + 2 + l_info + 1] = '\0';
      }

    }
  }

#endif /* defined(HAVE_GETPWUID) && defined(HAVE_GETEUID) */

  /* Mémoire vive */

#if defined(HAVE_SYS_SYSINFO_H) && defined(HAVE_SYSINFO)
  {
    struct sysinfo info;
    sysinfo(&info);
    ram = info.totalram / 1024;
    if (ram > 1)
      sprintf(str_ram, "%lu", (unsigned long)ram);
  }
#endif

  /* Repertoire courant */

  str_directory[0] = '\0';

#if defined(HAVE_GETCWD)
  if (getcwd(str_directory, PATH_MAX) == NULL)
    str_directory[0] = '\0';
#endif

  /* Affichage de la configuration du cas de lancement */
  /*---------------------------------------------------*/

  printf("\n\n%s\n", _("Case configuration\n"
                       "------------------\n"));

  if (strlen(str_date) > 0) {
    printf("  ");
    ecs_print_padded_str(_("Date"), 19);
    printf(" : %s\n", str_date);
  }

  if (strlen(str_system) > 0) {
    printf("  ");
    ecs_print_padded_str(_("System"), 19);
    printf(" : %s\n", str_system);
  }

  if (strlen(str_machine) > 0) {
    printf("  ");
    ecs_print_padded_str(_("Machine"), 19);
    printf(" : %s\n", str_machine);
  }

  printf("  ");
  ecs_print_padded_str(_("Processor"), 19);
  printf(" : %s\n", _sys_info_cpu());

  if (ram > 0) {
    printf("  ");
    ecs_print_padded_str(_("Memory"), 19);
    printf(" : %s\n", str_ram);
  }

  if (strlen(str_user) > 0) {
    printf("  ");
    ecs_print_padded_str(_("User"), 19);
    printf(" : %s\n", str_user);
  }

  if (strlen(str_directory) > 0) {
    printf("  ");
    ecs_print_padded_str(_("Directory"), 19);
    printf(" : %s\n", str_directory);
  }

  printf("\n");

  printf("  ");
  ecs_print_padded_str(_("Case name"), 19);
  printf(" : %s\n", cmd->nom_cas);

  printf("  ");
  ecs_print_padded_str(_("Mesh file"), 19);
  printf(" : %s\n", cmd->fic_maillage);

  printf("\n");
}

/*----------------------------------------------------------------------------
 *  Fonction qui diagnostique l'erreur de non accessibilite
 *   du fichier dont le nom est donne
 *----------------------------------------------------------------------------*/

static void
ecs_loc_cmd__diagnost_non_acces(const char  *fic_name,
                                const char  *msg_cmd_err)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if defined(ENOENT)
  if (errno == ENOENT)
    ecs_error(__FILE__, __LINE__, 0,
              msg_cmd_err, fic_name, _("The file does not exist."));
#endif

#if defined(EACCES)
  if (errno == EACCES)
    ecs_error(__FILE__, __LINE__, 0,
              msg_cmd_err, fic_name,
              _("Read permission on this file is refused."));
#endif

#if defined(ENAMETOOLONG)
  if (errno == ENAMETOOLONG)
    ecs_error(__FILE__, __LINE__, 0,
              msg_cmd_err, fic_name,
              _("The filename is too long (system limit)."));
#endif

  ecs_error(__FILE__, __LINE__, errno,
            msg_cmd_err, fic_name, "");
}

/*----------------------------------------------------------------------------
 *  Fonction qui verifie que le nom du fichier donne
 *   correspond bien a un fichier physiquement existant
 *----------------------------------------------------------------------------*/

static void
ecs_loc_cmd__teste_exist_fic(const char  *fic_name,
                             const char  *msg_cmd_err)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if defined(HAVE_STAT)

  struct stat buf;

  if (stat(fic_name, &buf) != 0) {

    ecs_loc_cmd__diagnost_non_acces(fic_name,
                                    msg_cmd_err);

  }
  else {

    if (S_ISREG(buf.st_mode) != true) {

      ecs_error(__FILE__, __LINE__, 0,
                msg_cmd_err,
                fic_name, _("The file is not a regular file."));

    }
  }

#else /* HAVE_STAT */

  FILE * fic;

  if ((fic = fopen(fic_name, "r")) == NULL) {

    ecs_loc_cmd__diagnost_non_acces(fic_name,
                                    msg_cmd_err);

  }
  else {

    fclose(fic);
  }

#endif /* HAVE_STAT */
}

/*============================================================================
 *                             Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Fonction qui lit la ligne de commande
 *----------------------------------------------------------------------------*/

ecs_cmd_t *
ecs_cmd__lit_arg(int    argc,
                 char  *argv[])
{
  int        iarg;
  size_t     lng;
  bool       bool_cmd_option_case = false;
  bool       bool_cmd_option_help = false;
  bool       bool_cmd_option_output_file = false;
  bool       bool_cmd_option_version = false;
  bool       bool_no_write = false;
  bool       bool_num = false;

  ecs_cmd_t  *cmd;

  const char  *cle_fmt = NULL;
  const char arg_err_keyword[]
    = N_("Error in command line specification.\n\n"
         "Keyword \"%s\" of option \"%s\" is not recognized.\n");

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /*----------------------------------------*/
  /* Initialisation des options de commande */
  /*----------------------------------------*/

  cmd = _cmd_initialize();

  /*---------------------------------------------*/
  /* Lecture des options de la ligne de commande */
  /*---------------------------------------------*/

  for (iarg = 1; iarg < argc; iarg++) {

    if (!strcmp ("--case", argv[iarg])) {

      if (bool_cmd_option_case == false)
        bool_cmd_option_case = true;
      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      if (argc - 1 > iarg && strncmp(argv[iarg + 1], "-", 1)) {

        lng = strlen(argv[++iarg]) + 1;
        ECS_REALLOC(cmd->nom_cas, lng, char);
        strcpy(cmd->nom_cas, argv[iarg]);

      }
      else {
        ecs_loc_cmd__aff_manque_arg(argc, argv, argv[iarg]);
      }
    }
    else if (!strcmp ("--dump", argv[iarg])) {

      /* Option non documentee */

      cmd->nbr_dump = 1;

      /* 1 argument optionnel : nombre n d'elements affiches en echo */

      if (argc -1 > iarg && strncmp(argv[iarg + 1], "-", 1))
        cmd->nbr_dump = (ecs_int_t) atoi(argv[++iarg]);

    }
    else if (!strcmp ("--no-write", argv[iarg]))
      bool_no_write = true;

    else if (!strcmp("-h", argv[iarg]) ||
             !strcmp("--help", argv[iarg]))
      bool_cmd_option_help = true;

    else if (!strcmp ("--post-error", argv[iarg]))
      _read_post_opt(argc, argv, cmd->post_err, &iarg);

    else if (!strcmp ("--post-volume", argv[iarg]))
      _read_post_opt(argc, argv, cmd->post_vol, &iarg);

    else if (!strcmp ("--log", argv[iarg])) {

      int have_error = 0;
      const char *outfic;
      char *outfic_err;
      FILE *ptrfic;

      if (bool_cmd_option_output_file == false)
        bool_cmd_option_output_file = true;
      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      if (argc - 1 > iarg && strncmp(argv[iarg + 1], "-", 1))
        outfic = argv[++iarg];
      else
        outfic = "preprocessor.log";

      ECS_MALLOC(outfic_err, strlen(outfic) + strlen(".err") + 1, char);

      sprintf(outfic_err, "%s.err", outfic);

#if defined(HAVE_DUP2)
      if ((ptrfic = freopen(outfic, "w", stdout)) == NULL ||
          dup2(fileno(ptrfic), fileno(stderr)) == -1)
        have_error = 1;
#else
      if (freopen(outfic, "w", stdout) == NULL ||
          freopen(outfic_err, "w", stderr) == NULL)
        have_error = 1;
#endif
      if (have_error) {
        _print_preamble(argc, argv);
        ecs_error(__FILE__, __LINE__, errno,
                  _("It is impossible to redirect the standard output "
                    "to file:\n%s"), outfic);
      }

      ECS_FREE(outfic_err);
    }

    else if (   !strcmp ("-o", argv[iarg])
             || !strcmp ("--out", argv[iarg])) {

      const char *outfic;

      if (cmd->nom_out != NULL)
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      if (argc - 1 > iarg && strncmp(argv[iarg + 1], "-", 1))
        outfic = argv[++iarg];
      else
        outfic = "mesh_input";

      ECS_MALLOC(cmd->nom_out, strlen(outfic) + 1, char);

      strcpy(cmd->nom_out, outfic);

    }

    else if (!strcmp ("--reorient", argv[iarg]))
      cmd->correct_orient = true;

    else if (!strcmp ("--version", argv[iarg]))
      bool_cmd_option_version = true;

    /* Option de choix de format */

    else if (!strcmp ("--format", argv[iarg])) {

      if (cle_fmt == NULL)
        cle_fmt = argv[iarg + 1];
      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      iarg++; /* Un argument lu -> avancer iarg */

    }

    /* Numéros de maillage */

    else if (!strcmp ("--num", argv[iarg])) {

      int iarg_num;

      cmd->n_num_maillage = 0;
      cmd->num_maillage = NULL;

      if (bool_num == false)
        bool_num = true;
      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      for (iarg_num  = iarg+1; iarg_num < argc; iarg_num++) {
        if (! (isdigit(argv[iarg_num][0]) && atoi(argv[iarg + 1]) > 0))
          break;
      }

      cmd->n_num_maillage = iarg_num - (iarg+1);

      if (cmd->n_num_maillage == 0) {
        _print_preamble(argc, argv);
        _print_help();
        ecs_error(__FILE__, __LINE__, 0,
                  _(arg_err_keyword), argv[iarg + 1], argv[iarg]);
      }

      ECS_MALLOC(cmd->num_maillage, cmd->n_num_maillage, int);
      for (iarg_num  = 0; iarg_num < cmd->n_num_maillage; iarg_num++)
        cmd->num_maillage[iarg_num] = atoi(argv[++iarg]);

    }

    /* Sous-options de création de groupes de cellules */

    else if (!strcmp ("--grp-cel", argv[iarg])) {

      if (   iarg + 1 < argc
          && !strcmp("section", argv[iarg + 1]))
        cmd->grp_cel_section = true;

      else if (   iarg + 1 < argc
          && !strcmp("zone", argv[iarg + 1]))
        cmd->grp_cel_zone = true;

      else {
        _print_preamble(argc, argv);
        _print_help();
        ecs_error(__FILE__, __LINE__, 0,
                  _(arg_err_keyword), argv[iarg + 1], argv[iarg]);
      }

      iarg++; /* Un sous-argument lu -> avancer iarg */
    }

    /* Création de groupes de faces */

    else if (!strcmp ("--grp-fac", argv[iarg])) {

      if (iarg + 1 < argc && !strcmp("section", argv[iarg + 1]))
        cmd->grp_fac_section = true;

      else if (iarg + 1 < argc && !strcmp("zone", argv[iarg + 1]))
        cmd->grp_fac_zone = true;

      else {
        _print_preamble(argc, argv);
        _print_help();
        ecs_error(__FILE__, __LINE__, 0,
                  _(arg_err_keyword), argv[iarg + 1], argv[iarg]);
      }

      iarg++; /* Un sous-argument lu -> avancer iarg */
    }

    /* Nom de fichier */

    else if (argv[iarg][0] != '-') {
      ECS_MALLOC(cmd->fic_maillage, strlen(argv[iarg]) + 1, char);
      strcpy(cmd->fic_maillage, argv[iarg]);
    }

    else {

      _print_preamble(argc, argv);

      _print_help();

      ecs_error(__FILE__, __LINE__, 0,
                _("Option \"%s\" is not recognized.\n"), argv[iarg]);

    }

  } /* Fin : boucle sur les arguments */

  /*---------------------------------------------------------------*/
  /* Affichage de la ligne de commande, du titre et de la version */
  /*---------------------------------------------------------------*/

  _print_preamble(argc, argv);

  if (cmd->nom_out == NULL && bool_no_write == false) {
    ECS_MALLOC(cmd->nom_out,
               strlen("mesh_input") + 1,
               char);
    strcpy(cmd->nom_out, "mesh_input");
  }

  /*-----------------------------------------*/
  /* Options devant etre traitees en premier */
  /*-----------------------------------------*/

  if (bool_cmd_option_help == true)
    _print_help();

  /*-------------------------------------------*/
  /* Options provoquant l'arret de l'execution */
  /*-------------------------------------------*/

  if (   (argc <= 1)
      || (             bool_cmd_option_help    == true)
      || (argc == 2 && bool_cmd_option_version == true)
      || (argc == 3 && bool_cmd_option_help    == true
                    && bool_cmd_option_version == true)) {

    ecs_cmd__detruit(cmd);
    ecs_exit(EXIT_SUCCESS);
  }

  /*-----------------------------------------------------------*/
  /* Verification que les donnees necessaires ont ete fournies */
  /*-----------------------------------------------------------*/

  if (cmd->n_num_maillage == 0) {
    cmd->n_num_maillage = 1;
    ECS_MALLOC(cmd->num_maillage, 1, int);
    cmd->num_maillage[0] = 0;
  }

  cmd->fmt_maillage = ecs_pre__type_format(cmd->fic_maillage, cle_fmt);

  if (cmd->fic_maillage == NULL) {

    _print_help();

    ecs_error(__FILE__, __LINE__, 0,
              _("Error in command line specification:\n\n"
                "Missing input mesh file name."));
  }

  /*--------------------------*/
  /* Verification des donnees */
  /*--------------------------*/

  /* Verification que le fichier de maillage existe */
  /*------------------------------------------------*/

  ecs_loc_cmd__teste_exist_fic(cmd->fic_maillage,
                               _("Mesh file \"%s\" "
                                 "is not accessible.\n%s"));

  /*---------------------------------------------------*/
  /* Affichage de la configuration du cas de lancement */
  /*---------------------------------------------------*/

  ecs_loc_cmd__aff_config(cmd);

  return cmd;
}

/*----------------------------------------------------------------------------
 *  Fonction liberant une structure `ecs_cmd_t' donnee en argument
 *  Elle renvoie un pointeur NULL
 *----------------------------------------------------------------------------*/

ecs_cmd_t *
ecs_cmd__detruit(ecs_cmd_t  *cmd)
{
  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  assert(cmd != NULL);

  /* Liberation du contenu de la structure `ecs_cmd_t' */
  /*===============================================*/

  if (cmd->nom_cas != NULL)
    ECS_FREE(cmd->nom_cas);

  if (cmd->nom_out != NULL)
    ECS_FREE(cmd->nom_out);

  ECS_FREE(cmd->fic_maillage);

  if (cmd->n_num_maillage > 0)
    ECS_FREE(cmd->num_maillage);

  /* Liberation de la structure `ecs_cmd_t' */
  /*====================================*/

  ECS_FREE(cmd);

  return cmd;
}

/*----------------------------------------------------------------------------*/
