/*============================================================================
 * Command line option parser and tracking.
 *============================================================================*/

/*
  This file is part of the Code_Saturne Preprocessor, element of the
  Code_Saturne CFD tool.

  Copyright (C) 1999-2010 EDF S.A., France

  contact: saturne-support@edf.fr

  The Code_Saturne Preprocessor is free software; you can redistribute it
  and/or modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The Code_Saturne Preprocessor is distributed in the hope that it will be
  useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with the Code_Saturne Preprocessor; if not, write to the
  Free Software Foundation, Inc.,
  51 Franklin St, Fifth Floor,
  Boston, MA  02110-1301  USA
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
#include "ecs_cmd_priv.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Definitions of file names and extensions
 *----------------------------------------------------------------------------*/

#define ECS_CMD_EXEC_NAME                             "cs_preprocess"

#define ECS_CMD_LOGFILE_NAME_DEFAULT               "preprocessor.log"
#define ECS_CMD_OUTFILE_NAME_DEFAULT            "preprocessor_output"

#define ECS_CMD_POST_CASE_DEFAULT                        "preprocess"

/*----------------------------------------------------------------------------
 * Définition des mots-cles pour les options de la ligne de commande
 *----------------------------------------------------------------------------*/

#define ECS_CMD_KEY_MESH_GRP_SECTION                        "section"
#define ECS_CMD_KEY_MESH_GRP_ZONE                              "zone"

/*----------------------------------------------------------------------------
 *  Définition des options de la ligne de commande
 *----------------------------------------------------------------------------*/

#define ECS_CMD_OPTION_CASE                                  "--case"

#define ECS_CMD_OPTION_DUMP                                  "--dump"
#define ECS_CMD_OPTION_NULL_COMM                         "--no-write"
#define ECS_CMD_OPTION_FMT_MESH_FILE                       "--format"
#define ECS_CMD_OPTION_NUM_MESH                               "--num"
#define ECS_CMD_OPTION_GRP_CEL_MESH                       "--grp-cel"
#define ECS_CMD_OPTION_GRP_FAC_MESH                       "--grp-fac"

#define ECS_CMD_OPTION_HELP                                  "--help"
#define ECS_CMD_OPTION_HELP_1                                    "-h"

#define ECS_CMD_OPTION_LOG_FILE                               "--log"

#define ECS_CMD_OPTION_OUTPUT_FILE                            "--out"
#define ECS_CMD_OPTION_OUTPUT_FILE_1                             "-o"

#define ECS_CMD_OPTION_ORIENT_CORREC                     "--reorient"

#if defined(HAVE_CGNS)
#define ECS_CMD_OPTION_POST_CGNS                             "--cgns"
#endif /* HAVE_CGNS */
#define ECS_CMD_OPTION_POST_ENS                           "--ensight"
#if defined(HAVE_MED)
#define ECS_CMD_OPTION_POST_MED                               "--med"
#endif /* HAVE_MED */

#define ECS_CMD_OPTION_POST_MAIN                           "--volume"
#define ECS_CMD_OPTION_POST_INFO                             "--info"

#define ECS_CMD_OPTION_VERSION                            "--version"

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
ecs_loc_cmd__aff_version(void)
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
ecs_loc_cmd__aff_titre(int    argc,
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

  ecs_loc_cmd__aff_version();
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
ecs_loc_cmd__aff_aide(void)
{
  char opt_str[81];

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf(_("\n\nUsage:  %s [<options>] <file>\n"),
         ECS_CMD_EXEC_NAME);

  /* General options */
  /*-----------------*/

  printf(_("\n\nGeneral options:\n\n"));

  _fct_prt(ECS_CMD_OPTION_HELP_1, "", _(": this help message"));
  _fct_prt(ECS_CMD_OPTION_HELP, "", _(": same"));

  printf("\n");

  _fct_prt(ECS_CMD_OPTION_LOG_FILE, _("[file]"),
           _(": redirect terminal output to a file"));

  sprintf(opt_str, _("  (default file: \"%s\")"),
          ECS_CMD_LOGFILE_NAME_DEFAULT);
  _fct_prt("", "", opt_str);

  printf("\n");

  _fct_prt(ECS_CMD_OPTION_NULL_COMM, "",
           _(": do not write preprocessor output"));

  printf("\n");

  _fct_prt(ECS_CMD_OPTION_OUTPUT_FILE_1, _("<file>"),
           _(": output file name"));
  sprintf(opt_str, _("  (default file: \"%s\")"),
          ECS_CMD_OUTFILE_NAME_DEFAULT);
  _fct_prt("", "", opt_str);

  _fct_prt(ECS_CMD_OPTION_OUTPUT_FILE, _("<file>"), _(": same"));

  printf("\n");

  _fct_prt(ECS_CMD_OPTION_ORIENT_CORREC, "",
           _(": if necessary, correct orientation of"));
  _fct_prt("", "", _("  cells and faces"));

  printf("\n");

  _fct_prt(ECS_CMD_OPTION_VERSION, "",
           _(": print version number"));

  /* Post-processing options */
  /*-------------------------*/

  printf(_("\n\nPostprocessing options:\n\n"));

  _fct_prt(ECS_CMD_OPTION_CASE,
           _("<name>"),
           _(": case name (without this option,"));

  sprintf(opt_str,
          _("  the default name is: \"%s\""),
          ECS_CMD_POST_CASE_DEFAULT);
  _fct_prt("", "", opt_str);

  printf("\n");

#if defined(HAVE_CGNS)

  sprintf(opt_str,
          _(": %s geometry output"), "CGNS");
  _fct_prt(ECS_CMD_OPTION_POST_CGNS, _("[<sub-options>]"), opt_str);

  printf("\n");

#endif /* HAVE_CGNS */

  sprintf(opt_str,
          _(": %s geometry output"), "EnSight Gold");
  _fct_prt(ECS_CMD_OPTION_POST_ENS, _("[<sub-options>]"), opt_str);


#if defined(HAVE_MED)

  printf("\n");

  sprintf(opt_str,
          _(": %s geometry output"), "MED");
  _fct_prt(ECS_CMD_OPTION_POST_MED, _("[<sub-options>]"), opt_str);

#endif /* HAVE_MED */

  /* Mesh selection sub-options */
  /*----------------------------*/

  printf(_("\n\nMesh selection options:\n\n"));

  _fct_prt(ECS_CMD_OPTION_FMT_MESH_FILE, _("<keyword>"),
           _(": selection of mesh file format"));

  _fct_prt(ECS_CMD_OPTION_NUM_MESH, "<n> [...]",
           _(": selection of mesh numbers in file"));
  _fct_prt("", "",
           _("  (if the format allows it)"));

  _fct_prt(ECS_CMD_OPTION_GRP_CEL_MESH, _("<keyword>"),
           _(": add groups of cells"));
  _fct_prt("", "", _("   * based on sections: keyword \""
                     ECS_CMD_KEY_MESH_GRP_SECTION"\""));
  _fct_prt("", "", _("   * based on zones:    keyword \""
                     ECS_CMD_KEY_MESH_GRP_ZONE"\""));
  _fct_prt("", "", _("  (based on format features/conventions)"));

  _fct_prt(ECS_CMD_OPTION_GRP_FAC_MESH, _("<keyword>"),
           _(": add groups of faces"));
  _fct_prt("", "", _("   * based on sections: keyword \""
                     ECS_CMD_KEY_MESH_GRP_SECTION"\""));
  _fct_prt("", "", _("   * based on zones:    keyword \""
                     ECS_CMD_KEY_MESH_GRP_ZONE"\""));
  _fct_prt("", "", _("  (based on format features/conventions)"));

  printf(_("\n\nAvailable mesh formats:\n"));
  printf(_("                                  extension:    keyword:\n"));

  ecs_pre__aff_formats();

  /* Post-processign sub-options */
  /*-----------------------------*/

  opt_str[0] = '\0';

#if defined(HAVE_CGNS)
  strcat(opt_str, ECS_CMD_OPTION_POST_CGNS);
  strcat(opt_str, ", ");
#endif

  strcat(opt_str, ECS_CMD_OPTION_POST_ENS);

#if defined(HAVE_MED)
  strcat(opt_str, ", ");
  strcat(opt_str, ECS_CMD_OPTION_POST_MED);
#endif

  printf(_("\n\nPostprocessing selection sub-options\n"
           " (%s):\n\n"), opt_str);

  _fct_prt(ECS_CMD_OPTION_POST_MAIN,
           "", _(": activate output of main mesh"));
  sprintf(opt_str, _("  (by default if %s not given)"),
          ECS_CMD_OPTION_POST_INFO);
  _fct_prt("", "", opt_str);

  _fct_prt(ECS_CMD_OPTION_POST_INFO,
              "", _(": activate output of information meshes"));
  sprintf(opt_str, _("  (by default if %s not given)"),
          ECS_CMD_OPTION_POST_MAIN);
  _fct_prt("", "", opt_str);

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

  ecs_loc_cmd__aff_titre(argc, argv);

  ecs_loc_cmd__aff_aide();

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

  ecs_loc_cmd__aff_titre(argc, argv);

  ecs_loc_cmd__aff_aide();

  ecs_error(__FILE__, __LINE__, 0,
            _("Error in command line specification.\n\n"
              "Option \"%s\" requires an argument."), opt);
}

/*----------------------------------------------------------------------------
 *  Fonction qui lit les sous-options d'un post traitement
 *----------------------------------------------------------------------------*/

static ecs_cmd_post_t  *
ecs_loc_cmd__lit_arg_post(int    argc,
                          char  *argv[],
                          int   *argpos)
{
  int         iarg;
  int         iarg_prec;

  bool        bool_fin = false;

  ecs_cmd_post_t  *cmd_post;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  iarg_prec = *argpos;

  cmd_post = NULL;

  ECS_MALLOC(cmd_post, 1, ecs_cmd_post_t);

  cmd_post->volume = false;
  cmd_post->info   = false;

  for (iarg = *argpos + 1; iarg < argc && bool_fin == false; iarg++) {

    if (!strcmp (ECS_CMD_OPTION_POST_MAIN, argv[iarg]))
      cmd_post->volume = true;

    else if (!strcmp (ECS_CMD_OPTION_POST_INFO, argv[iarg]))
      cmd_post->info = true;

    else {

      /* Autre option (pas une sous-option) -> on a fini */

      iarg--;
      bool_fin = true;

    }
  }

  /* Si aucune option de filtrage du post traitement n'est
     activée, tous les types post traitements sont actifs */

  if (cmd_post->volume == false && cmd_post->info   == false) {
    cmd_post->volume = true;
    cmd_post->info   = true;
  }

  *argpos = iarg - 1;

  return cmd_post;
}

/*----------------------------------------------------------------------------
 *  Fonction qui initialise les options de commande
 *----------------------------------------------------------------------------*/

static ecs_cmd_t *
ecs_loc_cmd__initialise(void)
{
  size_t   lng;
  ecs_cmd_t  *cmd;

  /*xxxxxxxxxxxxxxxxxxxxxxxxxxx Instructions xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ECS_MALLOC(cmd, 1, ecs_cmd_t);

  /* Initialisations */
  /*=================*/

  cmd->fic_maillage                  = NULL;

  lng =   strlen(ECS_CMD_POST_CASE_DEFAULT) + 1;
  ECS_MALLOC(cmd->nom_cas, lng, char);
  strcpy(cmd->nom_cas, ECS_CMD_POST_CASE_DEFAULT);

  cmd->nom_out                       = NULL;

  cmd->nbr_dump = 0;

  cmd->n_num_maillage = 0;
  cmd->num_maillage = NULL;
  cmd->fmt_maillage = ECS_PRE_FORMAT_NUL;
  cmd->grp_cel_section = false;
  cmd->grp_cel_zone    = false;
  cmd->grp_fac_section = false;
  cmd->grp_fac_zone    = false;

#if defined(HAVE_CGNS)

  cmd->post_cgns = NULL;

#endif /* HAVE_CGNS */

  cmd->post_ens = NULL;

#if defined(HAVE_MED)

  cmd->post_med = NULL;

#endif /* HAVE_MED */

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

  cmd = ecs_loc_cmd__initialise();

  /*---------------------------------------------*/
  /* Lecture des options de la ligne de commande */
  /*---------------------------------------------*/

  for (iarg = 1; iarg < argc; iarg++) {

    if (!strcmp (ECS_CMD_OPTION_CASE, argv[iarg])) {

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
    else if (!strcmp (ECS_CMD_OPTION_DUMP, argv[iarg])) {

      /* Option non documentee */

      cmd->nbr_dump = 1;

      /* 1 argument optionnel : nombre n d'elements affiches en echo */

      if (argc -1 > iarg && strncmp(argv[iarg + 1], "-", 1))
        cmd->nbr_dump = (ecs_int_t) atoi(argv[++iarg]);

    }
    else if (!strcmp (ECS_CMD_OPTION_NULL_COMM, argv[iarg]))
      bool_no_write = true;

    else if (!strcmp(ECS_CMD_OPTION_HELP_1, argv[iarg]) ||
             !strcmp(ECS_CMD_OPTION_HELP, argv[iarg]))
      bool_cmd_option_help = true;

#if defined(HAVE_CGNS)

    else if (!strcmp (ECS_CMD_OPTION_POST_CGNS, argv[iarg])) {

      if (cmd->post_cgns == NULL)
        cmd->post_cgns = ecs_loc_cmd__lit_arg_post(argc,
                                                   argv,
                                                   &iarg);
      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);
    }

#endif /* HAVE_CGNS */

    else if (!strcmp (ECS_CMD_OPTION_POST_ENS, argv[iarg])) {

      if (cmd->post_ens == NULL)
        cmd->post_ens = ecs_loc_cmd__lit_arg_post(argc,
                                                  argv,
                                                  &iarg);

      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);
    }

#if defined(HAVE_MED)

    else if (!strcmp (ECS_CMD_OPTION_POST_MED, argv[iarg])) {

      if (cmd->post_med == NULL)
        cmd->post_med = ecs_loc_cmd__lit_arg_post(argc,
                                                  argv,
                                                  &iarg);

      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);
    }

#endif /* HAVE_MED */

    else if (!strcmp (ECS_CMD_OPTION_LOG_FILE, argv[iarg])) {

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
        outfic = ECS_CMD_LOGFILE_NAME_DEFAULT;

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
        ecs_loc_cmd__aff_titre(argc, argv);
        ecs_error(__FILE__, __LINE__, errno,
                  _("It is impossible to redirect the standard output "
                    "to file:\n%s"), outfic);
      }

      ECS_FREE(outfic_err);
    }

    else if (   !strcmp (ECS_CMD_OPTION_OUTPUT_FILE_1, argv[iarg])
             || !strcmp (ECS_CMD_OPTION_OUTPUT_FILE, argv[iarg])) {

      const char *outfic;

      if (cmd->nom_out != NULL)
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      if (argc - 1 > iarg && strncmp(argv[iarg + 1], "-", 1))
        outfic = argv[++iarg];
      else
        outfic = ECS_CMD_OUTFILE_NAME_DEFAULT;

      ECS_MALLOC(cmd->nom_out, strlen(outfic) + 1, char);

      strcpy(cmd->nom_out, outfic);

    }

    else if (!strcmp (ECS_CMD_OPTION_ORIENT_CORREC, argv[iarg]))
      cmd->correct_orient = true;

    else if (!strcmp (ECS_CMD_OPTION_VERSION, argv[iarg]))
      bool_cmd_option_version = true;

    /* Option de choix de format */

    else if (!strcmp (ECS_CMD_OPTION_FMT_MESH_FILE, argv[iarg])) {

      if (cle_fmt == NULL)
        cle_fmt = argv[iarg + 1];
      else
        ecs_loc_cmd__aff_opt_en_double(argc, argv, argv[iarg]);

      iarg++; /* Un argument lu -> avancer iarg */

    }

    /* Numéros de maillage */

    else if (!strcmp (ECS_CMD_OPTION_NUM_MESH, argv[iarg])) {

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
        ecs_loc_cmd__aff_titre(argc, argv);
        ecs_loc_cmd__aff_aide();
        ecs_error(__FILE__, __LINE__, 0,
                  _(arg_err_keyword), argv[iarg + 1], argv[iarg]);
      }

      ECS_MALLOC(cmd->num_maillage, cmd->n_num_maillage, int);
      for (iarg_num  = 0; iarg_num < cmd->n_num_maillage; iarg_num++)
        cmd->num_maillage[iarg_num] = atoi(argv[++iarg]);

    }

    /* Sous-options de création de groupes de cellules */

    else if (!strcmp (ECS_CMD_OPTION_GRP_CEL_MESH, argv[iarg])) {

      if (   iarg + 1 < argc
          && !strcmp(ECS_CMD_KEY_MESH_GRP_SECTION, argv[iarg + 1]))
        cmd->grp_cel_section = true;

      else if (   iarg + 1 < argc
          && !strcmp(ECS_CMD_KEY_MESH_GRP_ZONE, argv[iarg + 1]))
        cmd->grp_cel_zone = true;

      else {
        ecs_loc_cmd__aff_titre(argc, argv);
        ecs_loc_cmd__aff_aide();
        ecs_error(__FILE__, __LINE__, 0,
                  _(arg_err_keyword), argv[iarg + 1], argv[iarg]);
      }

      iarg++; /* Un sous-argument lu -> avancer iarg */
    }

    /* Création de groupes de faces */

    else if (!strcmp (ECS_CMD_OPTION_GRP_FAC_MESH, argv[iarg])) {

      if (   iarg + 1 < argc
          && !strcmp(ECS_CMD_KEY_MESH_GRP_SECTION, argv[iarg + 1]))
        cmd->grp_fac_section = true;

      else if (   iarg + 1 < argc
          && !strcmp(ECS_CMD_KEY_MESH_GRP_ZONE, argv[iarg + 1]))
        cmd->grp_fac_zone = true;

      else {
        ecs_loc_cmd__aff_titre(argc, argv);
        ecs_loc_cmd__aff_aide();
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

      ecs_loc_cmd__aff_titre(argc, argv);

      ecs_loc_cmd__aff_aide();

      ecs_error(__FILE__, __LINE__, 0,
                _("Option \"%s\" is not recognized.\n"), argv[iarg]);

    }

  } /* Fin : boucle sur les arguments */

  /*---------------------------------------------------------------*/
  /* Affichage de la ligne de commande, du titre et de la version */
  /*---------------------------------------------------------------*/

  ecs_loc_cmd__aff_titre(argc, argv);

  if (cmd->nom_out == NULL && bool_no_write == false) {
    ECS_MALLOC(cmd->nom_out,
               strlen(ECS_CMD_OUTFILE_NAME_DEFAULT) + 1,
               char);
    strcpy(cmd->nom_out, ECS_CMD_OUTFILE_NAME_DEFAULT);
  }

  /*-----------------------------------------*/
  /* Options devant etre traitees en premier */
  /*-----------------------------------------*/

  if (bool_cmd_option_help == true)
    ecs_loc_cmd__aff_aide();

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

    ecs_loc_cmd__aff_aide();

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

#if defined(HAVE_CGNS)
  if (cmd->post_cgns != NULL)
    ECS_FREE(cmd->post_cgns);
#endif /* HAVE_CGNS */

  if (cmd->post_ens != NULL)
    ECS_FREE(cmd->post_ens);

#if defined(HAVE_MED)

  if (cmd->post_med != NULL)
    ECS_FREE(cmd->post_med);

#endif /* HAVE_MED */

  /* Liberation de la structure `ecs_cmd_t' */
  /*====================================*/

  ECS_FREE(cmd);

  return cmd;
}

/*----------------------------------------------------------------------------*/
