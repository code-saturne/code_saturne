<!--
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
-->

\page cs_ug_user_settings Setting up your environment

Shortcuts and command completion  {#sec_prg_environment_cs}
================================

It is recommended before running code_saturne to define an alias to the
`code_saturne` script (see [sec_prg_environement]), for example:

```
alias cs='${install_prefix}/bin/code_saturne'
```

where '${install_prefix} is the base directory where code_saturne and its components
have been installed.

This setting can be defined in the user's `.bashrc`,`.bash_alias`, or
equivalent files depending on the environment.
This step may be skipped if '${install_prefix} is in a standard location (such as
`/usr` or `/usr/local`), or if the code is already available as an environment module.

When using the *bash* shell, a completion file may be sourced so as to
allow for syntax auto-completion:
```
source ${install_prefix}/etc/bash_completion.d/code_saturne'.
```
This can alse  be defined in the user's `.bashrc` file.

When using multiple versions of the code, different aliases should be used for
each version.

Independently of these settings, using the absolute path to the `code_saturne`
command is always possible.

For more advanced settings, a [configuration file](@ref cs_user_configuration_file)
may be used.

Configuration file {#cs_user_configuration_file}
==================

A configuration file for code_saturne is available in `${install_prefix}/etc`.
This file can be useful as a post-install step for computing environments using a
batch system, for separate front-end and compute systems (such as some Cray systems),
or for coupling with Syrthes (see the installation documentation for more details).

A user may define a local configuration, by copying
`${install_prefix}/etc/code_saturne.cfg` (if present)
or `${install_prefix}/etc/code_saturne.cfg.template` to
`$HOME/.code_saturne.cfg`, then uncomment and define the applicable sections.

Note that this user configuration file's settings usually apply to all installed
code_saturne versions, so only the necessary options should be defined.

Two options in the `.code_saturne.cfg` file can be useful for the user:
* Set the optional temporary directory
  (see [case_structure_scratchdir](@ref case_structure_scratchdir))
  for more details on the temporary execution directory).
* Set the mesh database directory: it is possible to indicate a path where
  meshes are stored. In this case, the GUI will propose this directory
  automatically for mesh selection. Without the GUI, it is
  then possible to fill in the `cs_user_scripts.py` file (see
  [sec_prg_stepbystepcalculation) with the name of the desired mesh of the
  database directory and the code will find it automatically (be careful if you
  have the same name for a mesh in the database directory
  and in the `MESH` directory: the mesh in `MESH` will be used).
