<!--
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

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

\page cs_ug_studymanager The STUDYMANAGER tool

[TOC]

STUDYMANAGER
============

This document presents the STUDYMANAGER command. The aim of this command
is to drive code_saturne's cases automatically, to compare checkpoint files and
to display results.

STUDYMANAGER is a small framework to automate the launch of code_saturne
computations and do some operations on new results.

The script needs a source directory of code_saturne cases, called the
**repository**, which will be duplicated into a **destination**
directory, from which computations will be run.

For each duplicated case, STUDYMANAGER can compile the user
files, run the case, compare the obtained checkpoint file with the
previous one from the repository, and plot curves in order to
illustrate the computations.

For all these steps, STUDYMANAGER generates two reports:
- A _global_ report which summarizes the status of each case.
- A _detailed_ report which gives the differences between the new results
  and the previous ones in the repository, and displays
 the defined plots.

In the **repository**, previous results of computations are required only
for checkpoint files comparison purposes. They can be also useful, if the
user needs to run specific scripts.

Installation and prerequisites
------------------------------

STUDYMANAGER is available as a code_saturne command, and
does not need a specific installation: the related files
are installed with the other Python scripts of code_saturne. Nevertheless,
additional prerequisites which may be required are:
- `numpy`,
- `matplotlib`

Since these are used in a dynamic manner, they may be added after the code_saturne
installation, and do not require any re-installation.

Command line options
====================

The command line options can be listed with the following command:
`code_saturne studymanager -h`

- `-h, --help`: show the help message and exit
- `-f FILE, --file=FILE`: give the parameters file for STUDYMANAGER.
   This file is mandatory, and therefore this option must be completed
- `-q, --quiet`: do not print status messages to `stdout`
- `-u, --update`: update installation paths in scripts (i.e.
  `code_saturne`) only in the **repository**, reinitialize XML files of
  parameters and compile sources
- `-x, --update-xml`: update only XML files in the repository
- `-t, --test-compile`: compile all cases
- `-r, --run`: run all cases
- `-n N_ITERATIONS, --n-iterations=N_ITERATIONS`: maximum number of
  iterations for cases of the study
- `-c, --compare`: compare checkpoint files between **repository** and
  **destination**
- `-d REFERENCE, --ref-dir=REFERENCE`: absolute reference directory to
  compare dest with
- `-p, --post`: postprocess results of computations
- `-m ADDRESS1 ADDRESS2 ..., --mail=ADDRESS1 ADDRESS2 ...`: addresses for
  sending the reports
- `-l LOG_FILE, --log=LOG_FILE`: name of STUDYMANAGER log file (default
  value is `studymanager.log`)
- `-z, --disable-tex`: disable text rendering with \f$ \mbox{\LaTeX} \f$ when plotting
  with Matplotlib. It then uses Mathtext instead which is less complete
- `--rm`: remove all existing run directories in destination directory
- `--fow`: overwrite files in MESH and POST directories in destination
  directory
- `-s, --skip-pdflatex`: disable tex reports compilation with pdflatex
- `--fmt=DEFAULT_FMT`: set the global format for exporting Matplotlib
  figure (default is pdf)
- `--repo=REPO_PATH`: force the path to the repository directory
- `--dest=DEST_PATH`: force the path to the destination directory
- `-g, --debug`: activate debugging mode
- `--with-tags=WITH_TAGS`: only process runs with all specified tags
  (separated by commas)
- `--without-tags=WITHOUT_TAGS`: exclude any run with one of specified
  tags (separated by commas)
- `--create-xml`: create xml from study (current directory has to be
  a study)

Examples
--------

- read `sample.xml`, create **destination** directory and exit;
  ```
  $ code_saturne studymanager -f sample.xml
  ```
- copy all cases from the **repository** into the **destination**,
  compile all user files and run enabled cases:
  ```
  $ code_saturne studymanager -f sample.xml -r
  ```
- as above, and compare all new checkpoint files with those from the
  **repository** if defined in `sample.xml`
  ```
  $ code_saturne smgr -f sample.xml -r -c
  ```
- as above, and plots results if defined in `sample.xml`
  ```
  $ code_saturne smgr -f sample.xml -rcp
  ```
- as above, and send the two reports:
  ```
  $ code_saturne smgr -f sample.xml -r -c -p -m "dt@moulinsart.be dd@moulinsart.be"
  ```
- compare and plot results in the **destination** already computed
   ```
   $ code_saturne smgr -f sample.xml -c -p
   ```
- run cases tagged "coarse" (standing for coarse mesh for example)
  _and_ "hr" (standing for high Reynolds for example) only for 2 time
  iterations in destination directory of path `../RUNS/RIBS`
  (`RIBS} will be created, `RUNS` already exists). The command is
  launched from inside the study directory, so the repository
  containing the original study is simply indicated by `..`
  ```
  $ code_saturne smgr -f smgr_ribs.xml -r -n 2 --with-tags=coarse,hr --dest=../RUNS/RIBS --repo=..
  ```

### Note

The detailed report is generated only if the options `-c, --compare`
or `-p, --post` is present in the command line.

Parameters file
===============

The parameters file is an XML (text) file.

Begin and end of the parameters file
---------------------------------------

This example shows the four mandatory first lines of the parameters file.

```{.xml}
<?xml version="1.0"?>
<studymanager>
    <repository>/home/dupond/codesaturne/MyRepository</repository>
    <destination>/home/dupond/codesaturne/MyDestination</destination>
```

The third and fourth lines correspond to the definition of the
**repository** and **destination** directories.
Inside the markups `<repository>` and `<destination>` the user
must indicate the related directories. If the **destination** does
not exist, the directory is created.

The last line of the parameters file must be:

```{.xml}
</studymanager>
```

Case creation and compilation of the user files
-----------------------------------------------

When STUDYMANAGER is launched, the parameters file is parsed in order to
known which studies and cases from the **repository** should be copied
in the **destination**. The selection is done with the markups
`<study>` and `<case>` as in the following example:

```{.xml}
<?xml version="1.0"?>
<studymanager>
    <repository>/home/dupond/codesaturne/MyRepository</repository>
    <destination>/home/dupond/codesaturne/MyDestination</destination>

    <study label="MyStudy1" status="on">
        <case label="Grid1" run_id="Grid1" status="on" compute="on" post="off"/>
        <case label="Grid2" run_id="Grid2" status="off" compute="on" post="off"/>
    </study>
    <study label="MyStudy2" status="off">
        <case label="k-eps" status="on" compute="on" post="off"/>
        <case label="Rij-eps" status="on" compute="on" post="off"/>
    </study>
</studymanager>
```

The attributes are:
- `label`: the name of the file of the script;
- `status`: must be `on` or `off`,
  activate or deactivate the markup;
- `compute`: must be `on` or `off`,
  activate or deactivate the computation of the case;
- `post: must be `on` or `off`,
  activate or deactivate the post-processing of the case;
- `run_id`: name of the run directory (sub-directory of `RESU`) in
  which the result is stored. This attribute is optional. If it is
  not set (or if set to `run_id=""`), an automatic value will be
  proposed by the code (usually based on current date and time).
- `tags`: possible tags distinguishing the run from the others in
  the same XML parameter file (ex.: `tags="coarse,high-reynolds"`).

Only the attributes `label`, `status`, `compute`, and `post` are
mandatory.

If the directory specified by the attribute `run_id` already exists,
the computation is not performed again. For the post-processing step,
the existing results are taken into account only if no error file is
detected in the directory.

With the `status` attribute, a single case or a complete study can be
switched off. In the above example, only the case `Grid1` of the study
`MyStudy1` is going to be created.

After the creation of the directories in the **destination**, for each
case, all user files are compiled. The STUDYMANAGER stops if a compilation
error occurs: neither computation nor comparison nor plot will be performed,
even if they are switched on.

### Notes

- During the duplication (copy), all files are copied, except mesh files,
  for which a symbolic link is used.
- During the duplication, if a file already exists in the
  **destination**, this file is not overwritten by STUDYMANAGER
  by default.

Run cases {#sec_smgr_run}
---------

The computations are activated if the option `-r, --run` is present in
the command line.

All cases described in the parameters file with the attribute
`compute="on"`` are taken into account.

```{.xml}
<?xml version="1.0"?>
<studymanager>
    <repository>/home/dupond/codesaturne/MyRepository</repository>
    <destination>/home/dupond/codesaturne/MyDestination</destination>

    <study label="MyStudy1" status="on">
        <case label="Grid1" status="on" compute="on" post="off"/>
        <case label="Grid2" status="on" compute="off" post="off"/>
    </study>
    <study label="MyStudy2" status="on">
        <case label="k-eps" status="on" compute="on" post="off"/>
        <case label="Rij-eps" status="on" compute="on" post="off"/>
    </study>
</studymanager>
```

After the computation, if no error occurs, the attribute `compute} is set
to `"off"` in the copy of the parameters file in the
**destination**. This allows restarting STUDYMANAGER without re-running
successful previous computations.

Note that it is possible to run several times the same case in a given study.
The case has to be repeated in the parameters file:

```{.xml}
<?xml version="1.0"?>
<studymanager>
    <repository>/home/dupond/codesaturne/MyRepository</repository>
    <destination>/home/dupond/codesaturne/MyDestination</destination>

    <study label="MyStudy1" status="on">
        <case label="CASE1" run_id="Grid1" status="on" compute="on" post="on">
            <prepro label="grid.py" args="-m grid1.med -p cas.xml" status="on"/>
        </case>
        <case label="CASE1" run_id="Grid2" status="on" compute="on" post="on"/>
            <prepro label="grid.py" args="-m grid2.med -p cas.xml" status="on"/>
        </case>
    </study>
</studymanager>
```

If nothing is done, the case is repeated without modifications. In order to modify
the setup between two runs of the same case, an external script has to be used to
change the related setup (see [preprocessing](@ref sec_smgr_prepro)} and
[tricks](@ref sec_smgr_tricks) sections).

Compare checkpoint files
------------------------

The comparison is activated if the option `-c`, or `--compare` is present in
the command line.

In order to compare two checkpoint files for a given case, a markup
`<compare>` must be added as child of the considered case.
In the following example, a checkpoint file comparison is switched on for the
case _Grid1_ (for all
variables, with the default threshold), whereas no comparison is planned for
the case _Grid2_. The comparison is done by the same mechanism
as the `code_saturne bdiff` command.

```{.xml}
<study label='MyStudy1' status='on'>
    <case label='Grid1' status='on' compute="on" post="off">
        <compare dest="" repo="" status="on"/>
    </case>
    <case label='Grid2' status='on' compute="off" post="off"/>
</study>
```

The attributes are:

- `repo`: id of the results directory in the **repository** for
  example `repo="20110704-1116"`; if there is a single results directory
  in the `RESU` directory of the case, the id can be omitted:
  `repo=""`;

- `dest`: id of the results directory in the **destination**:
  * If the id is not known already because the case has not yet run,
    just leave the attribute empty `dest=""`, and the value will be
    updated after the run step in the **destination** directory
    (see section about [restart](@ref sec_smgr_restart));
  * if STUDYMANAGER is restarted without the run step (with the command
    line `code_saturne studymanager -f sample.xml -c` for example),
    the id of the results directory in the **destination** must be
    given (for example `dest="20110706-1523"`), but if there is a
    single results directory in the `RESU` directory of the case,
    the id can be omitted: with `dest=""`, the id will be
    completed automatically;

- `args`: additional options for the `code_saturne bdiff` command
   or underlying `cs_io_dump` tool:
  *  `--section`: name of a particular variable;
  *  `--threshold`: real value above which a difference is considered
     significant (default: *1e<sup>-30</sup>* for all variables);

- `status`: must be equal to `on` or `off`, used to activate or
   deactivate the markup.

Only the `repo`, `dest` and `status` attributes are mandatory.

Several comparisons with different options are permitted:
```{.xml}
<study label='MyStudy1' status='on'>
    <case label='Grid1' status='on' compute="on" post="off">
        <compare dest="" repo="" args="--section Pressure --threshold=1000" status="on"/>
        <compare dest="" repo="" args="--section VelocityX --threshold=1e-5" status="on"/>
        <compare dest="" repo="" args="--section VelocityY --threshold=1e-3" status="on"/>
    </case>
</study>
```

Comparison results will be summarized in a table in the file
`report_detailed.pdf`  (see [output](@ref sec_smgr_restart) section):

<table>
<tr><th> Variable Name <th> Diff. Max <th> Diff. Mean <th> Threshold
<tr><td> VelocityX     <td> 0.102701  <td> 0.00307058 <td> 1.0e-5
<tr><td> VelocityY     <td> 0.364351  <td> 0.00764912 <td> 1.0e-3
</table>

Alternatively, in order to compare all activated cases (status at on)
listed in a STUDYMANAGER parameter file, a reference directory can be
provided directly in the command line, as follows:

`code_saturne studymanager -f sample.xml -c -d /scratch/***/reference_destination_directory`

Run external preprocessing scripts with options {#sec_smgr_prepro}
-----------------------------------------------

The markup `<prepro>` can be added as a child of the considered case.

```{.xml}
<study label='STUDY' status='on'>
    <case label='CASE1' status='on' compute="on" post="on">
        <prepro label="mesh_coarse.py" args="-n 1" status="on"/>
    </case>
</study>
```

The attributes are:
- `label`: the name of the file of the considered script;
- `status`: must be `on` or `off`: activate or deactivate the markup;
- `args`: additional arguments to pass to the script.

Only the `label` and `status` attributes are mandatory.

An additional option `"-c"` (or `"--case"`) is given by
default with the path of the current case as argument (see example
in [tricks](@ref sec_smgr_tricks) for decoding options).

Note that all options must be processed by the script itself.

Several calls of the same script or to different scripts are permitted:
```{.xml}
<study label="STUDY" status="on">
    <case label="CASE1" status="on" compute="on" post="on">
        <prepro label="script_pre1.py" args="-n 1" status="on"/>
        <prepro label="script_pre2.py" args="-n 2" status="on"/>
    </case>
</study>
```

All preprocessing scripts are first searched in the `MESH` directory
from the current study in the **repository**. If a script is not found,
it is searched in the directories of the current case.
The main objective of running such external scripts is to create or modify
meshes or to modify the current setup of the related case (see
[tricks](@ref sec_smgr_tricks)).

Run external additional postprocessing scripts with options for a case
----------------------------------------------------------------------

The launch of external scripts is activated if the option `-p, --post`
is present in the command line.

The markup `<script>` has to be added as a child of the considered case.

```{.xml}
<study label='STUDY' status='on'>
    <case label='CASE1' status='on' compute="on" post="on">
        <script label="script_post.py" args="-n 1" dest="" repo="20110216-2147" status="on"/>
    </case>
</study>
```

The attributes are:
- `label`: the name of the file of the considered script;
- `status`: must be `on` or `off`: activate or deactivate the markup;
- `args`: the arguments to pass to the script;
- `repo` and `dest`: id of the results directory in the
  **repository** or in the **destination**;
  * If the id is not known already because the case has not yet run,
    just leave the attribute empty (`dest=""`), and the value will be
    updated after the run step in the **destination** directory
    (see [output](@ref sec_smgr_restart) section).
  * If there is a single results directory in the `RESU` directory
    (either in the **repository** or in the **destination**)
    of the case, the id can be omitted: `repo=""` or `dest=""`,
    the id will be completed automatically.
  If attributes `repo` and `dest` exist, their associated value
  will be passed to the script as arguments, with options `"-r"` and
  `"-d"` respectively.

Only the `label` and `status` attributes are mandatory.

Several calls of the same script or to different scripts are permitted:
```{.xml}
<study label="STUDY" status="on">
    <case label="CASE1" status="on" compute="on" post="on">
        <script label="script_post.py" args="-n 1" status="on"/>
        <script label="script_post.py" args="-n 2" status="on"/>
        <script label="script_post.py" args="-n 3" status="on"/>
        <script label="another_script.py" status="on"/>
    </case>
</study>
```

All postprocessing scripts must be in the `POST` directory from
the current study in the **repository**.
The main objective of running external scripts is to create or modify
results in order to plot them.

Example of a script which searches printed information in the listing,
(note the function to process the passed command line arguments):
```{.py}
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import string
from optparse import OptionParser

def process_cmd_line(argv):
    """Processes the passed command line arguments."""
    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-r", "--repo", dest="repo", type="string",
                      help="Directory of the result in the repository")

    parser.add_option("-d", "--dest", dest="dest", type="string",
                      help="Directory of the result in the destination")

    (options, args) = parser.parse_args(argv)
    return options

def main(options):
    m = os.path.join(options.dest, "run_solver.log")
    f = open(m)
    lines = f.readlines()
    f.close()

    g = open(os.path.join(options.dest, "water_level.dat"), "w")
    g.write("# time,   h_sim,   h_th\n")
    for l in lines:
       if l.rfind("time, h_sim, h_th") == 0:
           d = l.split()
           g.write("%s  %s  %s\n" % (d[3], d[4], d[5]))
    g.close()

if __name__ == '__main__':
    options = process_cmd_line(sys.argv[1:])
    main(options)
```

Run external additional postprocessing scripts with options for a study
-----------------------------------------------------------------------

The launch of external scripts is activated if the option `-p` or `--post`
is present in the command line.

The purpose of this functionality is to create new data based on several runs of
cases, and to plot them (see [2D plots](@ref sec_smgr_curves)) or to insert them
in the final detailed report (see [post-processing input](@ref sec_smgr_input)).

The `<postpro>` markup must to be added as a child of the considered study.

```{.xml}
<study label='STUDY' status='on'>
    <case label='CASE1' status='on' compute="on" post="on"/>
    <postpro label='Grid2.py' status="on" arg="-n 100">
        <data file="profile.dat">
            <plot fig="1" xcol="1" ycol="2" legend="Grid level 2" fmt='b-p'/>
            <plot fig="2" xcol="1" ycol="3" legend="Grid level 2" fmt='b-p'/>
        </data>
    <input file="output.dat" dest=""/>
    </postpro>
</study>
```

The attributes are:
- `label`: the name of the file of the considered script;
- `status`: must be `on}` or `off`: activate or deactivate the markup;
- `args`: the additional options to pass to the script;

Only the `label` and `status` attributes are mandatory.

The options given to the script in the command line are:

- `-s` or `--study`: label of the current study;
- `-c` or `--cases`: string which contains the list of the cases
- `-d` or `--directories`: string which contains the list
   of the directories of results.

Additional options can be passed to the script through the `args`
attributes.

Note that all options must be processed by the script itself.

Several calls of the same script or to different scripts are permitted.

Post-processing: 2D plots {#sec_smgr_curves}
-------------------------

The post-processing is activated if the option `-p` or `--post` is present
in the command line.

The following example shows the drawing of four curves (or plots, or 2D lines)
from two files of data (which have the same name `profile.dat`). There
are two subsets of curves (i.e. frames with axis and 2D lines), in a single
figure. The figure will be saved on the disk in a **pdf** (default)
or **png** format, in the `POST` directory of the related study
in the **destination**. Each drawing of a single curve is defined as a
markup child of a file of data inside a case. Subsets and figures are defined
as markup children of `<study>`.

```{.xml}
<study label='Study' status='on'>
    <case label='Grid1' status='on' compute="off" post="on">
        <data file="profile.dat" dest="">
            <plot fig="1" xcol="1" ycol="2" legend="Grid level 1" fmt='r-s'/>
            <plot fig="2" xcol="1" ycol="3" legend="Grid level 1" fmt='r-s'/>
        </data>
    </case>
    <case label='Grid2' status='on' compute="off" post="on">
        <data file="profile.dat" dest="">
            <plot fig="1" xcol="1" ycol="2" legend="Grid level 2" fmt='b-p'/>
            <plot fig="2" xcol="1" ycol="3" legend="Grid level 2" fmt='b-p'/>
        </data>
    </case>
    <subplot id="1" legstatus='on' legpos ='0.95 0.95' ylabel="U ($m/s$)" xlabel="Time ($s$)"/>
    <subplot id="2" legstatus='on' legpos ='0.95 0.95' ylabel="U ($m/s$)" xlabel="Time ($s$)"/>
    <figure name="velocity" idlist="1 2" figsize="(4,5)" format="png"/>
</study>
```

Define curves
-------------

The plots of computational data are built from data files. These data must be
ordered as column and the files should be in results directory in the
`RESU` directory (either in the **repository** or in the
**destination**). Lines starting with character `\#` are treated as
as comments.

In the parameters file, curves are defined with two markups:
`<data>` and `<plot>`:

- `<data>`: child of markup `<case>`, defines a file of data;
  * `file`: name of the file of data
  * `repo` or `dest`: id of the results directory either in the
    **repository** or in the **destination**;
    - If the id is not known already because the case has not yet run,
      just leave the attribute empty, with `dest=""`, and the value
      will be updated after the run step in the **destination**
      directory (see [output](@ref sec_smgr_restart) section).
    - If there is a single results directory in the `RESU` directory
      (either in the **repository** or in the **destination**)
      of the case, the id can be ommitted: `repo=""` or `dest=""`,
      and it will be completed automatically.
  The `file` attribute is mandatory, and either `repo` or
  `dest` must be present (but not the both), even if they are empty.

- `<plot>`: child of markup `<data>`, defines a single curve;
  the attributes are:
  * `fig}` id of the subset of curves (i.e. markup `<subplot>`)
    where the current curve should be plotted;
  * `xcol`: number of the column in the file of data for the abscissa;
  * `ycol`: number of the column in the file of data for the ordinate;
  * `legend`: add a label to a curve;
  * `fmt`: format of the line, composed from a symbol, a color and a
    linestyle, for example `fmt="r--"` for a dashed red line;
  * `xplus`: real to add to all values of the column `xcol`;
  * `yplus`: real to add to all values of the column `ycol`;
  * `xfois`: real to multiply to all values of the column `xcol`;
  * `yfois`: real to multiply to all values of the column `ycol`;
  * `xerr` or `xerrp`: draw horizontal error bar (see section
    on [error bars](@ref sec_smgr_err));
  * `yerr` or `yerrp`: draw vertical error bar (as above);
  * some standard options of 2D lines can be added, for example
    `markevery="2"` or `markersize="3.5"`. These options
    are summarized in the table @ref smgr_table_curves. Note that
    the options which are string of characters must be encased in
    quotes like this: `color="'g'"`.

<table>
<caption id="smgr_table_curves">Options authorized as attributes of the markup `plot`</caption>
<tr><th> Property <th> Value Type
<tr><td> alpha <td> float (0.0 transparent through 1.0 opaque)
<tr><td> antialiased or aa <td> `True` or `False`
<tr><td> color or c <td> any Matplotlib color
<tr><td> dash_capstyle <td> `butt`; `round`; `projecting`
<tr><td> dash_joinstyle <td> `miter`; `round`; `bevel`
<tr><td> dashes <td> sequence of on/off ink in points ex: `dashes="(5,3)"`
<tr><td> label <td> any string, same as legend
<tr><td> linestyle or ls <td>  `-`; `--`; `-.`; `:`; `steps`; ...
<tr><td> linewidth or lw <td> float value in points
<tr><td> marker <td>  `+`; `,`; `.`; `1`; `2`; `3`; `4`; ...
<tr><td> markeredgecolor or mec <td> any Matplotlib color
<tr><td> markeredgewidth or mew <td> float value in points
<tr><td> markerfacecolor or mfc <td> any Matplotlib color
<tr><td> markersize or ms <td> float
<tr><td> markevery <td> `None`; integer; (startind, stride)
<tr><td> solid_capstyle <td> `butt`; `round`; `projecting`
<tr><td> solid_joinstyle <td> `miter`; `round`; `bevel`
<tr><td> zorder <td> any number
</table>

The attributes `fig` and `ycol` are mandatory.

In case a column should undergo a transformation specified by the attributes
`xfois`,`yfois`,`xplus`,`yplus`, scale operations
take precedence over translation operations.

Details on 2D lines properties can be found in the
[Matplotlib documentation](https://matplotlib.org/contents.html).
For more advanced options see [Matplotlib raw commands](@ref sec_smgr_raw).

### Define subsets of curves

A subset of curves is a frame with two axis, axis labels, legend, title and
drawing of curves inside. Such subset is called subplot in
the Matplotlib vocabulary.

`<subplot>`: child of markup `<study>`, defines a frame with
several curves; the attributes are:
- `id`: id of the subplot, should be an integer;
- `legstatus`: if `"on"` display the frame of the legend;
- `legpos`: sequence of the relative coordinates of the center of
   the legend, it is possible to draw the legend outside the axis;
- `title`: set title of the subplot;
- `xlabel`: set label for the x axis;
- `ylabel`: set label for the y axis;
- `xlim`: set range for the x axis;
- `ylim`: set range for the y axis.

For more advanced options see [Matplotlib raw commands](@ref sec_smgr_raw).

### Define figures

Figure is a compound of subset of curves.

`<figure>`: child of markup `<study>`, defines a pictures
with a layout of frames; the attributes are:
- `name`: name of the file to be written on the disk;
- `idlist`: list of the subplot to be displayed in the figure;
- `title`: add a title on the top of the figure;
- `nbrow`: impose a number of rows of the layout of the subplots;
- `nbcol`: impose a number of columns of the layout of the subplots;
- `format`: format of the file to be written on the disk,
  `"pdf"` (default) or `"png"`;
  Other formats could be chosen (eps, ps, svg,...), but the pdf generation
  with pdflatex will not be possible in this case;
- `figsize`: width x height in inches; defaults to (4,4);
- `dpi`: resolution; defaults to 200 if format is set to pdf;
  or to 800 if format is set to png; only customizable for
  png format.

The `name` and `idlist` attributes are mandatory.

### Experimental or analytical data

A particular markup is provided for curves of experimental or analytical data:
`<measurement>`; the attributes are:
- `file`: name of the file to be read on the disk;
- `path`: path of the directory where the data file is located.
  The path can be omitted (`path=""`), and in this case, the file
  will be searched recursively in the directories of the considered study.

The `file` and `path` attributes are mandatory.

In order to draw curves of experimental or analytical data, the `<measurement`
markup should be used with the markup `<plot>` as illustrated below:

```{.xml}
<study label='MyStudy' status='on'>
    <measurement file='exp1.dat' path=''>
            <plot fig='1' xcol='1' ycol='2' legend='U Experimental data'/>
            <plot fig='2' xcol='3' ycol='4' legend='V Experimental data'/>
    </measurement>
    <measurement file='exp2.dat' path =''>
            <plot fig='1' xcol='1' ycol='2' legend='U Experimental data'/>
            <plot fig='2' xcol='1' ycol='3' legend='V Experimental data'/>
    </measurement>
    <case label='Grid1' status='on' compute="off" post="on">
        <data file="profile.dat" dest="">
            <plot fig="1" xcol="1" ycol="2" legend="U computed" fmt='r-s'/>
            <plot fig="2" xcol="1" ycol="3" legend="V computed" fmt='b-s'/>
        </data>
    </case>
</study>
<subplot id="1" legstatus='on'  ylabel="U ($m/s$)" xlabel= "$r$ ($m$)" legpos ='0.05 0.1'/>
<subplot id="2" legstatus='off' ylabel="V ($m/s$)" xlabel= "$r$ ($m$)"/>
<figure name="MyFigure" idlist="1 2"  figsize="(4,4)" />
```

### Curves with error bars {#sec_smgr_err}

In order to draw horizontal and vertical error bars, it is possible to
specify in the markup `<plot>` the attributes `xerr` and
`yerr` respectively (or `xerrp` and `yerrp`). The
value of these attributes could be:
- the number of the column, in the file of data, that contains the total
  absolute uncertainty spans:
  ```{.xml}
  <measurement file='axis.dat' path =''>
      <plot fig='1' xcol='1' ycol='3' legend='Experimental data' xerr='2' />
  </measurement>
   ```
- the numbers of the two columns, in the file of data, that contain the
  absolute low spans and absolute high spans of uncertainty:
  ```{.xml}
  <data file='profile.dat' dest="">
      <plot fig='1' xcol='1' ycol='2' legend='computation' yerr='3 4' />
  </data>
  ```
- a single real value equal to the percentage of uncertainty that should be
  applied to the considered data set:
  ```{.xml}
  <data file='profile.dat' dest="">
      <plot fig='1' xcol='1' ycol='2' legend='computation' yerrp='2.' />
  </data>
  ```
### Monitoring points or probes

A particular markup is provided for curves of probes data:
`<probes>`; the attributes are:

- `file`: name of the file to be read on the disk;
- `fig`: id of the subset of curves (i.e. markup `<subplot>`)
   where the current curve should be plotted;
- `dest`: id of the results directory in the **destination**:
  * If the id is not known already because the case has not yet run, just
    leave the attribute empty, with `dest=""`, and the value will be updated
    after the run step in the **destination** directory
    (see [output](@ref sec_smgr_restart) section).
  * If STUDYMANAGER is restarted without the run step (with the command
    line `code_saturne studymanager -f sample.xml -c` for example), the id of
    the results directory in the **destination** must be given (for example
    `dest="20110706-1523"`), but if there is a single results directory in
    the `RESU` directory of the case, the id can be omitted:
    with `dest=""`, the id will be completed automatically.

The `file`, `fig` and `dest` attributes are mandatory.

In order to draw curves of probe data, the `<probes>` markup
should be used as a child of a markup `<case>` as illustrated below:

```{.xml}
<study label='MyStudy' status='on'>
    <measurement file='exp1.dat' path=''>
        <plot fig='1' xcol='1' ycol='2' legend='U Experimental data'/>
    </measurement>
    <case label='Grid1' status='on' compute="off" post="on">
        <probes file="probes_U.dat" fig ="2" dest="">
        <data file="profile.dat" dest="">
            <plot fig="1" xcol="1" ycol="2" legend="U computed" fmt='r-s'/>
        </data>
    </case>
</study>
<subplot id="1" legstatus='on'  ylabel="U ($m/s$)" xlabel= "$r$ ($m$)" legpos ='0.05 0.1'/>
<subplot id="2" legstatus='on'  ylabel="U ($m/s$)" xlabel= "$time$ ($s$)" legpos ='0.05 0.1'/>
<figure title="Results" name="MyFigure" idlist="1"/>
<figure title="Grid1: probes for velocity"  name="MyProbes" idlist="2"/>
```

### Matplotlib raw commands {#sec_smgr_raw}

The parameters file allows executing additional Matplotlib commands (i.e
Python commands), for curves (2D lines), or subplot, or figures. For every object
drawn, `studymanager` associate a name to this object that can be reused in
standard Matplotlib commands. Therefore, children markup `<plt_command>`
could be added to `<plot>`, `<subplot>` or `<figure>`.

It is possible to add commands with the **Matlab style** or **Python
style**. For the Matlab style, commands are called as methods of the
`plt` module, and for the Python style commands or called as methods of the
instance of the graphical object.

Matlab style and Python style commands can be mixed.

- curves or 2D lines: when a curve is drawn, the associated name
  are `line` and `lines` (with `line = lines[0]`).

  ```{.xml}
  <plot fig="1" xcol="1" ycol="2" fmt='g^' legend="Simulated water level">
      <plt_command>plt.setp(line, color="blue")</plt_command>
      <plt_command>line.set_alpha(0.5)</plt_command>
  </plot>
  ```

- subset of curves (subplot): for each subset, the associated name is `ax`:

  ```{.xml}
  <subplot id="1" legend='Yes' legpos ='0.2 0.95'>
      <plt_command>plt.grid(True)</plt_command>
      <plt_command>plt.xlim(0, 20)</plt_command>
      <plt_command>ax.set_ylim(1, 3)</plt_command>
      <plt_command>plt.xlabel(r"Time ($s$)", fontsize=8)</plt_command>
      <plt_command>ax.set_ylabel(r"Level ($m$)", fontsize=8)</plt_command>
      <plt_command>for l in ax.xaxis.get_ticklabels(): l.set_fontsize(8)</plt_command>
      <plt_command>for l in ax.yaxis.get_ticklabels(): l.set_fontsize(8)</plt_command>
      <plt_command>plt.axis([-0.05, 1.6, 0.0, 0.15])</plt_command>
      <plt_command>plt.xticks([-3, -2, -1, 0, 1])</plt_command>
  </subplot>
```

Post-processing: input files {#sec_smgr_input}
----------------------------

The post-processing is activated if the option `-p, --post` is present
in the command line.

STUDYMANAGER can include files into the final detailed report. These
files must be in the directory of results either in the **destination** or
in the **repository**. The following example shows the inclusion of three
files: `performance.log` and `setup.log` from the
**destination**, and a `performance.log` from the **repository**:

```{.xml}
<case label='Grid1' status='on' compute="on" post="on">
    <input dest="" file="performance.log"/>
    <input dest="" file="setup.log"/>
    <input repo="" file="performance.log"/>
</case>
```

Text files, \f$ \mbox{\LaTeX} \f$ source files, or graphical (PNG, JPEG, or PDF) files
may be included.

In the parameters file, input files are defined with markups `<input>`
as children of a single markup `<case>`.
The attributes of `<input>` are:
- `file`: name of the file to be included
- `repo` or `dest`: id of the results directory either in the
  **repository** or in the **destination**;
  * If the id is not known already because the case has not yet run, just leave
    the attribute empty, with `dest=""`; the value will be updated after the run
    step in the **destination** directory
    (see [output](@ref sec_smgr_restart) section).
  * If there is a single results directory in the `RESU` directory
    (either in the **repository** or in the **destination**) of the case,
    the id can be omitted: with `repo=""` or `dest=""`, the id will be
    completed automatically.

The `file` attribute is mandatory, and either `repo` or
`dest` must be present (but not the both) even if it is empty.

Output and restart {#sec_smgr_restart}
==================

STUDYMANAGER produces several files in the **destination** directory:

- `report.txt`: standard output of the script;
- `auto_vnv.log`: log of the code and the `pdflatex` compilation;
- `report_global.pdf`: summary of the compilation, run, comparison,
   and plot steps;
- `report_detailed.pdf`: details the comparison and display the
   plot;
- `sample.xml`: udpated parameters file, useful for restart the
   script if an error occurs.

After the computation of a case, if no error occurs, the attribute
`compute` is set to `"off"` in the copy of the parameters file
in the **destination**. It is allow a restart of STUDYMANAGER without
re-run successful previous computations.
In the same manner, all empty attributes `repo=""` and `dest=""`
are completed in the updated parameters file.

Tricks {#sec_smgr_tricks}
======

## Syntax and additional markup

### How to comment markup in the parameters file?

The opening and closing markup for comments in XML are `<!--` and
`-->`. In the following example, nothing from the study
`MyStudy2` will be read:
```{.xml}
<?xml version="1.0"?>
<studymanager>
    <repository>/home/dupond/codesaturne/MyRepository</repository>
    <destination>/home/dupond/codesaturne/MyDestination</destination>

    <study label="MyStudy1" status="on">
        <case label="Grid1" status="on" compute="on" post="on"/>
        <case label="Grid2" status="on" compute="off" post="on"/>
    </study>
    <!--
    <study label="MyStudy2" status="on">
        <case label="k-eps" status="on" compute="on" post="on"/>
        <case label="Rij-eps" status="on" compute="on" post="on"/>
    </study>
    -->
</studymanager>
```

### How to add text in a figure?

It is possible to use raw commands:
```{.xml}
<subplot id='301' ylabel ='Location ($m$)' title='Before jet -0.885' legstatus='off'>
    <plt_command>plt.text(-4.2, 0.113, 'jet')</plt_command>
    <plt_command>plt.text(-4.6, 0.11, r'$\downarrow$', fontsize=15)</plt_command>
</subplot>
```

### Adjust margins for layout of subplots in a figure.

You have to use the raw command `subplots_adjust`:

```{.xml}
<subplot id="1" legend='Yes' legpos ='0.2 0.95'>
    <plt_command>plt.subplots_adjust(hspace=0.4, wspace=0.4, right=0.9,
                 left=0.15, bottom=0.2, top=0.9)</plt_command>
</subplot>
```

### How to find a syntax error in the XML file?

  When there is a typo in the parameters file,
  STUDYMANAGER indicates the location of the error
  with the line and the column of the file:
```
my_case.xml file reading error.

This file is not in accordance with XML specifications.

The parsing syntax error is:

my_case.xml:86:12: not well-formed (invalid token)
```

### How to render less-than and greater-than signs in legends, titles or axis labels?

The less-than < and greater-than > symbols are among the five predefined
entities of the XML specification that represent special characters.

In order to have one of the five predefined entities rendered in any legend, title or axis label,
use the string `&name;`. Refer to the following table for the name of the character to be rendered:

<table>
<caption id="smgr_table_xml_spe_sym">Special symbols of the XML specification</caption>
<tr><th> name <th> character <th> description
<tr><td> `quot` <td> " <td> double quotation mark
<tr><td> `amp` <td> & <td> ampersand
<tr><td> `apos` <td> ' <td> apostrophe
<tr><td> `lt` <td> < <td>  less-than sign
<tr><td> `gt` <td> > <td>  greater-than sign
</table>

For any of these predefined entities, the XML parser will first replace
the string `&name;` by the character, which will then allow \f$ \mbox{\LaTeX} \f$
(or Mathtext if \f$ \mbox{\LaTeX} \f$ is disabled) to process it.

For example, in order to write \f$ \lambda<1 \f$ in a legend, the following attribute will be used:
```
  <plot fig="4" fmt="k--" legend="solution for $\lambda &lt; 1$" xcol="1" ycol="2"/>
```

### How to set a logarithmic scale?

The following raw commands can be used:

```{.xml}
<subplot id="2" title="Grid convergence" xlabel="Number of cells" ylabel="Error (\%)">
    <plt_command>ax.set_xscale('log')</plt_command>
    <plt_command>ax.set_yscale('log')</plt_command>
</subplot>
```

### How to create a mesh automatically with SALOME?

The following example shows how to create a mesh with a SALOME command file:
```{.xml}
<study label="STUDY" status="on">
    <case label="CASE1" status="on" compute="on" post="on">
        <prepro label="salome" args="-t -u my_mesh.py" status="on"/>
    </case>
</study>
```
with the `salome` command (depending of the local installation of
SALOME):
```{.sh}
#!/bin/bash

export ROOT_SALOME=/home/salome/Salome-V9_3_0-x86_64-univ
source /home/salome/Salome-V9_3_0-x86_64-univ/salome_prerequisites.sh
source /home/salome//Salome-V9_3_0-x86_64-univ/salome_modules.sh

/home/salome/appli_V9_3_0/salome $*
```
and the script of SALOME commands `my_mesh.py`:
```{.py}
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import geompy
import smesh

# create a box
box = geompy.MakeBox(0., 0., 0., 100., 200., 300.)
idbox = geompy.addToStudy(box, "box")

# create a mesh
tetra = smesh.Mesh(box, "MeshBox")

algo1D = tetra.Segment()
algo1D.NumberOfSegments(7)

algo2D = tetra.Triangle()
algo2D.MaxElementArea(800.)

algo3D = tetra.Tetrahedron(smesh.NETGEN)
algo3D.MaxElementVolume(900.)

# compute the mesh
tetra.Compute()

# export the mesh in a MED file
tetra.ExportMED("./my_mesh.med")
```

## How to carry out a grid convergence study?

The following example shows how to carry out a grid convergence study by running
the same case three times and changing the parameters between each run with the
help of a preprocessing script.

Here the mesh, the maximum number of iterations, the reference time step and the
number of processes are modified, before each run, by the
`prepro.py` script.

The parameters file is as follows:

```{.xml}
<case compute="on" label="COUETTE" post="on" run_id="21_Theta_1" status="on">
    <prepro args="-m 21_Theta_1.med -p Couette.xml -n 4000 -a 1. -t 0.01024 -u 1"
            label="prepro.py" status="on"/prepro>
    <data dest="" file="profile.dat">
        <plot fig="5" fmt="r-+" legend="21 theta 1" markersize="5.5" xcol="1" ycol="5"/>
    </data>
</case>

<case compute="on" label="COUETTE" post="on" run_id="43_Theta_05" status="on">
    <prepro args="-m 43_Theta_05.med -p Couette.xml -n 8000 -a 0.5 -t 0.00512 -u 2"
            label="prepro.py" status="on"/prepro>
    <data dest="" file="profile.dat">
        <plot fig="5" fmt="b" legend="43 Theta 05" markersize="5.5" xcol="1" ycol="5"/>
    </data>
</case>

<case compute="on" label="COUETTE" post="on" run_id="86_Theta_025" status="on">
    <prepro args="-m 86_Theta_025.med -p Couette.xml -n 16000 -a 0.25 -t 0.00256 -u 4"
            label="prepro.py" status="on" /prepro>
    <data dest="" file="profile.dat">
        <plot fig="5" fmt="g" legend="86 Theta 025" markersize="5.5" xcol="1" ycol="5"/>
    </data>
</case>
```

Recall that the case attribute `run`_id` should be given a different
value for each run, while the `label` should stay the same and that the
preprocessing script should be copied in the `MESH` direcrtory of the study or in
the `DATA` directory of the case.

The preprocessing script is shown below. Note that it can be called inside
the parameters file without specifying a value for each option:

```{.py}
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys
import string
from optparse import OptionParser

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------
from model.ScriptRunningModel import ScriptRunningModel

#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """Processes the passed command line arguments."""
    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-c", "--case", dest="case", type="string",
                      help="Directory of the current case")

    parser.add_option("-p", "--param", dest="param", type="string",
                      help="Name of the parameters file")

    parser.add_option("-m", "--mesh", dest="mesh", type="string",
                      help="Name of the new mesh")

    parser.add_option("-n", "--iter-num", dest="iterationsNumber", type="int",
                      help="New iteration number")

    parser.add_option("-t", "--time-step", dest="timeStep", type="float",
                      help="New time step")

    parser.add_option("-a", "--perio-angle", dest="rotationAngle", type="float",
                      help="Periodicity angle")

    (options, args) = parser.parse_args(argv)

    return options

#-------------------------------------------------------------------------------

def main(options):
    from cs_package import package
    from model.XMLengine import Case
    from model.XMLinitialize import XMLinit
    from model.SolutionDomainModel import SolutionDomainModel
    from model.TimeStepModel import TimeStepModel

    fp = os.path.join(options.case, "DATA", options.param)
    if os.path.isfile(fp):
        try:
            case = Case(package = package(), file_name = fp)
        except:
            print("Parameters file reading error.\n")
            print("This file is not in accordance with XML specifications.")
            sys.exit(1)

        case['xmlfile'] = fp
        case.xmlCleanAllBlank(case.xmlRootNode())
        XMLinit(case).initialize()

        if options.mesh:
            s = SolutionDomainModel(case)
            l = s.getMeshList()
            s.delMesh(l[0])
            s.addMesh((options.mesh, None))

        if options.rotationAngle:
            s.setRotationAngle(0, options.rotationAngle)

        if (options.iterationsNumber):
            t = TimeStepModel(case)
            t.setStopCriterion('iterations', options.iterationsNumber)

        if (options.TimeStep):
            t = TimeStepModel(case)
            t.setTimeStep(options.TimeStep)

        case.xmlSaveDocument()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    options = process_cmd_line(sys.argv[1:])
    main(options)

#-------------------------------------------------------------------------------
```
