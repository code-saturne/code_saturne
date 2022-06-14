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

<!--

Colors
------

blueedf          rgb(0,91,187)
blueedf_light    rgb(0,91,187)
blueedf_dark     rgb(9,53,122)

greenedf         rgb(0,158,47)
greenedf_light   rgb(127,175,16)
greenbed_dark    rgb(48,119,16)

orangeedf        rgb(255,160,47)
orangeedf_light  rgb(255,122,0)
orangeedf_dark   rgb(254,88,21)

antiquewhite  rgb(0.98, 0.92, 0.84)
whitesmoke    rgb(0.96, 0.96, 0.96)

-->

\page cs_dg_debugging Debugging

[TOC]

Introduction
============

All non-trivial software has bugs, at least at first,
so debugging is a part of code development.

- Is is often easier to debug newly written code than older code

  * At least when done by the same person

  * Test and debug early, before you forget the details of your code

- Multiple tools and techniques are available

  * Mastering debugging tools and techniques make debugging less painful

- Bugs detected late are much more costly then those detected early

  * May require re-validation

  * Often harder to debug, because the code details need to be ``re-learned''

  * Be *proactive*, try to detect as many bugs as possible by testing

Debugging Tools
===============

Debugging methods
-----------------

When encountering or suspecting a bug, choosing the best debugging technique for
the situation can lead to one or more orders of magnitude in time savings.

- Choosing which tool to use is the difficult part, and
  can be based on several factors:

  * Comparison of experience to observed effects

  * Leveraging theoretical knowledge and experience

    - Guessing at whether bugs may be due to uninitialized values,
      out-of-bounds arrays accesses, bad option settings, numerical errors, ...

  * Sometimes, a bit of luck...

- code_saturne tries to help, so always check for
  <span style="color:rgb(254,88,21)"> **error messages** </span>

  * In code_saturne, <span style="color:rgb(48,119,16)"><b>`error*`</b></span>,
    <span style="color:rgb(48,119,16)"><b>`run_solver.log`</b></span>,
    messages in batch output logs (or in the console) should be checked.

    - For parallel runs, when both, <span style="color:rgb(48,119,16)">`error`</span>
      and <span style="color:rgb(48,119,16)">`error_r*`</span> are present,
      the latter are the ones which contain the useful information.

  * See [section in user guide](@ref sec_ug_troubleshhoting) for more details.

  * Some graphical checks with
    <span style="color:rgb(48,119,16)"><b>`postprocessing/error*`</b></span>
    outputs are also available for boundary conditions and linear solvers.

### Some debugging tools that should be considered

-  Source code proofreading

   * (Re-) <span style="color:rgb(255,160,47)"><b>checking for compiler warnings</b></span>.

   * Interactive debugging

   * Memory debugging

   * Checks/instrumentation in code...

     - Use <span style="color:rgb(48,119,16)">`--enable-debug`</span>
       to configure builds for debug.

       * Enables use of many `assert` checks in C code.

       * Enables arrays bounds-checking in Fortran.

   * Using recent versions of the [GCC](https://gcc.gnu.org/)
     or [clang](https://clang.llvm.org/) compiler, compile/run with
     [AddressSanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer),
     [UndefinedBehaviorSanitizer](https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html),
     and other tools of this [series](https://github.com/google/sanitizers) frequently.

     - Code built this way not compatible with runs under
       [Valgrind](https://valgrind.org/).

     - Not compatible either with some resource limits set on some clusters.

     - Overhead: usually about x3.

- When you known where to search,
   <span style="color:rgb(255,160,47)"><b>print</b></span>
   statements may be useful...

The GNU debugger
----------------

The GNU debugger
[<https://www.gnu.org/software/gdb>](<https://www.gnu.org/software/gdb>)
is a broadly available, interactive debugger for compiled languages
including C, C++, and Fortran.

- To debug an executable, run
  <span style="color:rgb(48,119,16)">`gdb <executable>`</span>

  * Under the gdb prompt, type
    <span style="color:rgb(48,119,16)">`help`</span> for built-in help,
    <span style="color:rgb(48,119,16)">`q`</span> to quit.

  * Help is grouped in categories

    - The most common options are:

      * <span style="color:rgb(48,119,16)">`b`</span> (set breakpoint),

      * <span style="color:rgb(48,119,16)">`c`</span>
        (continue to next statement),

      * <span style="color:rgb(48,119,16)">`s`</span>
        (step into function),

      * <span style="color:rgb(48,119,16)">`p`</span>
        (print).

  * Many front-ends are available, including:

    - A built-in user interface for terminals.

    - integration with text editors, especially

      * <span style="color:rgb(48,119,16)">Emacs</span> (built-in),

      * <span style="color:rgb(48,119,16)">vim</span>
        ([Conque GDB](https://www.vim.org/scripts/script.php?script_id=4582),
        Termdebug).

      * [Neovim](https://neovim.io)
        ([NeoDebug](https://github.com/cpiger/NeoDebug),
        [nvim-gdb](https://github.com/sakhnik/nvim-gdb/blob/master/test/prerequisites.sh)).

    - Standalone graphical interfaces:

      * [gdbgui](https://www.gdbgui.com/)

      * [KDbg](https://www.kdbg.org/),

      * [Nemiver](https://www.gnu.org/software/ddd),

      * [DDD](https://www.gnu.org/software/ddd),

    - integration in development environments:

      * [Kdevelop](https://www.kdevelop.org/)

      * [Anjuta](https://gitlab.gnome.org/GNOME/anjuta)

      * [Qt Creator](https://www.qt.io/download)

- GDB can provide some information on any compiled program, but provides
  more detailed and useful information when the program was compiled
  with debugging info. The matching compiler option is usually
  <span style="color:rgb(48,119,16)">`-g`</span>, and in the
  case of code_saturne, is provided using the
  <span style="color:rgb(48,119,16)">`--enable-debug`</span>
  configure option at installation.

### GDB basic interface

When used directly, GDB runs in a single terminal frame, as shown here.
Only the current line of code is shown, though the
<span style="color:rgb(48,119,16)">`list`</span>
command allows showing more.

\image html dg/gdb_screen.png "GDB in terminal mode" width=80%

When started with the
<span style="color:rgb(48,119,16)">`-tui`</span> option, GDB runs
in a split terminal, with source on top, commands on bottom.

- Using the <span style="color:rgb(255,160,47)">`CTRL+x+o`</span> key
  combination allows changing focus from one to the other.

- Using the <span style="color:rgb(255,160,47)">`CTRL+l`</span> key
  allows refreshing the display.

\image html dg/gdb_tui_screen.png "GDB with split screen" width=80%

GDB may also be run under Emacs, which provides syntax highlighting of
source code.

\image html dg/emacs_gud_screen.png "GDB under Emacs" width=55%

### Graphical front-end recommendations

Many graphical front-ends are available for gdb. When evaluating
a front-end, we recommend to check for the following features:

- **Must** provide a console to allow combining text-based commands
  with the graphical elements, or at least easily-accessible
  widgets in which watchpoints and expressions to print can be typed.

- Should allow some means (such as command-line options) to connect
  to a GDB server through a socket interface (more on this later)

The [DDD (Data Display Debugger)](https://www.gnu.org/software/ddd)
front-end is obsolete and uses a dated graphical toolkit, but has the
advantage of combining a command prompt with graphical tools, and is
very easy to use, so it might remain an option.

The [Nemiver](https://www.gnu.org/software/ddd) debugger also has a
GDB back-end. It offers a clean display, but lacks the possibility
of typing commands; everything must be done using the mouse and menus,
which is often tedious. The project seems abandoned.
[KDbg](https://www.kdbg.org/) is similar, slightly more practical,
but does not seem to have been very active since 2018.

The [gdbgui](https://www.gdbgui.com/) debugger seems promising, and
a good potential successor to DDD. It is based on a web-browser interface,

\image html dg/gdbgui_screen.png "gdbgui" width=80%

Full integrated development environments (including Qt Creator,
Visual Studio Code, Eclipse, Kdevelop, Anjuta) are outside the
scope of this documentation. Most members of the code_saturne development
team mostly use lighter, less integrated tools, so will not be able
to provide recommendations regarding their use.

GDB alternatives
----------------

The [Eclipse CDT](https://projects.eclipse.org/projects/tools.cdt)
and [Eclipse PTP (Parallel Tools Platform)](https://www.eclipse.org/ptp/downloads.php)
environments integrate debuggers, including a parallel debugger,
but may use a different syntax than "standalone" GDB, so they are not considered
here (though feedback and recommendations around these tools are welcome).

The [LLDB](https://lldb.llvm.org) debugger is an interesting competitor to GDB,
with a different (similar but more verbose) syntax. Is is not as widely
available yet, and is not yet handled by the code_saturne debug scripts,
though a user familiar with it could of course set it up.

The Valgrind tool suite
-----------------------

The [Valgrind](https://www.valgrind.org) tool suite allows the detection of
many memory management (and other) bugs.

- Dynamic instrumentation

  * No need for recompilation

    - Usable with any binary, but provides more info (i.e. code line numbers)
      with code compiled in debug mode

    - Depending on tool used, run time and memory overhead from 10-100x.

      * With default tool (<span style="color:rgb(48,119,16)">Memcheck</span>),
        10x30.

      * Use proactively, to detect bugs on small cases, before they
        become a problem in production cases.

Valgrind is easy to run:

- Prefix a standard command with <span style="color:rgb(48,119,16)">`valgrind`</span>

  * By default, uses the <span style="color:rgb(255,160,47)">`memcheck`</span> tool.

  * Tool may be changed using
    <span style="color:rgb(48,119,16)">`valgrind –tool=/cachegrind/callgrind/drd/massif/...`</span>

- Valgrind may be combined with GDB using its
  <span style="color:rgb(48,119,16)">`gdbserver`</span> mode.

  * To use this mode, call
    <span style="color:rgb(48,119,16)">`valgrind –vgdb-error=<number>`</span>

    - The number represents the number of errors after which the
       gdbserver is invoked (0 to start immediately).

\image html dg/valgrind_screen.png "Valgrind in a terminal" width=60%

GCC and clang sanitizers
========================

Recent versions of the LLVM clang and GCC compilers have additional
instrumentation options, allowing memory debugging with a lower overhead
than Valgrind.

Address Sanitizer
-----------------

For the most common errors, use
<span style="color:rgb(48,119,16)">AddressSanitizer</span>,
a fast memory error detector.

- For the code_saturne configure options, this means
  <span style="color:rgb(48,119,16)">`CFLAGS=-fsanitize=address`</span>
  <span style="color:rgb(48,119,16)">`FCFLAGS=-fsanitize=address`</span>
  <span style="color:rgb(48,119,16)">`LDFLAGS=-fsanitize=address`</span>

- This may sometimes require specifying
  <span style="color:rgb(48,119,16)">`export LD_LIBRARY_FLAGS=<path_to_compiler_libraries`</span>
  when the compiler is installed in a nonstandard path on older systems.

- On some machines, this may be unusable if memory resource limits
  are set (check using <span style="color:rgb(48,119,16)">`ulimit -c`</span>

- Note that the resulting code will not be usable under Valgrind.

- Uninitialized values are not detected by Address Sanitizer (but
  may be detected by UndefinedBehaviorSanitizer).

- Out-of-bounds errors for arrays on stack (fixed size, usually small)
  are not detected by Valgrind, but may be detected by AddressSanitizer.

- AddressSanitizer also includes a memory leak checker, which is useful
  but may also report errors due to system libraries, so to allow a "clean"
  exit, we may use:

    ```
    export ASAN_OPTIONS=detect_leaks=0
    ```

UndefinedBehaviorSanitizer
--------------------------

The <span style="color:rgb(48,119,16)">UndefinedBehaviorSanitizer</span>
instrumentation is also useful to detect other types of bugs, such
as division by zero, some memory errors, integer overflows, and more.

- This may sometimes require also specifying
  <span style="color:rgb(48,119,16)">`-lubsan`</span>
  and even in some cases specify
  <span style="color:rgb(48,119,16)">`LD_LIBRARY_FLAGS`</span>

- For the code_saturne configure options, this means
  <span style="color:rgb(48,119,16)">`CFLAGS=-fsanitize=undefined`</span>
  <span style="color:rgb(48,119,16)">`FCFLAGS=-fsanitize=undefined`</span>
  <span style="color:rgb(48,119,16)">`LDFLAGS=-fsanitize=undefined`</span>

- This may sometimes require specifying
  <span style="color:rgb(48,119,16)">`export LD_LIBRARY_FLAGS=<path_to_compiler_libraries`</span>
  as per AddressSanitizer.

- Note that only code compiled with those options is instrumented.

Application to code_saturne
===========================

Starting the code_saturne solver under a debugger  {#code_saturne_debug_launch_set}
-------------------------------------------------

Several ways of setting code_saturne to run under a debugger are possible:

- Using the GUI, set options in
    <span style="color: rgb(0,91,187)">Run computation/Advanced options`</span>

  * The associated help provides several examples

  * This sets a <span style="color: rgb(0,91,187)">`debug_args`</span> option
    under the current resources section in the case's
     <span style="color: rgb(0,91,187)">`run.cfg`</span> file, which
    can also be edited directly.

- The same options can be provided directly to
  <span style="color: rgb(48,119,16)">`code_saturne run`</span>
  using the <span style="color: rgb(0,91,187)">`--debug-args`</span> option

\image html dg/debug_wrapper.png "Example of use of debugger wrapper" width=60%

When the code is run, the debugger will then be launched automatically.

Arguments required by the debug wrapper
---------------------------------------

To run the execution under a debugger, a string with the following syntax
structure should be used:
  ```
  <debugger> [debugger_options] [valgrind [[valgrind_options]]
  ```

Or, for Valgrind only:

  ```
  <valgrind> [valgrind options]
  ```

where `< >` denote required arguments, and `[ ]` optional arguments.


The following debuggers and user interfaces are handled:
**gdb** (GNU Debugger), **cuda-gdb**, **cgdb** (console-front-end to gdb),
**gdbgui** (browser-based frontend to gdb), **ddd** (Data Display Debugger),
**emacs** (as gdb front-end), **kdbg** (KDbg),
**kdevelop** (KDE developement environment),
**gede** (simple Qt-based gdb GUI), **nemiver** (GNOME Nemiver debugger),
**valgrind** (Valgrind tools for memory debugging).

If a debugger may not be found in the PATH, an absolute path should be given instead.

The following debug wrapper options are handled:

<table>
<tr><td> `--asan-bp`          <td> adds a breakpoint for gcc's Address-Sanitizer
<tr><td> `--back-end=gdb`     <td> path to debugger back-end (for graphical front-ends)
<tr><td> `--breakpoints=LIST` <td> comma-separated list of breakpoints to insert
<tr><td> `--ranks=LIST`       <td> comma-separated list of MPI ranks to debug
<tr><td> `--terminal`         <td> terminal type to use for console debugger
</table>

Other, standard options specific to each debugger may also be used, as long as they
do not conflict with options in this wrapper

To combine a Valgrind tool with another debugger, Valgrind's `--vgdb-error=<num>`
option should be used, where `<num>` is the number of errors after which Valgrind's
gdb server should be invoked. This mode is only compatible with the following
debuggers: gdb, gdbgui, ddd, emacs, as the other debugger front-ends do not provide
the required startup options to connect with the gdb server.

Note that compared to a standalone use of the gdb debugger, running under
the debug wrapper automatically sets breakpoints at program start and end
and launch the program.

To define commands such as setting breakpoints, a small gdb script with
the `.gdb` extension is generated for each rank. It can be removed safely
after program start.

When running in parallel, several debugging windows will be opened if necessary.

### Example commands:

To simply run using the gdb debugger:
```
gdb
```

To run Valgrind's gdb server, stopping at the second error in the terminal:
```
--terminal=gnome-terminal gdb valgrind --vgdb-error=1
```

To run Valgrind's gdb server, stopping at the second error:
```
gdb valgrind --vgdb-error=1
```

Do do the same using the ddd front-end:
```
ddd valgrind --vgdb-error=1
```

To run under gdb with preset brekpoints at bft_error and exit functions:
```
gdb --asan-bp --breakpoints=bft_error,exit
```

To debug under gdbgui:
```
gdbgui
```

To run under Valgrind's default tool (Memcheck), with a user Valgrind build
```
<path_to_valgrind> --tool=massif
```

To run under Valgrind's Massif heap profiler:
```
valgrind --tool=massif
```

### Terminal settings

To allow for debugging parallel runs and combining GDB and Valgrind, GDB is
run under a new terminal.

- The type of terminal chosen can by defined using the `--terminal` option.

  * Known issue: on some Debian 10-based systems, running under `gnome-terminal`
    crashes GDB. Running under the default `xterm` or `konsole` works fine.

By default, <span style="color: rgb(48,119,16)">`xterm`</span> will be used.
This usually leads to very small, hard to read fonts. This can be fixed
by editing <span style="color: rgb(48,119,16)">`$HOME/.Xresources`</span>
such as in the following example:

``` {style="VBScriptstyle"}
!xterm*font:     *-fixed-*-*-*-18-*
xterm*faceName: Liberation Mono:size=10:antialias=false
xterm*font: 7x13
xterm*VT100.geometry: 120x60
URxvt*geometry:  120x60
URxvt.font: xft:Terminus:antialias=false:size=10
```

Starting code_saturne under a debugger manually
-----------------------------------------------

Starting the debugger manually in an execution directory avoids creating
many directories and waiting for pre-processing before each run.

- <span style="color: rgb(48,119,16)">`cd`</span> to the run directory under
  <span style="color: rgb(0,91,187)">`RESU/<run_id>`</span>.

- To determine the code options already configured, run
  <span style="color: rgb(48,119,16)">`cat run_solver`</span> to view the
  execution commands.

- Add the debugger commands to this to run (unless already done
  through the GUI or user script).

  * To make this easier, code_saturne provides a
    <span style="color: rgb(0,91,187)">`cs_debug_wrapper.py`</span> script, in the
    <span style="color: rgb(0,91,187)">`python/code_saturne/base`</span> directory
    of the source tree (and in the
    <span style="color: rgb(0,91,187)">`lib/python<version>/site-packages/code_saturne`</span>
    directory of an installed build).

  * Run <span style="color: rgb(48,119,16)">`cs_debug_wrapper.py --help`</span> for
    instructions.

- The XML file may be modified directly using
  <span style="color: rgb(48,119,16)">`code_saturne gui <file>`</span>
  (ignoring the directory warning).

  * If mathematical expressions are modified, and additional step is required.
    in this case, it is simpler to generate a new run.

- When modifying user-defined functions, do not forget to run
  <span style="color: rgb(48,119,16)">`code_saturne compile -s src_saturne`</span>
  to update the <span style="color: rgb(0,91,187)">`cs_solver`</span> executable.

### Running under Vim, Neovim, or Atom

The code_saturne debug wrapper does not yet support launching GDB under Vim,
Neovim, or Atom.

Various examples of use of debugging with Vim are found on the web, explaining
how [Termdebug](https://www.dannyadam.com/blog/2019/05/debugging-in-vim/)
for example can be used.

Support for these tools could be added if users can provide a command-line
example for launching a debugging session with their favorite editor
and help test this.

Parallel Debugging
==================

Parallel Debugging: MPI
-----------------------

Debugging parallel code_saturne runs is not very different from
debugging serial runs.

- If a true parallel debugger such as
  <span style="color: rgb(48,119,16)">TotalView</span> or
  <span style="color: rgb(48,119,16)">Arm DDT</span>
  is available, do not hesitate to use it (by adapting the
  <span style="color: rgb(48,119,16)">`run_solver`</span> script in the
  exection directory), and ignore the rest of this slide.

- When no true parallel debugger is available, serial debuggers may be
  used.

  * Usually one for each process, though using multiple program
    features allows running only selected ranks under a debugger.

    - For example:
      <span style="color: rgb(48,119,16)">`mpiexec -n 2 <program> : - n 1 <debug_wrapper> <program> : -n 3 <program>`</span>
      to debug rank 2 of 6

  * The execution may not be restarted from the debugger; the whole
    parallel run must be restarted.

    - Very painful if not automated.

    - This is where the
      <span style="color: rgb(0,91,187)">`cs_debug_wrapper.py`</span> script
      really becomes useful.

      + This script also includes a `--ranks` filter option so as to call the
        debugger only on selected ranks. For example, using `--ranks=2,5`
        will launch MPI ranks 2 and 5 under a debugger, while other ranks will
        be run normally.

- For code_saturne under GDB, to determine a given process's rank, type:
  <span style="color: rgb(48,119,16)">`print cs_glob_rank_id`</span>

Parallel Debugging: OpenMP
--------------------------

Debugging OpenMP data races is much more tricky.

- Most errors are due to missing <span style="color: rgb(0,91,187)">`private`</span>
  attributes in OpenMP pragmas.

  * In C, using local variable declarations avoids most of these, as
    those variables are automatically thread-private.

  * Valgrind's <span style="color: rgb(48,119,16)">DRD</span> (Data Race Detector) tool
    is quite useful here.

    - <span style="color: rgb(48,119,16)">`valgrind –tool=drd –check-stack-var=yes –read-var-info=yes`</span>

    - Check
      [https://valgrind.org/docs/manual/drd-manual.html#drd-manual.openmp](https://valgrind.org/docs/manual/drd-manual.html#drd-manual.openmp>) for more information.

  * GCC's or clang's <span style="color: rgb(48,119,16)">ThreadSanitizer</span>
    is also very useful here.

- In both cases, to avoid false positives, GCC must be built with the
  <span style="color: rgb(48,119,16)">`–disable-linux-futex`</span>
  configure option, so this requires a special build of GCC.

  * With more recent versions of GCC, this may not be sufficient to
    avoid false positives...

    - probably due to some optimizations in thread management.

Miscellaneous tips
==================

When having or suspecting issues with loading or selection of shared libraries, use
```
export LD_DEBUG=libs
```
before running the code from a terminal. This will log many operations of
the dynamic loader.

To obtain more info on available options, use `LD_DEBUG=help` with any
program, for example:
```
LD_DEBUG=help cat
```

