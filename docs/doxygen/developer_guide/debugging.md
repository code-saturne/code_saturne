% All settings, including handout. See style directory.
\input{cs_beamer_settings}
\geometry{paperwidth=140mm,paperheight=105mm}

%%%%%%%%%%%%%%%%%%%%%%%
% Background image if needed
%%%%%%%%%%%%%%%%%%%%%%%
\setbeamertemplate{background canvas}{
\includegraphics[width=\paperwidth]{bg_title_EDF_green.pdf}
}

%%%%%%%%%%%%%%%%%%%%%%%
% To debug
%%%%%%%%%%%%%%%%%%%%%%%
%\includeonlyframes{current}

\title[\CS: advanced developer training]{
\vspace{1cm}\\
\Huge \textbf{\CS advanced user and developer training: debugging}
}

\author[\bf{\CS~dev. team}]{
 \textcolor{white}{\mbox{\CS~development team
 %\inst{1}
 }}
}

\date{
  \flushleft \textcolor{white}{\tiny \the\year}
\vspace{-1cm}\\
\flushright \includegraphics[width=0.15\textwidth]{LOGO-blanc-HD.png}
}

\institute{
\textcolor{white}{\inst{1}Fluid Mechanics, Energy and Environment, EDF R\&D, Chatou, France,\\
}}


%%%%%%%%%%%%
% LOGO TITLE PAGE
%%%%%%%%%%%%
%\logo{}%\includegraphics[width=0.1\textwidth]{Docs/EDF_Logo_PMS_v_F.pdf}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO PUT THE TABLE OF CONTENTS
% AT THE BEGINNING OF EACH SECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\AtBeginSection[] %
{ \begin{frame}<beamer>
\frametitle{Overview}
\begin{small}
\tableofcontents[currentsection]
\end{small}
\end{frame} }

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BEGIN DOC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{document}

%%%%%%%%%%%%%
%TITLE
%%%%%%%%%%%%%
\begin{frame}[plain]
	\titlepage
\end{frame}

%Cancel the logos for the following and the background
\logo{}
\setbeamertemplate{background canvas}{}

\mainfootline

\lstdefinestyle{LFrame}{frame=single,frameround=tttf,xrightmargin=3pt,fillcolor=\color{whitesmoke},backgroundcolor=\color{whitesmoke},rulecolor=\color{greenedf}}

\lstdefinestyle{FFrame}{frame=single,frameround=tttf,xrightmargin=3pt,fillcolor=\color{whitesmoke},backgroundcolor=\color{whitesmoke},rulecolor=\color{blueedf}}

\lstdefinestyle{Cstyle}{language=C,basicstyle=\footnotesize\ttfamily,keywordstyle=\color{blue}\bfseries,commentstyle=\color{red},stringstyle=\color{greenedf}\ttfamily,showstringspaces=false,labelstyle=\tiny,escapeinside={!@}{@!}}

\lstdefinestyle{Fstyle}{language={[95]Fortran},basicstyle=\footnotesize\ttfamily,keywordstyle=\color{blue}\bfseries,commentstyle=\color{red},stringstyle=\color{greenedf}\ttfamily,showstringspaces=false,labelstyle=\tiny,escapeinside={!@}{@!}}

\lstdefinestyle{VBScriptstyle}{language=VBScript, basicstyle=\footnotesize\ttfamily, keywordstyle=\color{blueedf}\ttfamily, stringstyle=\color{orangededf}\ttfamily, commentstyle=\color{greenedf}\ttfamily, morecomment=[l][\color{magenta}]{\#}}

\def\inlinec{\lstinline[style=Cstyle]}

%\lstset{Cstyle}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\begin{frame}
\frametitle{Overview}
\begin{small}
\tableofcontents
\end{small}
\end{frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toolchain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\section{Introduction}

\begin{frame}[t]{\CS Debugging}

All non-trivial software has bugs, at least at first,
so debugging is a part of code development.

\begin{block}{}

\begin{itemize}
  \blueitem<+-> Is is often easier to debug newly written code than older code
  \begin{itemize}
    \item At least when done by the same person
    \item Test and debug early, before you forget the details of your code
  \end{itemize}
  \blueitem<+-> Multiple tools and techniques are available
  \begin{itemize}
    \item Mastering debugging tools and techniques make debugging less painful
  \end{itemize}
  \blueitem<+-> Bugs detected late are much more costly then those detected early
  \begin{itemize}
    \item May require re-validation
    \item Often harder to debug, because the code details need to be ``re-learned''
    \item Be proactive, try to detect as many bugs as possible by testing
  \end{itemize}
\end{itemize}

\end{block}

\end{frame}

\section{Tools}

\subsection{Debugging methods}

\begin{frame}[t]{Debugging methods}

When encountering or suspecting a bug, choosing the best debugging technique for
the situation can lead to one or more orders of magnitude in time savings.

\begin{block}{}

\begin{itemize}
  \blueitem<+-> Choosing which tool to use is the difficult part, and
  can be based on several factors:
  \begin{itemize}
    \item<+-> Comparison of experience to observed effects
    \item<+-> Leveraging theoretical knowledge and experience
    \begin{itemize}
      \item<+-> Guessing at whether bugs may be due to uninitialized values,
        out-of-bounds arrays accesses, bad option settings, numerical errors, ...
    \end{itemize}
    \item<+-> Sometimes, a bit of luck...
  \end{itemize}
  \blueitem<+-> \CS tries to help, so always check for \textcolor{orangeedf}{error messages}
  \begin{itemize}
    \item<+-> In \CS, \textcolor{greenedf}{error*}, \textcolor{greenedf}{listing},
      and error logs in batch run logs (or in the console) should be checked
    \begin{itemize}
      \item<+-> For parallel runs, when both \textcolor{greenedf}{error*} and \textcolor{greenedf}{error\_r*} are present, the \textcolor{greenedf}{error\_r*} files are the ones which contain the useful information
    \end{itemize}
    \item<+-> Some graphical checks with
      \textcolor{greenedf}{postprocessing/error.*} are also available for
      boundary conditions and linear solvers.
  \end{itemize}
\end{itemize}

\end{block}

\end{frame}

\begin{frame}[t]{}

Some debugging tools that should be considered

\begin{block}{}

\begin{itemize}
   \blueitem<+-> Source code proofreading
  \begin{itemize}
    \item \textcolor{orangeedf}{Re}-checking for compiler warnings
  \end{itemize}
  \blueitem<+-> Interactive debugging
  \blueitem<+-> Memory debugging
  \blueitem<+-> Checks/instrumentation in code...
  \begin{itemize}
    \item Use \alerttt{--enable-debug} to configure builds
      for debug
    \begin{itemize}
      \item Enables use of many \alerttt{assert} checks in C code
      \item Enables arrays bounds-checking in Fortran
    \end{itemize}
    \item Using recent versions of the \toolname{GCC} or \toolname{clang} compiler, compile/run with
      \toolname{AddressSanitizer}, \toolname{UndefinedBehaviorSanitizer}, and \toolname{ThreadSanitizer} occasionally
    \begin{itemize}
      \item Code built this way not compatible with runs under \toolname{Valgrind}
      \item Not compatible either with some resource limits set on some clusters
      \item Overhead: usually x3
    \end{itemize}
  \end{itemize}
  \blueitem<+-> When you known where to search, \textcolor{orangeedf}{print} statements may be useful...
\end{itemize}

\end{block}

\end{frame}

\subsection{The GNU debugger}

\begin{frame}[t]{The GNU debugger}

The GNU debugger \textcolor{blueedf}{\url{http://www.gnu.org/software/gdb}} is a broadly available, interactive debugger for compiled languages including C, C++, and Fortran.

\begin{block}{}

\begin{itemize}
  \blueitem<+-> To debug an executable, run \alerttt{gdb <executable>}
  \begin{itemize}
    \item Under the gdb prompt, type \alerttt{help} for built-in help, \alerttt{q} to quit
    \begin{itemize}
      \item Help is grouped in categories
      \item The most common options are \alerttt{b} (set breakpoint), \alerttt{n} (continue to next statement), \alerttt{c} (continue to next breakpoint), and \alerttt{p} (print)
    \end{itemize}
  \end{itemize}
  \blueitem<+-> Many front-ends are available, such as \toolname{Emacs}, \toolname{vim} (3 modes), and the graphical \toolname{ddd}, \toolname{Nemiver}, \toolname{Kdevelop}, an \toolname{gdbgui}
  \blueitem<+-> GDB can provide some information on any program, but can provide more detailed and useful information when the program was compiled with debugging info. The matching compiler option is usually \textcolor{orangeedf}{-g}, and in the case of \CS, is provided using the \alerttt{--enable-debug} configure option at installation
\end{itemize}

\end{block}

\end{frame}

\begin{frame}[t]{}

When used directly, GDB runs in a single terminal frame, as shown here.
Only the current line of code is shown, though the \alerttt{list} command allows showing more.
\vskip5pt
\centering
\includegraphics[width=0.8\textwidth]{graphics/gdb_screen.png}

\end{frame}

\begin{frame}[t]{}

When started with the \alerttt{--tui} option, GDB runs in a split terminal, with source on top, commands on bottom. Using the \textcolor{orangeedf}{CTRL+x+o} key combination allows changing focus from one to the other. Using \textcolor{orangeedf}{CTRL+l} allows refreshing the display.

\vskip5pt
\centering
\includegraphics[width=0.8\textwidth]{graphics/gdb_tui_screen.png}

\end{frame}

\begin{frame}[t]{}

GDB may also be run under Emacs, which provides syntax highlighting of source code.

\vskip5pt
\centering
\includegraphics[width=0.65\textwidth]{graphics/emacs_gud_screen.png}

\end{frame}

\begin{frame}[t]{}

The DDD (Data Display Debugger) front-end uses a dated graphical toolkit, but has the advantage of combining a command prompt with graphical tools, and is very easy to use.

\vskip5pt
\centering
\includegraphics[width=0.6\textwidth]{graphics/ddd_screen.png}

\end{frame}

\begin{frame}[t]{}

The GNOME Nemiver debugger also has a GDB backend. This debugger offers a clean display, but lacks the possibility of typing commands; everything must be done using the mouse and menus, which is often tedious.

\vskip5pt
\centering
\includegraphics[width=0.7\textwidth]{graphics/nemiver_screen.png}

\end{frame}

\subsection{The GNU Valgrind tool suite}

\begin{frame}[t]{The Valgrind tool suite}

The Valgrind tools \textcolor{blueedf}{\url{http://www.valgrind.org}} allows the
detection of many memory management (and other) bugs.

\begin{block}{}

\begin{itemize}
  \blueitem<+-> Dynamic instrumentation
  \begin{itemize}
    \item No need for recompilation
    \begin{itemize}
      \item Provides more info (i.e. code line numbers) with code compiled in
        debug mode
    \end{itemize}
    \item Depending on tool used, run time and memory overhead from 10-100x
    \begin{itemize}
      \item With default tool (\textcolor{blueedf}{Memcheck}), 10x30
      \item Use proactively, to detect bugs on small cases, before they become a
        problem in production cases
    \end{itemize}
  \end{itemize}
  \blueitem<+-> May be combined with an interactive debugger
\end{itemize}

\end{block}

\end{frame}

\begin{frame}[t]{}

Valgrind is easy to run:

\begin{block}{}

\begin{itemize}
  \blueitem<+-> Prefix a standard command with \alerttt{valgrind}
  \begin{itemize}
    \item By default, uses the \textcolor{orangeedf}{\texttt{memcheck}} tool
    \item Tool may be changed using \alerttt{valgrind --tool=/cachegrind/callgrind/drd/massif/...}
  \end{itemize}
  \blueitem<+-> Valgrind can be combined with the \textcolor{orangeedf}{\texttt{gdb}} debugger, using its \textcolor{orangeedf}{\texttt{gdbserver}} mode
  \begin{itemize}
    \item To use this mode, call \alerttt{valgrind --vgdb-error=<number>} where the number represents the number of errors after which the gdbserver is invoked (0 to start immediatedly)
  \end{itemize}
\end{itemize}

\end{block}

\end{frame}

\begin{frame}[t]{}

\vskip5pt
\centering
\includegraphics[width=0.8\textwidth]{graphics/valgrind_screen.png}

\end{frame}

\begin{frame}[t]{GCC and clang sanitizers}

Recent versions of the LLVM clang and GCC compilers have additional instrumentation options, allowing memory debugging with a lower overhead than Valgrind.

\begin{block}{}

\begin{itemize}
  \blueitem<+-> For the most common errors, use \toolname{AddressSanitizer}, a fast memory error detector
  \begin{itemize}
    \item For the \CS configure options, this means \alerttt{CFLAGS=-fsanitize=address}, \alerttt{FCFLAGS=-fsanitize=address}, and \alerttt{LDFLAGS=-fsanitize=address}
    \item This may sometimes require specifying \alerttt{export LD\_LIBRARY\_FLAGS=<path\_to\_compiler\_libraries} when the compiler is installed in a nonstandard path on older systems
    \item On some machines, this may be unusable if memory resource limits are set (check using \alerttt{ulimit -c})
    \item Note that the resulting code will not be usable under Valgrind
    \item Uninitialized values are not detected by Address Sanitizer (but may be detected by UndefinedBehaviorSanitizer)
    \item Out-of-bounds errors for arrays on stack (fixed size, usually small) are not detected by Valgrind, but may be detected by AddressSanitizer
  \end{itemize}
\end{itemize}

\end{block}

export ASAN_OPTIONS=detect_leaks=0

\end{frame}

\begin{frame}[t]{}

\begin{block}{}

\begin{itemize}
  \blueitem<+-> The \toolname{UndefinedBehaviorSanitizer} instrumentation is also useful to detect other types of bugs, such as division by zero, some memory errors, integer overflows, and more
  \begin{itemize}
    \item This may sometimes require also specifying \alerttt{-lubsan} and even in some cases specify \alerttt{LD\_LIBRARY\_FLAGS} as above
    \item For the \CS configure options, this means \alerttt{CFLAGS=-fsanitize=undefined}, \alerttt{FCFLAGS=-fsanitize=undefined}, and \alerttt{LDFLAGS=-fsanitize=undefined}
    \item This may sometimes require specifying \alerttt{LD\_LIBRARY\_FLAGS} as per AddressSanitizer
  \end{itemize}
  \blueitem<+-> Note that only code compiled with those options is instrumented
  \blueitem<+-> The \toolname{ThreadSanitizer} is useful to detect OpenMP errors, but requires a specific build of GCC to avoid false errors; its use is detailed in the ``Parallel Debugging'' section
\end{itemize}

\end{block}

\end{frame}

\section{Application to \CS}

\begin{frame}[t]{Starting \CS under a debugger}

Several ways of running \CS under a debugger are possible:

\begin{itemize}
  \blueitem<+-> Using the GUI or the \exmplett{domain.debug} setting in \exmplett{cs\_user\_scripts.py} to automatically run the code under a debugger
  \begin{itemize}
    \item Set options in \exmplett{Run computation/Advanced options}
    \item As for regular runs, this will create a new directory under \exmplett{RESU} for each run and test
  \end{itemize}
  \blueitem<+-> Preparing a run directory using \alerttt{code\_saturne~run~[options]~--initialize}, then running the debugger manually from the run directory
  \begin{itemize}
    \item If the code has crashed during a previous run, this is not necessary, as the matching run directory remains in a initialized state
  \end{itemize}
  \blueitem<+-> Combining both approaches:
  \begin{itemize}
    \item Prepare a first run using the GUI or user script to handle the debugger syntax, then (re-)run the debugger manually
  \end{itemize}
\end{itemize}

\end{frame}

\begin{frame}{Example of use of debugger wrapper}
  \centering
  \includegraphics[width=1.0\textwidth]{graphics/debug_wrapper.png}
\end{frame}

\begin{frame}[fragile]{A XTerm configuration example}{file \alerttt{.Xresources} on your home}
  %
  \begin{lstlisting}[style=VBScriptstyle]
!xterm*font:     *-fixed-*-*-*-18-*
xterm*faceName: Liberation Mono:size=10:antialias=false
xterm*font: 7x13
xterm*VT100.geometry: 120x60
URxvt*geometry:  120x60
URxvt.font: xft:Terminus:antialias=false:size=10
  \end{lstlisting}
\end{frame}


\begin{frame}[t]{Starting \CS under a debugger manually}

Starting the debugger manually in an execution directory avoids creating many directories and waiting for preprocessing before each run.

\begin{itemize}
  \blueitem<+-> \alerttt{cd} to the run directory under \exmplett{RESU/<run\_id>}
  \item To determine the code options already configured, run \alerttt{cat~run\_solver} To view the execution commands
  \blueitem<+-> Add the debugger commands to this to run (unless already done through the GUI or user script)
  \begin{itemize}
    \item To make this easier, \CS provides a \exmplett{cs\_debug\_wrapper.py} script, in the \exmplett{bin} directory of the source tree (and in the \exmplett{lib/python<version>/site-packages/code\_saturne} directory of an installed build)
    \item Run  \alerttt{cs\_debug\_wrapper.py~--help} for instructions
  \end{itemize}
  \blueitem<+-> The XML file may be modified directly using \alerttt{code\_saturne gui <file>} (ignoring the directory warning)
  \blueitem<+-> When modifying user-defined functions, do not forget to run \alerttt{code\_saturne~compile~-s~src\_saturne} to update the \exmplett{cs\_solver} executable
\end{itemize}

\end{frame}

\section{Parallel Debugging}

\begin{frame}[t]{Parallel Debugging: MPI}

\begin{block}{}

Debugging parallel \CS runs is not very different from debugging serial runs.

\begin{itemize}
  \blueitem<+-> If a true parallel debugger such as \toolname{TotalView} or \toolname{Arm DDT} is available, do not hesitate to use it, and ignore the rest of this slide
  \blueitem<+-> When no true parallel debugger is available, serial debuggers may be used
  \begin{itemize}
    \item Usually one for each process, though using multiple program features allows running only selected ranks under a debugger
    \begin{itemize}
      \item For example: \alerttt{mpiexec -n 2 <program> : - n 1 <debug\_wrapper> <program> : -n 3 <program>} to debug rank 2 of 6
    \end{itemize}
    \item The execution may not be restarted from the debugger; the whole parallel run must be restarted
    \begin{itemize}
      \item Very painful if not automated
      \item This is where the \exmplett{cs\_debug\_wrapper.py} script really becomes useful
    \end{itemize}
  \end{itemize}
  \blueitem<+-> For \CS under \toolname{GDB}, to determine a given process's rank, type: \alerttt{print cs\_glob\_rank\_id}
\end{itemize}

\end{block}

\end{frame}

\begin{frame}[t]{Parallel Debugging: OpenMP}

\begin{block}{}

Debugging OpenMP data races is much more tricky.

\begin{itemize}
  \blueitem<+-> Most errors are due to missing \exmplett{private} attributes in OpenMP pragmas
  \begin{itemize}
    \item In C, using local variable declarations avoids most of these, as those variables are automatically thread-private
  \end{itemize}
  \blueitem<+-> Valgrind's \toolname{DRD} (Data Race Detector) tool is quite useful here
  \begin{itemize}
     \item \alerttt{valgrind --tool=drd --check-stack-var=yes --read-var-info=yes}
     \item Check \textcolor{blueedf}{\url{http://valgrind.org/docs/manual/drd-manual.html\#drd-manual.openmp}} for more information
  \end{itemize}
  \blueitem<+-> GCC's \toolname{ThreadSanitizer} is also very useful here
  \blueitem<+-> In both cases, to avoid false positives, GCC must be built with the \alerttt{--disable-linux-futex} configure option, so this requires a special build of GCC
  \begin{itemize}
     \item With more recent versions of GCC, this may not be sufficient to avoid false positives...
  \end{itemize}
\end{itemize}

\end{block}

\end{frame}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\section[]{Questions}

\begin{frame}[label=end]
\begin{center}
\begin{Huge}
\alert{
\Huge Thank you for your attention...

Any question?
}
\end{Huge}
\end{center}
\end{frame}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The end of the document
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\end{document}
%
