%
% LaTeX class for Eureka documents
% Original author  : Wim van Ackooij
% edited by Brice Chardin
%

\ProvidesClass{edf}[2010/05/20]

% pass options to the article class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax

\LoadClass[a4paper,10pt]{report}

% \RequirePackage{uarial} % may be required for arial font: \fontfamily{ua1}
\RequirePackage{fancyhdr}
%\RequirePackage{cleveref} %ye renuda moving
\RequirePackage{xcolor}
\RequirePackage{array}
\RequirePackage{textcomp} % for copyright symbol in arial
\RequirePackage{xspace}
\RequirePackage{graphicx}
\RequirePackage{subfigure}
\RequirePackage{rotating} %ds adding
\RequirePackage{afterpage} %yf adding
\RequirePackage{multirow} %yf adding
\RequirePackage{float} %yf adding
\RequirePackage{lastpage} %yf adding
\RequirePackage[colorlinks=true,hyperfootnotes=false,linkcolor=black,citecolor=black,urlcolor=black]{hyperref}
\RequirePackage{amsmath}  %ye renuda adding
\RequirePackage{cleveref}  %ye renuda adding
\RequirePackage{siunitx}  %ye renuda adding
%\RequirePackage{totpages} %yf adding

%\RequirePackage{draftwatermark}
%\SetWatermarkLightness{0.8}
%\SetWatermarkAngle{25}
%\SetWatermarkScale{2}
%\SetWatermarkFontSize{2cm}
%\SetWatermarkText{   DRAFT}

%PAGE LAYOUT
%margins
\setlength{\voffset}{-1in} %top margin (page offset) = 0
\setlength{\topmargin}{1.25cm} %top margin = 1.25cm
\setlength{\hoffset}{-1in} %left margin (page offset) = 0
\setlength{\oddsidemargin}{2.3677cm} %odd left margin = 2.5cm -3.75pt
\setlength{\evensidemargin}{\oddsidemargin} %even left margin = odd left margin
%header
\setlength{\headheight}{2cm}
\setlength{\headsep}{0.5cm}
%footer
\setlength{\footskip}{1cm}
%text area : height
\setlength{\textheight}{\paperheight}  % should be 29.7cm
\addtolength{\textheight}{-5.55cm}
%text area : width
\setlength{\textwidth}{\paperwidth} % should be 21cm
\addtolength{\textwidth}{-4.815cm} %5cm -5.25pt
%no indentation, parskip 0.25cm
\parindent=0in
\setlength{\parskip}{0.25cm}

%EDF headers color
\definecolor{clsheaderblue}{RGB}{31,61,114}

\fancyhf{} % deletes original headers and footers
\renewcommand{\headrulewidth}{0pt}
\newcommand{\pageend}{\pageref*{LastPage}}
\pagestyle{fancyplain}

% Header lengths
\newlength{\clsgauche}
\setlength{\clsgauche}{1.75cm}
\addtolength{\clsgauche}{-4pt}
\newlength{\clsdroite}
\setlength{\clsdroite}{3.5cm}
\addtolength{\clsdroite}{-3.25pt}
\newlength{\clsmilieu}
\setlength{\clsmilieu}{\textwidth}
\addtolength{\clsmilieu}{-\clsgauche}
\addtolength{\clsmilieu}{-\clsdroite}
\addtolength{\clsmilieu}{-2pt}
% header vertical span
\newlength{\clsblankheaderspace}
\setlength{\clsblankheaderspace}{-2cm}
\addtolength{\clsblankheaderspace}{11.5pt}
\newlength{\clsblankheaderskip}
\setlength{\clsblankheaderskip}{2cm}
\addtolength{\clsblankheaderskip}{-1.5pt}

% Footer lengths ; these lengths are not exactly the same as in the .doc template
\newlength{\clsgauchef}
\setlength{\clsgauchef}{5.25cm}
\newlength{\clsdroitef}
\setlength{\clsdroitef}{\clsgauchef}
\newlength{\clsmilieuf}
\setlength{\clsmilieuf}{\textwidth}
\addtolength{\clsmilieuf}{-\clsgauchef}
\addtolength{\clsmilieuf}{-\clsdroitef}
\addtolength{\clsmilieuf}{-2pt}

\newcommand{\header}[3]{%
  \fancyhead[C]{%
    \hfuzz=3pt % disable warnings for this header, which has a 3pt overfull hbox
    \begin{tabular}{!{\color{clsheaderblue}\vrule width 0.75pt}@{}p{\clsgauche}@{}
                    !{\color{clsheaderblue}\vrule width 0.75pt}
                    @{\hspace{1pt}}p{\clsmilieu}@{\hspace{1pt}}
                    !{\color{clsheaderblue}\vrule width 0.75pt}@{}p{\clsdroite}@{}
                    !{\color{clsheaderblue}\vrule width 0.75pt}@{}p{0pt}@{}}
      \noalign{{\color{clsheaderblue}\hrule height 0.75pt}}
      \parbox[t]{\clsgauche}{
        \fontsize{8}{10}
        \fontfamily{ua1}\selectfont
        \raggedright
        \hspace{1pt}
        \textbf{EDF R\&D}
      }
      &
      \vskip-10pt
      \parbox[t]{\clsmilieu}{
        \centering
        \fontsize{6}{9}
        \fontfamily{ua1}\selectfont
        #1
      }
      \vskip10pt
      &
      \parbox[t]{\clsdroite}{
        \fontsize{8}{10}
        \fontfamily{ua1}\selectfont
        \raggedleft
        \textbf{#2} \hspace*{1.75pt} \\[6pt]
        \textbf{Version #3} \hspace*{1.75pt}
      }
      &
      \rule[\clsblankheaderspace]{0pt}{\clsblankheaderskip}
      \cr
      \noalign{{\color{clsheaderblue}\hrule height 0.75pt}}
    \end{tabular}\@gobble % disable fancyhdr automatic header resizing
  }
}

\newcommand{\footer}[2]{%
  \fancyfoot[C]{%
    \hfuzz=3pt % disable warnings for this footer, which has a 3pt overfull hbox
    \begin{tabular}{!{\color{clsheaderblue}\vrule width 0.75pt}@{}p{\clsgauchef}@{}
                    !{\color{clsheaderblue}\vrule width 0.75pt}@{\hspace{1pt}}p{\clsmilieuf}@{\hspace{1pt}}
                    !{\color{clsheaderblue}\vrule width 0.75pt}@{}p{\clsdroitef}@{}
                    !{\color{clsheaderblue}\vrule width 0.75pt}}
      \noalign{{\color{clsheaderblue}\hrule height 0.75pt}}
      \parbox[t]{\clsgauchef}{
        \raggedright
        \fontsize{9}{9}
        \fontfamily{ua1}\selectfont
        \hspace*{1pt}
        Accessibility : \textbf{#1}
      }
      &
      \parbox[t]{\clsmilieuf}{
        \centering
        \fontsize{9}{9}
        \fontfamily{ua1}\selectfont
        Page \thepage\ on \pageend
      }
      &
      \parbox[t]{\clsdroitef}{
        \raggedleft
        \fontsize{9}{9}
        \fontfamily{ua1}\selectfont
        \copyright EDF SA #2 \hspace*{1pt}
      }
      \cr
      \noalign{{\color{clsheaderblue}\hrule height 0.75pt}}
    \end{tabular}
  }
}

\renewenvironment{thebibliography}[1]
     {\list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}

% Condensed TOC
\renewcommand\tableofcontents{%
  \renewcommand{\contentsname}{Sommaire / Summary}
  {%
    \newpage
    \setlength{\parskip}{0pt}
%    \begin{center}
      \section*{\contentsname
        \@mkboth{%
          \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
%    \end{center}
    \@starttoc{toc}%
    \newpage
  }
}

% default header/footer values, should be redefined
\header{Title}{H-000-0000-00000-FR}{0.0}
\footer{EDF}{\textcolor{red}{En projet}}

\newcommand{\caution}{%
\newpage
\begin{center}
\section*{AVERTISSEMENT / CAUTION}
\end{center}
\addcontentsline{toc}{section}{AVERTISSEMENT / CAUTION}

\bigskip

L'acc\`{e}s \`{a} ce document, ainsi que son utilisation, sont strictement limit\'{e}s aux personnes express\'{e}ment habilit\'{e}es par EDF.

EDF ne pourra \^{e}tre tenu responsable, au titre d'une action en responsabilit\'{e} contractuelle, en responsabilit\'{e} d\'{e}lictuelle ou de tout autre action, de tout dommage direct ou indirect, ou de quelque nature qu'il soit, ou de tout pr\'{e}judice, notamment, de nature financier ou commercial, r\'{e}sultant de l'utilisation d'une quelconque information contenue dans ce document.

Les donn\'{e}es et informations contenues dans ce document sont fournies "en l'\'{e}tat" sans aucune garantie expresse ou tacite de quelque nature que ce soit.

Toute modification, reproduction, extraction d'\'{e}l\'{e}ments, r\'{e}utilisation de tout ou partie de ce document sans autorisation pr\'{e}alable \'{e}crite d'EDF ainsi que toute diffusion externe \`{a} EDF du pr\'{e}sent document ou des informations qu'il contient est strictement interdite sous peine de sanctions.

\bigskip

\begin{center}
\--\--\--\--\--\--\--
\end{center}

The access to this document and its use are strictly limited to the persons expressly authorized to do so by EDF.

EDF shall not be deemed liable as a consequence of any action, for any direct or indirect damage, including, among others, commercial or financial loss arising from the use of any information contained in this document.

This document and the information contained therein are provided "as are" without any warranty of any kind, either expressed or implied.

Any total or partial modification, reproduction, new use, distribution or extraction of elements of this document or its content, without the express and prior written consent of EDF is strictly forbidden. Failure to comply to the above provisions will expose to sanctions.}

\newcommand{\synthese}[1]{%
\newpage
\begin{center}
\section*{Synth\`{e}se}
\end{center}
\addcontentsline{toc}{section}{Synth\`{e}se}
#1}

\newcommand{\summary}[1]{%
\newpage
\begin{center}
\section*{Executive Summary}
\end{center}
\addcontentsline{toc}{section}{Executive Summary}
#1}

