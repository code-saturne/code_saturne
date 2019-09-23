# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.studymanager.cs_studymanager_run import run_studymanager_command

#-------------------------------------------------------------------------------

class TexWriter(object):
    """
    """
    def __init__(self, dest, filename, log, pdflatex):
        self.__dest = dest
        self.__filename = os.path.join(self.__dest, filename)
        self.__doc = []
        self.__log = log
        self.__pdflatex = pdflatex


    def rawLine(self, line):
        self.__doc.append(line)


    def appendLine(self, line):
        line = line.replace("_", "\_")
        self.__doc.append("%s \n" % (line))


    def addFigure(self, g):
        self.__doc.append("\\begin{center}\n")
        self.__doc.append("\\includegraphics[width=0.99\\textwidth]{%s}\n" % g)
        self.__doc.append("\\end{center}\n")


    def addInput(self, filename):
        self.appendLine("\n\\begin{verbatim}")
        f = open(filename)
        self.rawLine(f.read())
        f.close()
        self.appendLine("\\end{verbatim}\n")


    def addTexInput(self, filename):
        f = open(filename)
        self.rawLine(f.read())
        f.close()


    def tabCreate(self, columns):
        assert type(columns) == list

        self.__doc.append("\\begin{center}\n")
        s = "l|"*len(columns)
        self.__doc.append("\\begin{longtable}{|%s}\n" % s)
        self.__doc.append("\hline\n")

        for i in range(len(columns)):
            if i != len(columns)-1:
                self.__doc.append("\\textbf{%s} &" % (columns[i]))
            else:
                self.__doc.append("\\textbf{%s} \\\ \n" % (columns[i]))

        self.__doc.append("\hline\n")
        self.__doc.append("\hline\n")


    def tabWrite(self, columns):
        assert type(columns) == list

        for i in range(len(columns)):
            if columns[i] == "OK":
                columns[i] = "\\textcolor{green}{OK}"
            elif columns[i] == "KO":
                columns[i] = "\\textcolor{red}{KO}"
            elif columns[i] == None:
                columns[i] = "\\textit{default}"

            if i != len(columns)-1:
                self.__doc.append("%s &" % (columns[i]))
            else:
                self.__doc.append("%s \\\ \n" % (columns[i]))

        self.__doc.append("\hline\n")


    def tabClose(self):
        self.__doc.append("\\end{longtable}\n")
        self.__doc.append("\\end{center}\n")


    def write(self):
        # header
        head = []
        head.append("\\documentclass[a4paper, 10pt]{article}\n")
        head.append("\\usepackage[latin1]{inputenc}\n")
        head.append("\\usepackage[T1]{fontenc}\n")
        head.append("\\usepackage[normalem]{ulem}\n")
        head.append("\\usepackage[french]{babel}\n")
        head.append("\\usepackage{verbatim}\n")
        head.append("\\usepackage{color}\n")
        head.append("\\usepackage{longtable}\n")
        head.append("\\usepackage{listings}\n")
        head.append("\\usepackage{ifpdf}\n")
        head.append("\\usepackage{tikz}\n")
        head.append("\\usepackage{pgfplots}\n")
        head.append("\\pgfplotsset{\n")
        head.append("compat=newest,\n")
        head.append("xlabel near ticks,\n")
        head.append("ylabel near ticks\n")
        head.append("}\n")
        head.append("\\setlength{\\voffset}{0pt}")
        head.append("\\setlength{\\topmargin}{0pt}")
        head.append("\\addtolength{\\topmargin}{-13mm}")
        head.append("\\setlength{\\headheight}{15mm}")
        head.append("\\setlength{\\headsep}{6mm}")
        head.append("\\setlength{\\textheight}{233mm}")
        head.append("\\setlength{\\footskip}{15mm}")
        head.append("\\setlength{\\hoffset}{0pt}")
        head.append("\\setlength{\\evensidemargin}{0mm}")
        head.append("\\setlength{\\oddsidemargin}{0mm}")
        head.append("\\setlength{\\textwidth}{162mm}")
        head.append("\\setlength{\\parindent}{0mm}")
        head.append("\\setlength{\\parskip}{6pt}")
        head.append("\\setlength{\\tabcolsep}{1mm}")


        head.append("\\begin{document}\n")

        # end
        tail = []
        tail.append("\\end{document} \n")

        # write
        f = open(self.__filename + ".tex", mode='w')
        f.writelines(head + self.__doc + tail)
        f.close()


    def make_pdf(self):
        """
        Buld the pdf file, and clean the temporary files.
        """
        cmd = "pdflatex " + self.__filename + ".tex"
        r, t = run_studymanager_command(cmd, self.__log)

        for suffixe in ["tex", "log", "aux"]:
            f = self.__filename + "." + suffixe
            if os.path.isfile(f):
                os.remove(self.__filename + "." + suffixe)

        return self.__filename + ".pdf"


    def tex_finalize(self):
        """
        Finalize by making pdf if pdflatex is enabled,
        returns pdf file name (or tex file name if pdflatex disabled).
        """
        if self.__pdflatex:
            filename = self.make_pdf()
        else:
            filename = self.__filename + ".tex"

        return filename

#-------------------------------------------------------------------------------

class Report1(TexWriter):
    """
    Global report.
    """
    def __init__(self, dest, label, log, report, xml, pdflatex):
        TexWriter.__init__(self, dest, label, log, pdflatex)
        self.appendLine("\\section{Summary}")
        self.tabCreate(["Study / Case", "Compilation", "Run", "Time (s)", "Difference"])
        self.xml = xml
        self.report = report


    def add_row(self, studyLabel, caseLabel, is_compil, is_run, is_time, is_compare, is_diff):
        if is_compare == "not done":
            threshold = "Not used"
            is_diff   = "Not used"

        label = "%s / %s" % (studyLabel, caseLabel)
        label = label.replace("_", "\_")
        self.tabWrite([label, is_compil, is_run, is_time, is_diff])


    def close(self):
        self.tabClose()

        self.appendLine("\\section{Log}")
        self.appendLine("\\tiny\n\\begin{verbatim}")
        f = open(self.report)
        self.rawLine(f.read())
        f.close()
        self.appendLine("\\end{verbatim}\n\\normalsize")

        self.appendLine("\\section{File of commands}")
        self.appendLine("\\tiny\n\\begin{verbatim}")
        self.rawLine(self.xml)
        self.appendLine("")
        self.appendLine("\\end{verbatim}\n\\normalsize")

        self.write()

        return self.tex_finalize()

#-------------------------------------------------------------------------------

class Report2(TexWriter):
    """
    Detailed report.
    """
    def __init__(self, dest, label, log, pdflatex):
        """
        """
        TexWriter.__init__(self, dest, label, log, pdflatex)


    def add_row(self, values, studyLabel, caseLabel):
        nbvalue = len(values)
        row_max = 40

        if nbvalue:
            self.tabCreate(["Variable Name", "Diff. Max", "Diff. Mean", "Threshold"])
            for j in range(len(values)):
                self.tabWrite([values[j][0], values[j][1], values[j][2], values[j][3]])

            self.tabClose()
            self.appendLine("\n \\newpage \n")


    def close(self):
        self.write()
        return self.tex_finalize()

#-------------------------------------------------------------------------------

def test():
    """
    Test function.
    The pupose is to build 2 documents:
    1. a global report for all cases:
    2. a detailled report for each case:
    """

    dest      = os.getcwd()
    doc       = Report1(dest, "r1")
    compil    = True
    run       = True
    compare   = True
    diff      = True
    repo      = True
    threshold = 1.e-15
    doc.add_row("MYSTUDY", "MYCASE", compil, run, compare, diff, repo, threshold)
    doc.close()

    msg = """
  .----------------------------.
  |   Code_Saturne file dump   |
  `----------------------------'

Opening input file: "RESU/20110217-2231/checkpoint/main"

  File type: Checkpoint / restart, R0

  Base header size: 128
  Header alignment: 64
  Body alignment:   64

Opening input file: "RESU/20110217-2233/checkpoint/main"

  File type: Checkpoint / restart, R0

  Base header size: 128
  Header alignment: 64
  Body alignment:   64

  "nbre_pas_de_temps               "; Location:  0; Type: i4    ; Size: 1
    Differences: 1

  "instant_precedent               "; Location:  0; Type: r8    ; Size: 1
    Differences: 1; Max: 200; Mean: 200

  "pression_ce_phase01             "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14574; Max: 3.09509; Mean: 0.159934

  "vitesse_u_ce_phase01            "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14574; Max: 2.00929; Mean: 0.118064

  "vitesse_v_ce_phase01            "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14574; Max: 1.84415; Mean: 0.065676

  "vitesse_w_ce_phase01            "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14574; Max: 1.15005; Mean: 0.093984

  "k_ce_phase01                    "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14574; Max: 2.49588; Mean: 0.0191537

  "eps_ce_phase01                  "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14574; Max: 7438.13; Mean: 12.6477

  "scalaire_ce_0001                "; Location:  1; Type: r8    ; Size: 14574
    Differences: 14500; Max: 0.11862; Mean: 0.00345369

    """

    doc  = Report2(dest, "r2")
    v    = []
    info = msg.replace("\"", " ").replace(";", " ").replace(":", " ").split()
    print(info)
    repo = "RESU/20110217-2231/checkpoint/main"
    dest = "RESU/20110217-2233/checkpoint/main"

    for i in range(len(info)):
        if info[i][:4] == 'Diff':
            if info[i-3] not in ['i4', 'u4']:
                v.append([info[i-7].replace("_", "\_"), info[i+3], info[i+5]])
    print(v)
    doc.add_row(v, "MYSTUDY", "MYCASE", threshold)
    doc.close()


if __name__ == '__main__':
    test()

#-------------------------------------------------------------------------------
