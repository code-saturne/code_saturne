# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
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
    def __init__(self, dest, filename, pdflatex):
        self.__dest = dest
        self.__filename = os.path.join(self.__dest, filename)
        self.__doc = []
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
            elif columns[i] is None:
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
        Build the pdf file, and clean the temporary files.
        """
        cmd = "pdflatex " + self.__filename + ".tex"
        log_name = "make_pdf.log"
        log_file = open(log_name, mode='w')
        error, time = run_studymanager_command(cmd, log_file)
        log_file.close()

        if not error:
            os.remove(log_name)
            for suffixe in ["tex", "log", "aux"]:
                f = self.__filename + "." + suffixe
                if os.path.isfile(f):
                    os.remove(self.__filename + "." + suffixe)
        else:
            print(" /!\ ERROR during pdf generation. See %s\n",
                  log_name)

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

class Report(TexWriter):
    """
    Figures report.
    """
    def __init__(self, dest, label, pdflatex):
        """
        """
        TexWriter.__init__(self, dest, label, pdflatex)

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
