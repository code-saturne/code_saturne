#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

import os, sys
import string
import getpass
import platform

from optparse import OptionParser
from datetime import datetime, date

import smtplib

try: # email version 3.0 (Python2 up to 2.6)
    from email.Utils import COMMASPACE, formatdate
    from email import Encoders
    from email.MIMEMultipart import MIMEMultipart
    from email.MIMEBase import MIMEBase
    from email.MIMEText import MIMEText
except Exception: # email version 4.0 (Python2 from Python 2.5)
    from email.utils import COMMASPACE, formatdate
    from email import encoders
    from email.mime.multipart import MIMEMultipart
    from email.mime.base import MIMEBase
    from email.mime.text import MIMEText

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

# Do not import studymanager yet, as it pulls Python packages such as
# matplotlib which may not be in the standard path, and may need
# sourcing of a specific environment (which itself is delayed in case
# the main and coputation packages are not the same).

#-------------------------------------------------------------------------------
# Processes the passed command line arguments
#-------------------------------------------------------------------------------

def process_cmd_line(argv, pkg):
    """
    Processes the passed command line arguments.
    """
    if sys.argv[0][-3:] == '.py':
        usage = "usage: %prog <options>"
    else:
        usage = "usage: %prog studymanager <options>"

    parser = OptionParser(usage=usage)

    parser.add_option("-f", "--file", dest="filename", type="string",
                      metavar="FILE", help="xml FILE of parameters")

    parser.add_option("-q", "--quiet",
                      action="store_true", dest="quiet", default=False,
                      help="don't print status messages to stdout")

    parser.add_option("-u", "--update",
                      action="store_true", dest="update", default=False,
                      help="update scripts in the repository")

    parser.add_option("-x", "--update-xml",
                      action="store_true", dest="update_xml", default=False,
                      help="update only xml files in the repository")

    parser.add_option("-t", "--test-compilation",
                      action="store_true", dest="test_compilation",
                      default=False, help="compile all cases")

    parser.add_option("-r", "--run",
                      action="store_true", dest="runcase", default=False,
                      help="run all cases")

    parser.add_option("--n-procs",  dest="n_procs", default=None, type="int",
                      help="Optional number of processors requested for the computations")

    parser.add_option("-n", "--n-iterations", dest="n_iterations",
                      type="int", help="maximum number of iterations for cases of the study")

    parser.add_option("-c", "--compare",
                      action="store_true", dest="compare", default=False,
                      help="compare results between repository and destination")

    parser.add_option("-d", "--ref-dir", dest="reference", type="string",
                      help="absolute reference directory to compare dest with")

    parser.add_option("-p", "--post",
                      action="store_true", dest="post", default=False,
                      help="postprocess results of computations")

    parser.add_option("-m", "--mail", dest="addresses", default="",
                      type="string", metavar="ADDRESS1 ADDRESS2 ...",
                      help="addresses for sending the reports")

    parser.add_option("-l", "--log", dest="log_file", default="studymanager.log",
                      type="string",
                      help="name of studymanager log file (default value is 'studymanager.log'")

    parser.add_option("-z", "--disable-tex",
                      action="store_true", dest="disable_tex", default=False,
                      help="disable text rendering with LaTex in Matplotlib (use Mathtext)")

    parser.add_option("--rm",
                      action="store_true", dest="remove_existing", default=False,
                      help="remove existing run directories")

    parser.add_option("--fow",
                      action="store_true", dest="force_overwrite", default=False,
                      help="overwrite files in MESH and POST directories")

    parser.add_option("-s", "--skip-pdflatex", default=False,
                      action="store_true", dest="disable_pdflatex",
                      help="disable tex reports compilation with pdflatex")

    parser.add_option("--fmt", type="string",
                      dest="default_fmt", default="pdf",
                      help="set the global format for exporting matplotlib figure (default is pdf)")

    parser.add_option("--repo", type="string",
                      dest="repo_path", default="",
                      help="force the path to the repository directory")

    parser.add_option("--dest", type="string",
                      dest="dest_path", default="",
                      help="force the path to the destination directory")

    parser.add_option("-g", "--debug", action="store_true",
                      dest="debug", default=False,
                      help="activate debugging mode")

    parser.add_option("--with-tags", type="string",
                      dest="with_tags", default="",
                      help="only process runs with all specified tags (separated by commas)")

    parser.add_option("--without-tags", type="string",
                      dest="without_tags", default="",
                      help="exclude any run with one of specified tags (separated by commas)")

    (options, args) = parser.parse_args(argv)

    return  options

#-------------------------------------------------------------------------------
# Send the report.
#-------------------------------------------------------------------------------

def send_report(code_name, report, labels, to, files):
    """
    Build the mail to be send.
    """
    # Note that the TO variable must be a list, and that you have
    # to add the From, To, and Subject headers to the message yourself.
    # The TO argument to the sendmail method is only used to route the
    # message to the right recipients, and doesn't have to match the headers
    # in the message itself. This can be used to implement Bcc headers; to send
    # a blind copy, add the address to the TO argument, but do not include
    # it in the To header.

    if code_name == "Code_Saturne":
        FROM    = "saturne-support@edf.fr"
    else:
        FROM    = "neptune-cfd-support@edf.fr"

    SERVER  = "localhost"
    TO      = to
    SUBJECT = "%s. Study Manager %s: %s " % (code_name, date.today(), labels)
    TEXT    = report
    RETOUR  = "An error occurs during the sending of the log mail"
    FILES   = files
    MESSAGE = """From: %s\nTo: %s\nSubject: %s\n\n%s""" % (FROM, ", ".join(TO), SUBJECT, TEXT)
    retour = "\n"

    try:
        send_mail(FROM, TO, SUBJECT, MESSAGE, FILES, SERVER)
    except:
        send_mail(FROM, TO, "[ERROR] %s" % SUBJECT, retour, [], SERVER)

#-------------------------------------------------------------------------------
# Send an email
#-------------------------------------------------------------------------------

def send_mail(send_from, send_to, subject, text, files=[], server="localhost"):
    """
    Send the report by mail.
    """
    assert type(send_to) == list
    assert type(files)   == list

    msg = MIMEMultipart()
    msg['From']    = send_from
    msg['To']      = COMMASPACE.join(send_to)
    msg['Date']    = formatdate(localtime=True)
    msg['Subject'] = subject
    msg.attach( MIMEText(text) )

    for f in files:
        part = MIMEBase('application', "octet-stream")
        part.set_payload( open(f,"rb").read())
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition',
                        'attachment; filename="%s"' % os.path.basename(f))
        msg.attach(part)

    smtp = smtplib.SMTP(server)
    smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.close()

#-------------------------------------------------------------------------------
# Helper to find the release of Linux
#-------------------------------------------------------------------------------

def release():
    p = "/etc/issue"
    if os.path.isfile(p):
        f = open(p, 'r')
        issue = f.read()[:-1]
        f.close()
    else:
        issue = ""
    return issue

#-------------------------------------------------------------------------------
# Start point of studymanager script
#-------------------------------------------------------------------------------

def run_studymanager(pkg, options):
    """
    Main function
      1. parse the command line,
      2. read the file of parameters
      3. create all studies,
      4. compile sources
      5. run all cases
      6. compare results
      7. plot result
      8. reporting by mail
    """

    # Source environment if required before importing studymanager modules, as
    # it pulls Python packages such as matplotlib which may not be in
    # the standard path.
    from cs_exec_environment import set_modules, source_rcfile, enquote_arg
    set_modules(pkg)
    source_rcfile(pkg)

    from studymanager.cs_studymanager_study import Studies

    # Scripts

    exe = os.path.join(pkg.get_dir('bindir'), pkg.name)
    if sys.platform.startswith('win'):
        exe = exe + ".com"
    exe = enquote_arg(exe)

    dif = pkg.get_io_dump()

    for p in exe, dif:
        if not os.path.isfile(p):
            print("Error: executable %s not found." % p)
            if not sys.platform.startswith('win'):
                return 1

    dif += " -d"

    # Read the file of parameters

    studies = Studies(pkg, options, exe, dif)
    if options.debug:
        print(" run_studymanager() >> Studies are initialized")
    if options.update_xml == False:
        os.chdir(studies.getDestination())
    else:
        os.chdir(studies.getRepository())

    # Print header

    studies.reporting(" -------------")
    studies.reporting(" Study Manager")
    studies.reporting(" -------------\n")
    studies.reporting(" Code name:         " + pkg.name)
    studies.reporting(" Kernel version:    " + pkg.version)
    studies.reporting(" Install directory: " + pkg.get_dir('exec_prefix'))
    studies.reporting(" File dump:         " + dif)
    studies.reporting(" Repository:        " + studies.getRepository())
    studies.reporting(" Destination:       " + studies.getDestination())
    studies.reporting("\n Informations:")
    studies.reporting(" -------------\n")
    studies.reporting(" Date:               " + datetime.now().strftime("%A %B %d %H:%M:%S %Y"))
    studies.reporting(" Platform:           " + platform.platform())
    studies.reporting(" Computer:           " + platform.uname()[1] + "  " + release())
    studies.reporting(" Process Id:         " + str(os.getpid()))
    studies.reporting(" User name:          " + getpass.getuser())
    studies.reporting(" Working directory:  " + os.getcwd())
    studies.reporting("\n")

    # Update repository if needed

    if options.update:
        studies.updateRepository()

    # Update only xml data if needed

    if options.update_xml:
        studies.updateRepository(True)

    # Test sources compilation for all cases

    if options.test_compilation:
        studies.test_compilation()

    # Check if xml for result directories in the repository are OK

    if options.compare:
        studies.check_compare(destination=False)

    if options.post:
        studies.check_script(destination=False)
        studies.check_plots_and_input(destination=False)

    # Create all studies and all cases

    if options.compare or options.post or options.runcase:
        studies.create_studies()

    # Preprocessing and run all cases
    if options.debug:
        print(" run_studymanager() >> Starts running...")

    if options.runcase:
        studies.run()

    if options.debug:
        print(" run_studymanager() >> Exits runs")

    # Compare checkpoint files

    if options.compare:
        studies.check_compare()
        studies.compare()

    # Postprocess results and probes

    if options.post:
        checked_scripts = studies.check_script()
        if checked_scripts:
            studies.scripts()
        studies.check_plots_and_input()
        studies.postpro()
        studies.plot()

    studies.reporting("\n --------------------")
    studies.reporting(" End of Study Manager")
    studies.reporting(" --------------------")

    # Reporting - attached files are either pdf or
    # raw tex files if pdflatex is disabled
    attached_file = studies.build_reports("report_global",
                                          "report_detailed")

    if len(options.addresses.split()) > 0:
        send_report(pkg.code_name, studies.logs(),
                    studies.getlabel(),
                    options.addresses.split(),
                    attached_file)

    return 0

def studymanager_usage():

    usage = \
        """Usage: %(prog) studymanager <options>

Options:
  -h, --help            show this help message and exit
  -f FILE, --file=FILE  xml FILE of parameters
  -q, --quiet           don't print status messages to stdout
  -u, --update          update scripts in the repository
  -x, --update-xml      update only xml files in the repository
  -t, --test-compile    compile all cases
  -r, --run             run all cases
  -n N_ITERATIONS, --n-iterations=N_ITERATIONS
                        maximum number of iterations for cases of the study
  -c, --compare         compare results between repository and destination
  -d REFERENCE, --ref-dir=REFERENCE
                        absolute reference directory to compare dest with
  -p, --post            postprocess results of computations
  -m ADDRESS1 ADDRESS2 ..., --mail=ADDRESS1 ADDRESS2 ...
                        addresses for sending the reports
  -l LOG_FILE, --log=LOG_FILE
                        name of studymanager log file (default value is
                        'studymanager.log'
  -z, --disable-tex     disable text rendering with LaTex in Matplotlib (use
                        Mathtext)
  --rm                  remove existing run directories
  -s, --skip-reports    disable the generation of reports
  --fmt=DEFAULT_FMT     Set the default format for exporting matplotlib figuer
  --repo=REPO_PATH      Force the path to the repository directory
  --dest=DEST_PATH      Force the path to the destination directory
  -g, --debug           Activate debugging mode
"""

    print(usage % {'prog':sys.argv[0]})

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    if argv == []:
        studymanager_usage()
        sys.exit(1)
    else:
        # Process command line
        options = process_cmd_line(argv, pkg)

        retcode = run_studymanager(pkg, options)
        sys.exit(retcode)

#-------------------------------------------------------------------------------
