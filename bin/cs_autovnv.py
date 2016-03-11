#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2016 EDF S.A.
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

# Do not import Autovnv yet, as it pulls Python packages such as
# matplotlib or vtk which may not be in the standard path, and may need
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
        usage = "usage: %prog [options]"
    else:
        usage = "usage: %prog autovnv [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-f", "--file", dest="filename", type="string",
                      metavar="FILE", help="xml FILE of parameters")

    parser.add_option("-q", "--quiet",
                      action="store_true", dest="verbose", default=False,
                      help="don't print status messages to stdout")

    parser.add_option("-u", "--update",
                      action="store_true", dest="update", default=False,
                      help="update scripts in the repository")

    parser.add_option("-x", "--update-xml",
                      action="store_true", dest="updatexml", default=False,
                      help="update only xml files in the repository")

    parser.add_option("-t", "--test-compile",
                      action="store_true", dest="compile", default=False,
                      help="compile all cases")

    parser.add_option("-r", "--run",
                      action="store_true", dest="runcase", default=False,
                      help="run all cases")

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

    parser.add_option("-l", "--log", dest="log_file", default="auto_vnv.log",
                      type="string",
                      help="name of autovnv log file (default value is 'auto_vnv.log'")

    parser.add_option("-z", "--disable-tex",
                      action="store_true", dest="disable_tex", default=False,
                      help="disable text rendering with LaTex in Matplotlib (use Mathtext)")

    (options, args) = parser.parse_args(argv)

    return  options.filename, options.verbose, \
        options.update, options.updatexml, options.compile, \
        options.runcase, options.n_iterations, options.compare, \
        options.reference, options.post, options.log_file, options.addresses, \
        options.disable_tex

#-------------------------------------------------------------------------------
# Send the report.
#-------------------------------------------------------------------------------

def sendmail(code_name, report, labels, to, files):
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
    SUBJECT = "%s. Auto V&V %s: %s " % (code_name, date.today(), labels)
    TEXT    = report
    RETOUR  = "An error occurs during the sending of the log mail"
    FILES   = files
    MESSAGE = """From: %s\nTo: %s\nSubject: %s\n\n%s""" % (FROM, ", ".join(TO), SUBJECT, TEXT)
    retour = "\n"

    try:
        send_mail(FROM, TO, SUBJECT, MESSAGE, FILES, SERVER)
    except:
        send_mail(FROM, TO, "[ERROR] %s" % SUBJECT, retour, [], SERVER)


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
# Start point of Auto V & V script
#-------------------------------------------------------------------------------

def runAutoverif(pkg, opt_f, opt_v, opt_u, opt_x, opt_t, opt_r, opt_n, opt_c, opt_d, opt_p, opt_l, opt_to, opt_z):
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

    # Source environment if required before importing Autovnv modules, as it
    # pulls Python packages such as matplotlib or vtk which may not be in the
    # standard path.
    from cs_exec_environment import set_modules, source_rcfile, enquote_arg
    set_modules(pkg)
    source_rcfile(pkg)

    from autovnv.Study import Studies

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
                sys.exit(1)

    dif += " -d"

    # Read the file of parameters

    studies = Studies(pkg, opt_f, opt_v, opt_x, opt_r, opt_n, opt_c, opt_d, opt_p, exe, dif, opt_l, opt_z)
    if opt_x == False:
        os.chdir(studies.getDestination())
    else:
        os.chdir(studies.getRepository())

    # Print header

    studies.reporting(" ----------")
    studies.reporting(" Auto V & V")
    studies.reporting(" ----------\n")
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

    if opt_u:
        studies.updateRepository()
        sys.exit(0)

    # Update only xml data if needed

    if opt_x:
        studies.updateRepository(True)
        sys.exit(0)

    # Check if xml for result directories in the repository are OK

    if opt_c:
        studies.check_compare(destination=False)

    if opt_p:
        studies.check_script(destination=False)
        studies.check_plot(destination=False)

    # Create all studies and all cases

    studies.createStudies()

    # Compile sources of all cases

    if opt_t:
        studies.compilation()

    # Preprocessing and run all cases

    if opt_r:
        studies.run()

    # Compare checkpoint files

    if opt_c:
        studies.check_compare()
        studies.compare()

    # Postprocess results and probes

    if opt_p:
        studies.check_script()
        studies.scripts()
        studies.check_plot()
        studies.postpro()
        studies.plot()

    studies.reporting("\n -----------------")
    studies.reporting(" End of Auto V & V")
    studies.reporting(" -----------------")

    # Reporting

    attached_file = studies.build_reports("report_global", "report_detailed")
    if opt_to:
        sendmail(pkg.code_name, studies.logs(), studies.getlabel(), opt_to, attached_file)

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

def main(argv, pkg):
    """
    Main function.
    """

    # Command line

    opt_f, opt_v, opt_u, opt_x, opt_t, opt_r, opt_n, opt_c, opt_d, opt_p, opt_l, addresses, opt_z = process_cmd_line(argv, pkg)
    opt_to  = addresses.split()

    retcode = runAutoverif(pkg, opt_f, opt_v, opt_u, opt_x, opt_t, opt_r, opt_n, opt_c, opt_d, opt_p, opt_l, opt_to, opt_z)
    sys.exit(retcode)

#-------------------------------------------------------------------------------
