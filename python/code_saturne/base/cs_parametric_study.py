#!/usr/bin/env python3

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

"""
Update XML file with common parametric arguments.
"""

import os, sys
from argparse import ArgumentParser

from code_saturne.base import cs_parametric_setup

#-------------------------------------------------------------------------------
# Process the command line
#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """

    setup_parser = cs_parametric_setup.arg_parser(argv)
    prog = os.path.basename(sys.argv[0]) + " " + sys.argv[1]

    parser = ArgumentParser(parents=[setup_parser],
                            prog=prog,
                            description="Run a parametric study",
                            conflict_handler='resolve')

    parser.add_argument("--update", dest="updateXml", default = False,
                        action= 'store_true',
                        help="Update the xml file.")

    parser.add_argument("--create-from-ref", dest="createFromRef", type=str,
                        help="Reference case to create a copy of")

    parser.add_argument("-c", "--case", dest="case", type=str,
                        help="Directory of the current case")

    parser.add_argument("-p", "--param", dest="param", type=str,
                        help="Name f the file of parameters")

    (options, args) = parser.parse_known_args(argv)

    if args != []:
        print(" -------------------------- ")
        print(" Ignored arguments : ")
        print(" ", args)
        print(" -------------------------- ")

    return options

#-------------------------------------------------------------------------------
# Load case from XML file
#-------------------------------------------------------------------------------

def load_case_setup_from_xml(xml_file, pkg=None):
    """
    Load case from XML file.
    @param xml_file: path to xml file
    """

    pkg = None
    if not pkg:
        from code_saturne.base.cs_package import package
        pkg = package()

    from code_saturne.model.XMLengine import Case

    try:
        case_setup = Case(package=pkg, file_name=xml_file)
        case_setup['xmlfile'] = xml_file
    except:
        print("Error while reading parameters files.")
        print("This file is not in accordance with XML specifications.")
        sys.exit(1)

    module_name = case_setup.module_name()

    if module_name == 'code_saturne':
        from code_saturne.model.XMLinitialize import XMLinit
    elif module_name == 'neptune_cfd':
        from code_saturne.model.XMLinitializeNeptune import XMLinitNeptune as XMLinit

    # internal functions of code_saturne to properly initialize the
    # xml structure
    case_setup.xmlCleanAllBlank(case_setup.xmlRootNode())
    XMLinit(case_setup).initialize()

    return case_setup

# ------------------------------------------------------------------------------
# Class which allows the execution of code_saturne script commands
# ------------------------------------------------------------------------------

class cs_exec(object):

    def __init__(self, pkg=None):

        # package
        self.pkg = None
        if pkg:
            self.pkg = pkg
        else:
            from code_saturne.base.cs_package import package
            self.pkg = package()

    #---------------------------------------------------------------------------

    def run_cs_command(self, args):
        """
        Execute code_saturne with a given command line.
        @param args: list of arguments composing the command line
        """

        from code_saturne.base.cs_script import master_script
        script = master_script(argv=args, package=self.pkg)

        script.execute()
        pass

    #---------------------------------------------------------------------------

    def run_cs_case(self, case_dir, xmlfile="setup.xml", nprocs=1):
        """
        Run a code_saturne case.
        @param case_dir: path to the case folder
        @param xmlfile : name of the xml file (default is setup.xml)
        @param nprocs  : number of cores to uses (default is 1)
        """

        # Check where we are before launching
        origin = os.getcwd()

        # Go to case folder
        os.chdir(os.path.join(case_dir, "DATA"))

        # Set the arguments list
        args = ["run", "--param", xmlfile, "--nprocs", str(nprocs)]

        self.run_cs_command(args)

        os.chdir(origin)

    #---------------------------------------------------------------------------

    def generate_new_case(self, new_case=None, ref_case=None):
        """
        Generate a new case based on a reference case.
        @param new_case: name of new case
        @param ref_case: path to reference case
        """

        args = ["create", "-c"]

        # Set case name
        if new_case:
            args.append(new_case)
        else:
            args.append("CASE1")

        # Set ref case if needed
        if ref_case:
            args.extend(["--copy-from", ref_case])

        self.run_cs_command(args)

#-------------------------------------------------------------------------------
# Function which will update an xml file based on the cs_modify_xml class
#-------------------------------------------------------------------------------

def update_xml_file(pkg, filepath, options):

    # Sanity check
    if not os.path.isfile(filepath):
        print("File %s does not exist!" % filepath)
        sys.exit(1)

    case_setup = load_case_from_xml(xml_file, pkg)
    cs_parametric_setup.update_case_model(case_setup, options, pkg)

    case_setup.xmlSaveDocument()

#-------------------------------------------------------------------------------
# Run main function which modifies the case parameters
#-------------------------------------------------------------------------------

def main(args, pkg):

    # Read input arguments
    # --------------------

    if args == []:
        args.append("--help")
    opts = process_cmd_line(args)

    # Check if one needs to create a case from a reference one
    # --------------------------------------------------------

    if opts.createFromRef:
        cmd_mgr = cs_exec(pkg)
        cmd_mgr.generate_new_case(opts.case, opts.createFromRef)

        # Set default name if needed
        if opts.case is None:
            opts.case = "CASE1"

    # Sanity check of case or xml in input parameters
    # -----------------------------------------------

    if opts.case is None and opts.param is None:
        print("No case name nor parameters file provided")
        sys.exit(1)

    fp = ""
    if opts.case:
        fp = os.path.join(opts.case, "DATA")
        if opts.param:
            fp = os.path.join(fp, opts.param)
        else:
            fp = os.path.join(fp, "setup.xml")

    elif opts.param:
        fp = opts.param

    if fp == "":
        print("No case name nor parameters file provided")
        sys.exit(1)

    # Update xml if needed
    # --------------------

    if opts.updateXml:
        update_xml_file(pkg, fp, opts)

#===============================================================================
# Main function for standalone use
#===============================================================================

if __name__=="__main__":

    args = None

    if len(sys.argv) > 1:
        args = sys.argv[1:]
    else:
        args = ["--help"]

    main(args=args, pkg=None)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
