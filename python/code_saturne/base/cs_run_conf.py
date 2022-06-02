#!/usr/bin/env python3
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

import sys, os

from collections import OrderedDict
import configparser

from code_saturne.base import cs_batch

#-------------------------------------------------------------------------------
# Get info on locally installed version
#-------------------------------------------------------------------------------

def get_install_config_info(pkg):
    """
    Return info based on install configuration:
    dictionnary with 'batch', 'submit_command', 'resource', 'compute_builds'
    """

    c = OrderedDict({'batch': None,
                     'submit_command': None,
                     'resource_name': None,
                     'compute_builds': None})

    # Use alternate compute (back-end) package if defined

    config = configparser.ConfigParser()
    config.read(pkg.get_configfiles())

    for o in ('batch', 'submit_command', 'resource_name'):
        if config.has_option('install', o):
            c[o] = config.get('install', o)

    if config.has_option('install', 'compute_versions'):
        c['compute_builds'] = config.get('install',
                                         'compute_versions').split(':')

    if c['resource_name']:
        c['resource_name'] = c['resource_name'].lower()

    return c

#-------------------------------------------------------------------------------
# Get local resource name based on installation info
#-------------------------------------------------------------------------------

def get_resource_name(i_c):
    """
    Return name of current resource, based on dictionnary returned
    from get_install_config_info.
    """

    # Resources: try to find a matching section, using
    # resource_name, batch, and job_defaults in decreasing priority.

    resource_name = i_c['resource_name']
    if not resource_name:
        resource_name = i_c['batch']
        if resource_name:
            resource_name = os.path.basename(resource_name).lower()
    if not resource_name:
        resource_name = 'job_defaults'

    return resource_name

#===============================================================================
# Class used to manage run configuration files
#===============================================================================

class run_conf(object):

    def __init__(self,
                 path,
                 package=None,
                 create_if_missing=True,
                 rebuild=False):
        """
        Initialize run configuration info object.
        """

        self.path = path
        self.package = package
        self.sections = OrderedDict()
        self.comments = {}

        if not path:
            return

        need_rebuild = rebuild
        force_rebuild = False

        lines = []

        try:
            f = open(self.path, mode = 'r')
            lines = f.readlines()
            f.close()

        except IOError:
            if create_if_missing or rebuild:
                need_rebuild = True
                force_rebuild = True
                lines = []
            else:
                print("Error: can not open or read %s\n" % self.path)
                sys.exit(1)

        for i in range(len(lines)):
            lines[i] = lines[i].rstrip()

        self.__parse__(lines)

        if need_rebuild:
            self.rebuild_resource(force_rebuild)

    #---------------------------------------------------------------------------

    def rebuild_resource(self, force_rebuild=False):
        """
        Rebuild dictionnary using reference templates
        """

        i_c = get_install_config_info(self.package)
        resource_name = get_resource_name(i_c)

        if not resource_name in self.sections:
            self.sections[resource_name] = OrderedDict()

        # Determine or build default job name for batch

        try:
            topdir, scriptdir = os.path.split(os.path.split(self.path)[0])
            if scriptdir == 'DATA':
                studydir, casedir = os.path.split(topdir)
                studydir = os.path.split(studydir)[1]
            else:
                casedir = ''
                studydir = scriptdir
        except Exception:
            casedir = ''
            studydir = ''

        job_name = studydir.lower()
        if casedir:
            job_name += '_' + casedir.lower()

        job_lines = cs_batch.generate_header(job_name=job_name,
                                             package=self.package)

        if job_lines:
            if force_rebuild or not 'job_header' in self.sections[resource_name]:
                self.sections[resource_name]['job_header'] \
                    = os.linesep.join(job_lines)
        else:
            for k in ('n_procs', 'n_threads'):
                if not k in self.sections[resource_name]:
                    self.sections[resource_name][k] = ''

        return

    #---------------------------------------------------------------------------

    def __strip_extra_lines__(self, value):
        """
        Remove start and end empty lines
        """

        br = os.linesep
        lbr = len(br)

        while value[:lbr] == br:
            value = value[lbr:]
        while value[-lbr:] == br:
            value = value[:-lbr]

        return value

    #---------------------------------------------------------------------------

    def __parse__(self, lines):
        """
        Parse the file.
        """

        br = os.linesep
        lbr = len(br)

        sections = OrderedDict()
        comments = {}

        section_name = 'job_defaults'
        section_dict = OrderedDict()
        comment_dict = {}
        sub_section = False
        key_indent = -1  # not inside key if < 0
        key = ''
        value = ''
        comment = ''

        for line in lines:

            if sub_section:
                l = line.lstrip()
                if l[:1] == '[' and l.find(']') > 0:
                    sub_section = False
                    section_dict[key] = self.__strip_extra_lines__(value)
                    key = ''
                    value = ''
                else:
                    value += br + line
                    continue

            # determine indent and effective content
            l = line.lstrip()
            indent = len(line) - len(l)

            # handle multiline
            if key_indent > -1:
                if indent > key_indent or len(line) == 0:
                    if len(value) > 0:
                        value += br + l
                    else:
                        value += l
                    continue
                else:
                    key_indent = -1

            # ignore empty lines
            if not line:
                if comment != '':
                    comment += br
                continue

            if key_indent < 0 and key:
                section_dict[key] = self.__strip_extra_lines__(value)
                key = ''
                value = ''

            if l[:1] in ('#', ';'):  # comment line
                comment += line + br
                continue

            if l[:1] == '[':
                i = l.find(']')
                if i > 0:
                    # New section; save previous one
                    sections[section_name] = section_dict
                    comments[section_name] = comment_dict
                    section_name = l[1:i].strip().lower()
                    # Possible subsection
                    j = section_name.find(':')
                    if j > -1:
                        key = section_name[j+1:]
                        section_name = section_name[:j]
                        sub_section = True
                        value = ''
                    if section_name in sections:
                        section_dict = sections[section_name]
                        comment_dict = comments[section_name]
                    else:
                        section_dict = OrderedDict()
                        comment_dict = {}
                    if comment != '':
                        if '' in comment_dict:
                            comment_dict[''] += br + comment.rstrip()
                        else:
                            comment_dict[''] = comment.rstrip()
                        comment = ''
            else:
                for i, c in enumerate(l):
                    if c in (':', '='):
                        key = l[:i].rstrip().lower()
                        value = l[i+1:].lstrip()
                        key_indent = indent
                        if comment != '':
                            comment_dict[key] = comment.rstrip()
                            comment = ''
                        break

        if sub_section or key:
            section_dict[key] = self.__strip_extra_lines__(value)

        sections[section_name] = section_dict
        comments[section_name] = comment_dict

        self.sections = sections
        self.comments = comments

    #---------------------------------------------------------------------------

    def __rebuild_lines__(self):
        """
        Update lines based on the current dictionnary values.
        """

        br = os.linesep
        lbr = len(br)

        # The following keys will use multiline sections by default
        multiline_section_keys = ('job_header',
                                  'run_prologue', 'run_epilogue',
                                  'compute_prologue', 'compute_epilogue')

        # Default sections we prefer to place first
        default_sections = ['setup', 'job_defaults', 'run']

        other_sections = []
        for sn in self.sections:
            if sn not in default_sections:
                other_sections.append(sn)
        for sn in self.comments:
            if sn not in default_sections:
                other_sections.append(sn)
        other_sections = list(set(other_sections))
        other_sections.sort()

        section_queue = default_sections + other_sections

        lines = []

        prev_is_multiline = False

        # now loop on sections

        for sn in section_queue:

            sn_comment = None
            if sn in self.comments:
                if '' in self.comments[sn]:
                    sn_comment = self.comments[sn][''].split(br)
                    if prev_is_multiline:
                        lines.append('[]')
                        lines.append('')
                        prev_is_multiline = False
                    lines += sn_comment
                    lines.append('')

            keys = []
            if sn in self.sections:
                for k in self.sections[sn]:
                    if self.sections[sn][k] != None:
                        keys.append(k)
            if sn in self.comments:
                for k in self.comments[sn]:
                    if self.comments[sn][k] != None:
                        keys.append(k)

            keys = list(set(keys))
            keys.sort()

            if sn == 'run':
                for k in ('stage', 'initialize', 'compute', 'finalize'):
                    if k in keys:
                        keys.remove(k)
                        keys.append(k)

            # Loop on keys with at least a value or comment

            simple_keys = []
            multiline_keys = []

            for k in keys:
                if k == '':
                    continue
                if k in multiline_section_keys:
                    multiline_keys.append(k)
                else:
                    simple_keys.append(k)

            has_simple = False

            if len(simple_keys) > 0 or sn_comment != None:
                lines.append('[' + sn.lower() + ']')
                lines.append('')
                prev_is_multiline = False
                has_simple = True

            k_value_l = 0

            # Standard keys first

            for k in simple_keys:

                k_comment = None
                k_value = None
                k_value_l = 0

                if sn in self.comments:
                    if k in self.comments[sn]:
                        k_comment = self.comments[sn][k].split(br)
                if sn in self.sections:
                    if k in self.sections[sn]:
                        k_value = str(self.sections[sn][k]).split(br)
                        k_value_l = len(k_value)

                if k_comment:
                    lines += k_comment
                if k_value_l > 0:
                    if k_value_l == 1:
                        v = k_value[0]
                        if v in ('True', 'False'):
                            v = v.lower()
                        lines.append(k + ': ' + v)
                    else:
                        if k_comment:
                            lines.append('')
                        lines.append(k + ': ')
                        for l in k_value:
                            lines.append('    ' + l)
                        lines.append('')
                elif k_comment:
                    lines.append(k + ': ')

            if has_simple:
                if k_value_l < 2:
                    lines.append('')

            # Multiline keys next

            for k in multiline_keys:

                k_comment = None
                k_value = None

                if sn in self.comments:
                    if k in self.comments[sn]:
                        k_comment = self.comments[sn][k].split(br)
                if sn in self.sections:
                    if k in self.sections[sn]:
                        k_value = str(self.sections[sn][k]).split(br)

                if k_comment:
                    if prev_is_multiline:
                        lines.append('[]')
                        lines.append('')
                    lines += k_comment
                    lines.append('')

                lines.append('[{0}:{1}]'.format(sn, k))
                lines.append('')

                if k_value:
                    lines += k_value
                    lines.append('')

                prev_is_multiline = True

        # End of loop on sections

        return lines

    #---------------------------------------------------------------------------

    def get(self, section, key):
        """
        Get a section:value combination (as a string, or None if not present)
        """

        value = None
        if section in self.sections:
            if key in self.sections[section]:
                value = str(self.sections[section][key])

        return value

    #---------------------------------------------------------------------------

    def get_bool(self, section, key):
        """
        Get a section:value combination as a boolean type, or None if not present
        """

        value = None
        if section in self.sections:
            if key in self.sections[section]:
                value = str(self.sections[section][key])

        if value != None:
            if value in ('0', 'f', 'false', 'n', 'no'):
                value = False
            elif value in ('1', 't', 'true', 'y', 'yes'):
                value = True
            elif value == '':
                value = None
            else:
                if self.path:
                    err_str = 'In file: ' + self.path + ':' + os.linesep
                else:
                    err_str = 'In run_conf object (' + str(self) + '):' + os.linesep
                    err_str += '  ['+section+']' + os.linesep
                    err_str += '  '+key+': ' + str(value) + os.linesep
                    err_str += 'is not a boolean value.'
                    raise TypeError(err_str)

        return value

    #---------------------------------------------------------------------------

    def get_int(self, section, key):
        """
        Get a section:value combination as an integer, or None if not present
        """

        value = None
        if section in self.sections:
            if key in self.sections[section]:
                value = str(self.sections[section][key])

        if value:
            try:
                value = int(value)
            except Exception:
                if self.path:
                    err_str = 'In file: ' + self.path + ':' + os.linesep
                else:
                    err_str = 'In run_conf object (' + str(self) + '):' + os.linesep
                err_str += '  ['+section+']' + os.linesep
                err_str += '  '+key+': ' + str(value) + os.linesep
                err_str += 'is not an integer value.'
                raise TypeError(err_str)

        return value

    #---------------------------------------------------------------------------

    def set(self, section, key, value):
        """
        Set a section:key value combination
        """

        if not section in self.sections:
            self.sections[section] = OrderedDict()
        if value != None:
            self.sections[section][key] = str(value)
        else:
            self.sections[section][key] = None

    #---------------------------------------------------------------------------

    def save(self, path=None, new=False):
        """
        Save run configuration, optionally to a different path.
        """

        br = os.linesep

        if path:
            self.path = path
            if not new:
                self.lines = []
                try:
                    f = open(self.path, mode = 'r')
                    self.lines = f.readlines()
                    f.close()
                    for i in range(len(self.lines)):
                        self.lines[i] = self.lines[i].rstrip()
                except IOError:
                    pass

        lines = self.__rebuild_lines__()

        f = open(self.path, mode='w')
        for line in lines:
            f.write(line + br)
        f.close()

    #---------------------------------------------------------------------------

    def get_coupling_parameters(self):
        """
        Return the list of coupled domains defined inside a run.cfg.
        Empty list if single-case run.
        """

        domains = []

        domain_names = self.get('setup', 'coupled_domains')

        if domain_names:
            for dom in domain_names.split(":"):
                d = self.sections[dom.lower()]
                domains.append(OrderedDict())

                for key in d.keys():
                    if d[key] and d[key] != 'None':
                        if key in ('n_procs_max', 'n_procs_min', 'n_procs_weight'):
                            domains[-1][key] = int(d[key])
                        else:
                            domains[-1][key] = str(d[key])
                    else:
                        domains[-1][key] = None

        return domains


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
