#!/usr/bin/env python
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2013 EDF S.A.
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

import sys
import os.path
import subprocess

#===============================================================================
# Utility functions
#===============================================================================

#-------------------------------------------------------------------------------

def trunk_major_minor(prev_major, prev_minor):
    """
    For trunk, determine next version based on previous branches
    """

    # Current release cycle plans for 4 minor releases per major release

    major = int(prev_major)
    minor = int(prev_minor) + 1
    if minor > 3:
        major += 1
        minor = 0

    return str(major), str(minor)

#-------------------------------------------------------------------------------

def svn_version_from_branch(branch):
    """
    Determine branch version number from Subversion branch
    """

    major = None
    minor = None

    if branch[0:7] == 'Version':
        v = branch[7:]
        sep = v.find('_')
        if (sep > -1):
            major = v[0:sep]
            minor = v[sep+1:]

    return major, minor

#-------------------------------------------------------------------------------

def svn_version_from_tag(tag):
    """
    Determine branch version number from Subversion tag
    """

    major = None
    minor = None
    release = None
    extra = None

    if tag[0] == 'V':
        v = tag[1:].split('_', 3)
        if len(v) >= 1:
            major = v[0]
        if len(v) >= 2:
            minor = v[1]
        if len(v) >= 3:
            if v[2] > '-1':
                release = v[2]
        if len(v) >= 4:
            extra = v[3]

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def svn_version(srcdir):
    """
    Determine version information from Subversion.
    """

    major = ''
    minor = ''
    release = ''
    extra = ''

    # Use the XML output from Subversion, to avoid formatting/language
    # environment parsing issues

    try:
        from xml.dom import minidom
    except Exception:
        return

    # Get base info from 'svn info'

    url = None
    revision = None

    cmd = ['svn', 'info', '--xml', srcdir]
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output = p.communicate()
    if p.returncode != 0:
        return major, minor, release, extra
    else:
        try:
            svn_info = minidom.parseString(output[0]).documentElement
            url = svn_info.getElementsByTagName("url").item(0).firstChild.data
            entry_node = svn_info.getElementsByTagName("entry")[0]
            revision = entry_node.getAttribute("revision")
        except Exception:
            pass

    if not url:
        return major, minor, release, extra

    # If we are in trunk, use branches info to determine next version

    if url[-6:] == '/trunk':
        cmd = ['svn', 'ls', '--xml', url[:-5] + 'branches']

        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output = p.communicate()
        if p.returncode != 0:
            return major, minor, release, extra
        else:
            try:
                version_max = 'Version3_0'
                svn_info = minidom.parseString(output[0]).documentElement
                for name in svn_info.getElementsByTagName("name"):
                    branch = name.firstChild.data
                    if branch[0:7] == 'Version':
                        if branch > version_max:
                            version_max = branch
                v_maj, v_min = svn_version_from_branch(version_max)
                major, minor = trunk_major_minor(v_maj, v_min)
                release = ''
                extra += '-alpha-r' + revision
            except Exception:
                pass

    # If we are in a branch, use tags to determine next version

    elif url.find('/branches/') > -1:
        url_id_0 = url.find('/branches/')
        url_id_1 = url_id_0 + len('/branches/')
        url_root = url[:url_id_0]

        if url[url_id_1:url_id_1+7] == 'Version':
            major, minor = svn_version_from_branch(url[url_id_1:])

        cmd = ['svn', 'ls', '--xml', url_root + '/tags']

        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output = p.communicate()
        if p.returncode != 0:
            return major, minor, release, extra
        else:
            try:
                tag_ref = 'V' + str(major) + '_' + str(minor) + '_'
                tag_ref_l = len(tag_ref)
                name_max = tag_ref + '-1'
                release = ''
                svn_info = minidom.parseString(output[0]).documentElement
                for tag in svn_info.getElementsByTagName("name"):
                    name = tag.firstChild.data
                    if name[0:tag_ref_l] == tag_ref and name > name_max:
                        rev_node = tag.parentNode.getElementsByTagName("commit")[0]
                        tag_rev = rev_node.getAttribute("revision")
                        if (tag_rev == revision):
                            release = tag[tag_ref_l:]
                            break
                        elif tag_rev < revision:
                            name_max = name
                major, minor, release, extra_0 = svn_version_from_tag(name_max)
                if not release:
                    extra += '-beta'
                if tag_rev != revision:
                    extra += '-r' + revision
            except Exception:
                pass

    # If we are at a tag, simply use if

    elif url.find('/tags/') > -1:
        url_id_0 = url.find('/tags/')
        url_id_1 = url_id_0 + len('/tags/')
        tag = url[url_id_1:]
        major, minor, release, extra = svn_version_from_tag(url[url_id_1:])

    # If version is modified, indicate it

    cmd = ['svn', 'status', '--xml', srcdir]

    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output = p.communicate()
    if p.returncode != 0:
        return major, minor, release, extra
    else:
        if output[0].find('entry') > -1:
            extra += '-m'

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def git_svn_version(srcdir, url, revision):
    """
    Determine partial version information from Subversion.
    """

    major = ''
    minor = ''
    release = ''
    extra = ''

    # Use the XML output from Subversion, to avoid formatting/language
    # environment parsing issues

    # If we are in trunk, use branches info to determine next version

    if url[-6:] == '/trunk':
        release = ''
        extra += '-alpha-svn' + revision

        cmd = ['git', '--no-pager', '--git-dir='+os.path.join(srcdir, '.git'),
               '--work-tree='+srcdir, 'branch', '-r']
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        output = p.communicate()
        if p.returncode != 0:
            return major, minor, release, extra
        else:
            version_max = 'Version3_0'
            for branch in output[0].split():
                if branch[0:7] == 'Version':
                    if branch > version_max:
                        version_max = branch
                v_maj, v_min = svn_version_from_branch(version_max)
                major, minor = trunk_major_minor(v_maj, v_min)
                release = ''

    # If we are in a branch, use tags to determine next version

    elif url.find('/branches/') > -1:
        url_id_0 = url.find('/branches/')
        url_id_1 = url_id_0 + len('/branches/')
        url_root = url[:url_id_0]

        if url[url_id_1:url_id_1+7] == 'Version':
            major, minor = svn_version_from_branch(url[url_id_1:])
            extra += '-svn' + revision

    # If we are at a tag, simply use if

    elif url.find('/tags/') > -1:
        url_id_0 = url.find('/tags/')
        url_id_1 = url_id_0 + len('/tags/')
        tag = url[url_id_1:]
        major, minor, release, extra = svn_version_from_tag(url[url_id_1:])

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def git_version(srcdir):
    """
    Determine version information from Git.
    """

    major = ''
    minor = ''
    release = ''
    extra = ''

    # Use the XML output from Subversion, to avoid formatting/language
    # environment parsing issues

    try:
        from xml.dom import minidom
    except Exception:
        return

    # Get base info from 'svn info'

    url = None
    revision = None
    head_is_svn = False

    cmd = ['git', '--no-pager', '--git-dir='+os.path.join(srcdir, '.git'),
           '--work-tree='+srcdir, 'log', '-n', '200']
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output = p.communicate()
    if p.returncode != 0:
        return major, minor, release, extra
    else:
        svn_id = output[0].find('git-svn-id:')
        if svn_id > 1:
            svn_rev = output[0][svn_id:].split(' ', 2)[1]
            svn_url, svn_revision = svn_rev.split('@')
            major, minor, release, extra \
                = git_svn_version(srcdir, svn_url, svn_revision)
            # Try to count commits before the one with svn
            svn_log = output[0][:svn_id].split()
            c1 = svn_log.count('commit')
            c2 = svn_log.count('Author')
            if c1 == 1 and c2 == 2:
                head_is_svn = True

    if not head_is_svn:
        commit_id = output[0].find('commit')
        commit_rev = output[0][commit_id:].split(' ', 2)[1].split('\n')[0]
        extra += '-git' + commit_rev

    cmd = ['git', '--no-pager', '--git-dir='+os.path.join(srcdir, '.git'),
           '--work-tree='+srcdir, 'status', '-s', '-uno']
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    output = p.communicate()
    if p.returncode != 0:
        return major, minor, release, extra
    else:
        if len(output[0].split('\n')) > 1:
            extra += '-m'

    return major, minor, release, extra

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    # Simple parse for command line

    cmd = ''
    if len(sys.argv) == 1:
        cmd = '--help'
    elif sys.argv[1] in ['-h', '--help']:
        cmd = '--help'
    else:
        cmd = '--version'
        srcdir = sys.argv[1]

    # Call tools

    if cmd == '--version':
        if os.path.isdir(os.path.join(srcdir, ".svn")):
            major, minor, release, extra = svn_version(srcdir)
        elif os.path.isdir(os.path.join(srcdir, ".git")):
            major, minor, release, extra = git_version(srcdir)

        version = major + '.' + minor
        if release:
            version += '.' + release
        if extra:
            version += extra
        print(version)

    else:
        usage = \
            """Usage: %(prog)s <srcdir>

Options:
  -h, --help  show this help message and exit
"""
        sys.stderr.write(usage)
        sys.exit(1)

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
