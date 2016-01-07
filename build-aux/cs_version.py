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

import sys
import os.path

#-------------------------------------------------------------------------------

# Generate version strings, using local and optionally version control info.

# Notes:
#-------
# Currently, basic version information is determined from the NEWS file.
# This is not quite standard (though methods and automation level used by other
# tools vary widely), but we want to minimize files that need to be edited
# when releasing a version, but use a robust mechanism which will work
# for versions pulled by different version control tools (at least
# Subversion, used for the main repository, git svn clones, which may or
# may not have up-to-date fetched remotes information, and clones of those
# clones...). As we need to maintain the NEWS file anyways to provide a clear
# user-readable summary of changes between versions, using a consistent
# presentation, using this file to determine version information is a
# practical solution, which requires enforcing that the NEWS file is
# maintained correctly, which is good practice.

#===============================================================================
# Utility functions
#===============================================================================

#-------------------------------------------------------------------------------

def version_from_news(srcdir):
    """
    Determine version information from NEWS file.
    """

    major = 'unknown'
    minor = ''
    release = ''
    extra = ''

    trunk = -1
    version = [-1]
    release = [-1, '', '', '']

    f = open(os.path.join(srcdir, "NEWS"))
    lines = f.readlines()
    f.close()

    for i in range(1, len(lines)):

        # Release names are underlined by '==='

        if lines[i][0:7] == '=======':

            j = i-1
            l = lines[j]

            if l[0:6] == 'Trunk ' and trunk < 0:
                trunk = j

            # If release appears, at which line, and is there a date,
            # or is it not yet released ?

            elif l[0:8] == 'Release ' or l[0:8] == 'Version ':

                v, n, info = l.split(' ', 2)

                if version[0] < 0:
                    version = [j] + n.split('.')

                if release[0] < 0:
                    s = ''.join(c for c in info if c.isdigit())
                    if len(s) > 0:
                            release = [j] + n.split('.')
                            break

    if trunk > -1:
        if version[0] > -1:
            major, minor = trunk_major_minor(version[1], version[2])
            release = ''
            extra = '-alpha'

    elif version[0] > -1:
        major, minor, release, extra \
            = branch_release_extra(version[1], version[2], version[3],
                                   release[1], release[2], release[3])

    return major, minor, release, extra

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

def branch_release_extra(cur_major, cur_minor, cur_release,
                         prev_major, prev_minor, prev_release):
    """
    For branch, determine next version based on released versions
    """

    major = cur_major
    minor = cur_minor
    release = cur_release
    extra = ''

    # Same version, next release

    if major == prev_major and minor == prev_minor:
        release = prev_release
        if cur_release > prev_release:
            extra = '-patch'

    # Different version, or, no previous release

    else:
        try:
            if int(cur_release) == 0:
                release = ''
                extra = '-beta'
        except Exception:
            pass

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def svn_version_from_branch(branch):
    """
    Determine branch version number from Subversion branch
    """

    major = ''
    minor = ''

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

    major = ''
    minor = ''
    release = ''
    extra = ''

    if tag[0] == 'V':
        v = tag[1:].split('_', 3)
        if len(v) >= 1:
            major = v[0]
        if len(v) >= 2:
            minor = v[1]
        if len(v) >= 3:
            if v[2] > '-1':
                release = v[2]
            elif v[2][0:3] == '-rc':
                extra += v[3]
        if len(v) >= 4:
            extra += v[3]

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def svn_version(srcdir, defaults):
    """
    Determine version information from Subversion.
    """

    # Use default info to pre-determine next version

    major, minor, release, extra = defaults

    import subprocess

    # Use the XML output from Subversion, to avoid formatting/language
    # environment parsing issues

    try:
        from xml.dom import minidom
    except Exception:
        return major, minor, release, extra, revision

    # Get base info from 'svn info'

    url = None
    revision = ''

    cmd = ['svn', 'info', '--xml', srcdir]
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
    output = p.communicate()
    if p.returncode != 0:
        return major, minor, release, extra, revision
    else:
        try:
            svn_info = minidom.parseString(output[0]).documentElement
            url = svn_info.getElementsByTagName("url").item(0).firstChild.data
            entry_node = svn_info.getElementsByTagName("entry")[0]
            revision = entry_node.getAttribute("revision")
        except Exception:
            pass

    if not url:
        return major, minor, release, extra, revision

    # If we are in trunk, use previous local info to determine next version

    # If we are in a branch, use tags to determine next version

    if url.find('/branches/') > -1:
        url_id_0 = url.find('/branches/')
        url_id_1 = url_id_0 + len('/branches/')
        url_root = url[:url_id_0]

        if url[url_id_1:url_id_1+7] == 'Version':
            major, minor = svn_version_from_branch(url[url_id_1:])

        cmd = ['svn', 'ls', '--xml', url_root + '/tags']

        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
        output = p.communicate()
        if p.returncode != 0:
            return major, minor, release, extra, revision
        else:
            try:
                tag_ref = 'V' + str(major) + '_' + str(minor) + '_'
                tag_ref_l = len(tag_ref)
                name_max = tag_ref + '-1'
                release = ''
                svn_info = minidom.parseString(output[0]).documentElement
                for tag in svn_info.getElementsByTagName("name"):
                    name = tag.firstChild.data
                    if (name[0:tag_ref_l] == tag_ref and name > name_max
                        and name.find(name_max) < 0):
                        rev_node = tag.parentNode.getElementsByTagName("commit")[0]
                        tag_rev = rev_node.getAttribute("revision")
                        if (tag_rev == revision):
                            release = tag[tag_ref_l:]
                            break
                        elif tag_rev < revision:
                            name_max = name
                major, minor, release, extra = svn_version_from_tag(name_max)
                if not release:
                    extra += '-beta'
                if tag_rev == revision:
                    revision = ''
            except Exception:
                pass

    # If we are at a tag, simply use if

    elif url.find('/tags/') > -1:
        url_id_0 = url.find('/tags/')
        url_id_1 = url_id_0 + len('/tags/')
        tag = url[url_id_1:]
        major, minor, release, extra = svn_version_from_tag(url[url_id_1:])
        revision = ''

    if revision:
        revision = '-r' + revision
    else:
        revision = ''

    return major, minor, release, extra, revision

#-------------------------------------------------------------------------------

def svn_version_is_modified(srcdir):
    """
    Determine if version is modifed relative to Subversion revision
    (requires a local operation only).
    """

    import subprocess

    p = subprocess.Popen(['svn', 'status', '--xml', srcdir],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
    output = p.communicate()

    if p.returncode == 0:
        n_entries = output[0].count('<entry')
        n_po_entries = output[0].count('.po"')
        if n_entries - n_po_entries > 0:
            return True

    return False

#-------------------------------------------------------------------------------

def git_svn_version(srcdir, url, revision, defaults):
    """
    Determine partial version information from Subversion.
    """

    # Use known info to pre-determine next version

    major, minor, release, extra = defaults

    import subprocess

    # If we are in a branch, use tags if available to determine next version

    if url.find('/branches/') > -1:
        url_id_0 = url.find('/branches/')
        url_id_1 = url_id_0 + len('/branches/')
        url_root = url[:url_id_0]

        if url[url_id_1:url_id_1+7] == 'Version':
            major, minor = svn_version_from_branch(url[url_id_1:])
            revision = '-svn' + revision

    # If we are at a tag, simply use if

    elif url.find('/tags/') > -1:
        url_id_0 = url.find('/tags/')
        url_id_1 = url_id_0 + len('/tags/')
        tag = url[url_id_1:]
        major, minor, release, extra = svn_version_from_tag(url[url_id_1:])
        revision = ''

    # If we are in trunk, branch data might not be available to determine
    # next version, so use defaults

    else:
        revision = '-svn' + revision

    return major, minor, release, extra, revision

#-------------------------------------------------------------------------------

def git_version(srcdir, defaults):
    """
    Determine version information from Git.
    """

    # Use known info to pre-determine next version

    major, minor, release, extra = defaults

    import subprocess

    # Get base info from 'svn info'

    url = None
    revision = ''
    head_is_svn = False
    newline_failed = False

    cmd = ['git', '--no-pager', '--git-dir='+os.path.join(srcdir, '.git'),
           '--work-tree='+srcdir, 'log', '-n', '200']

    # Universal newlines options may fail in some cases with
    # Python 3, so check for newlines using '\n' and '\\n' later.

    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=False)
    output = p.communicate()
    if p.returncode != 0:
        return major, minor, release, extra
    else:
        o0 = str(output[0])
        svn_id = o0.find('git-svn-id:')
        if svn_id > 1:
            svn_rev = o0[svn_id:].split(' ', 2)[1]
            svn_url, svn_revision = svn_rev.split('@')
            major, minor, release, extra, revision \
                = git_svn_version(srcdir, svn_url, svn_revision, defaults)
            # Try to count commits before the one with svn
            svn_log = o0[:svn_id].split()
            c1 = svn_log.count('commit')
            c2 = svn_log.count('Author:')
            if c1 == 1 and c2 == 1:
                head_is_svn = True

    if not head_is_svn:
        commit_id = o0.find('commit')
        commit_rev = o0[commit_id:].split(' ', 2)[1].split('\n')[0]
        commit_rev = commit_rev.split('\\n')[0] # newline detection robustness
        revision += '-git' + commit_rev

    return major, minor, release, extra, revision

#-------------------------------------------------------------------------------

def git_version_is_modified(srcdir):
    """
    Determine if version is modifed relative to Git commit
    (requires a local operation only).
    """

    import subprocess

    cmd = ['git', '--no-pager', '--git-dir='+os.path.join(srcdir, '.git'),
           '--work-tree='+srcdir, 'status', '-s', '-uno']
    p = subprocess.Popen(cmd,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         universal_newlines=True)
    output = p.communicate()
    o0 = str(output[0])
    if p.returncode == 0:
        changes = o0.split('\n')
        n_changes = len(changes) - 1
        for c in changes:
            if c[-3:] == '.po':
                n_changes -= 1
        if n_changes > 0:
            return True

    return False

#-------------------------------------------------------------------------------

def build_version_string(major, minor, release, extra, revision = None):
    """
    Build version string from substrings.
    """

    version = major + '.' + minor
    if release:
        version += '.' + release
    if extra:
        version += extra
    if revision:
        version += revision
    version = version.replace('\r', '')
    return version

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    srcdir = os.path.split(os.path.abspath(os.path.dirname(sys.argv[0])))[0]

    # Default command: return simple string based on local information

    if len(sys.argv) == 1:
        major, minor, release, extra = version_from_news(srcdir)
        sys.stdout.write(build_version_string(major, minor, release, extra))

    elif sys.argv[1] == '--short':
        major, minor, release, extra = version_from_news(srcdir)
        sys.stdout.write(build_version_string(major, minor, '', extra))

    # Call tools

    elif sys.argv[1] in ['--full', '--revision-only', '--verbose']:

        defaults = version_from_news(srcdir)
        release = ''
        modified = ''
        revision = ''
        extra = ''

        if os.path.isdir(os.path.join(srcdir, ".svn")):
            major, minor, release, extra, revision \
                = svn_version(srcdir, defaults)
            modified = svn_version_is_modified(srcdir)
        elif os.path.isdir(os.path.join(srcdir, ".git")):
            major, minor, release, extra, revision \
                = git_version(srcdir, defaults)
            modified = git_version_is_modified(srcdir)
        elif os.path.isfile(os.path.join(srcdir, "build-aux", "version")):
            f = open(os.path.join(srcdir, "build-aux", "version"))
            lines = f.readlines()
            f.close()
            v = lines[0].strip()
            for e in ['-alpha', '-beta', '-rc']:
                i = v.find(e)
                if i > -1:
                    j = v[i+1:].find('-')
                    if j < 0:
                        j = v[i+1:].find('_')
                    if j > -1:
                        extra = v[i:i+j+1]
                    else:
                        extra = v[i:]
                    v = v.replace(extra, '')
            vl = v.split(".")
            l = len(vl)
            vt = vl.pop(l-1)
            vl += vt.split("-", 1)
            major = vl[0]
            minor = vl[1]
            if l > 2:
                release = vl[2]
            if len(vl) > l:
                revision = '-' + vl[l]
        else:
            major, minor, release, extra = version_from_news(srcdir)

        if sys.argv[1] == '--full':
            if modified:
                revision += '-m'
            sys.stdout.write(build_version_string(major, minor, release, extra, revision))
        elif sys.argv[1] == '--revision-only':
            if modified:
                revision += '-m'
            sys.stdout.write(revision[1:]) # remove first '-' character
        elif sys.argv[1] == '--verbose':
            print('major:    ' + major)
            print('minor:    ' + minor)
            print('release:  ' + release)
            print('extra:    ' + extra)
            print('revision: ' + revision[1:])
            print('modified: ' + str(modified))

    else:
        usage = \
            """Usage: %(prog)s [option]

Default:              print version string

Options:
  --full              print full version string with available version control info
  --short             print short version string with major, minor, and extras
  --revision-only     print version string with version control info only
  --verbose           print fields on different lines

Options:
  -h, --help  show this help message and exit"""

        print(usage % {'prog':sys.argv[0]})
        sys.exit(1)

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
