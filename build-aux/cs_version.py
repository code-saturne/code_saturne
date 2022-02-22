#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2021 EDF S.A.
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
# for versions either pulled by Git or simple release tarballs.
# As we need to maintain the NEWS file anyways to provide a clear
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

    master = -1
    version = [-1]
    release = [-1, '', '', '']

    f = open(os.path.join(srcdir, "NEWS.md"))
    lines = f.readlines()
    f.close()

    for i in range(1, len(lines)):

        # Release names are underlined by '---'

        if lines[i][0:7] == '-------':

            j = i-1
            l = lines[j]

            if l[0:7] == 'Master ' and master < 0:
                master = j

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

    if master > -1:
        if version[0] > -1:
            major, minor = master_major_minor(version[1], version[2])
            release = ''
            extra = '-alpha'

    elif version[0] > -1:
        while len(version) < 4:
            version.append('')
        major, minor, release, extra \
            = branch_release_extra(version[1], version[2], version[3],
                                   release[1], release[2], release[3])

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def master_major_minor(prev_major, prev_minor):
    """
    For master, determine next version based on previous branches
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
            if not cur_release or int(cur_release) == 0:
                release = ''
                extra = '-beta'
        except Exception:
            pass

    return major, minor, release, extra

#-------------------------------------------------------------------------------

def git_version(srcdir, defaults):
    """
    Determine version information from Git.
    """

    # Use known info to pre-determine next version

    major, minor, release, extra = defaults

    import subprocess

    # Get base info

    url = None
    revision = ''
    newline_failed = False

    cmd = ['git', '--no-pager', '--git-dir='+os.path.join(srcdir, '.git'),
           '--work-tree='+srcdir, 'log', '--pretty=format:%h','-n', '1']

    # Universal newlines options may fail in some cases with
    # Python 3, so check for newlines using '\n' and '\\n' later.

    try:
        p = subprocess.Popen(cmd,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=False)
        output = p.communicate()
        if p.returncode == 0:
            commit_rev = output[0]
            if sys.version[0] == '3':
                commit_rev = commit_rev.decode('utf-8')
            revision += '-' + commit_rev
    except Exception:
        pass

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
    try:
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
    except Exception:
        pass

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

def replace_in_file(src, dest, major, minor, release, extra, revision, modified):
    """
    Replace expressions in file based on version numbers.
    """

    f = open(src, 'r')
    lines = f.readlines()
    f.close()

    cs_version = build_version_string(major, minor, release, extra)
    cs_version_full = build_version_string(major, minor, release, extra,
                                           revision)
    cs_version_short = build_version_string(major, minor, '', extra)

    if modified:
        cs_version_full += '-m'
    cs_revision = revision[1:]
    if modified:
        revision += '-m'

    for i, l in enumerate(lines):
        l = l.replace("@cs_version@", cs_version)
        l = l.replace("@cs_version_full@", cs_version_full)
        l = l.replace("@cs_version_short@", cs_version_short)
        l = l.replace("@cs_revision@", cs_revision)
        lines[i] = l

    if not dest:
        for l in lines:
            sys.stdout.write(l)

    else:
        changed = True

        cmp_lines = []
        if os.path.isfile(dest):
            f = open(dest, 'r')
            cmp_lines = f.readlines()
            f.close()
            changed = False

        if len(cmp_lines) != len(lines):
            changed = True
        else:
            for i, l in enumerate(lines):
                if l != cmp_lines[i]:
                    changed = True
                    break

        if changed:
            f = open(dest, 'w')
            for l in lines:
                f.write(l)
            f.close()

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import sys

    srcdir = os.path.split(os.path.abspath(os.path.dirname(sys.argv[0])))[0]

    # Default command: return simple string based on local information

    n_args = len(sys.argv) - 1

    if n_args == 0:
        major, minor, release, extra = version_from_news(srcdir)
        sys.stdout.write(build_version_string(major, minor, release, extra))

    elif sys.argv[1] == '--short':
        major, minor, release, extra = version_from_news(srcdir)
        if not extra in ['-alpha', '-beta']:
            extra = ''
        sys.stdout.write(build_version_string(major, minor, '', extra))

    # Call tools

    elif sys.argv[1] in ['--full', '--revision-only', '--verbose',
                         '--replace']:

        defaults = version_from_news(srcdir)
        release = ''
        modified = ''
        revision = ''
        extra = ''

        if os.path.isdir(os.path.join(srcdir, ".git")):
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
        elif sys.argv[1] == '--replace' and n_args in (2, 3):
            src = sys.argv[2]
            dest = None
            if n_args == 3:
                dest = sys.argv[3]
            replace_in_file(src, dest, major, minor, release,
                            extra, revision, modified)

    else:
        usage = \
            """Usage: %(prog)s [option]

Default:              print version string

Options:
  --full                      print full version string with available
                              version control info
  --short                     print short version string with major, minor,
                              and extras
  --revision-only             print version string with version control info only
  --update <template> [path]  build or rebuild output file from template
  --verbose                   print fields on different lines

Options:
  -h, --help  show this help message and exit"""

        print(usage % {'prog':sys.argv[0]})
        sys.exit(1)

    sys.exit(0)

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
