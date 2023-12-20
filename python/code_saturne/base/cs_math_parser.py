#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2023 EDF S.A.
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

import os
import re

#===============================================================================
# Utility functions
#===============================================================================

def create_req_field(name, dim=0):

    r = {'name':name,
         'dim':dim,
         'components':[]}

    return r

#-------------------------------------------------------------------------------

def rfield_add_comp(rf, c):

    rf['components'].append(c)
    rf['dim'] += 1

#-------------------------------------------------------------------------------

def split_req_components(req_list):
    """
    Look at a list of field names used in the formula.
    Check if its a component (f[X], f[XY], ..) or a field (f).
    return a list with it
    """
    req_fields = []
    for r in req_list:

        rf = r
        if bool(re.search('\[[A-Za-z0-9]\]', r)):
            rf = re.sub('\[[A-Za-z0-9]\]', '', r)
        elif bool(re.search('\[[A-Za-z0-9][A-Za-z0-9]\]', r)):
            rf = re.sub('\[[A-Za-z0-9][A-Za-z0-9]\]', '', r)

        if rf == r:
            req_fields.append(create_req_field(r,1))

        else:
            new_field = True
            for f in req_fields:
                if f['name'] == rf:
                    rfield_add_comp(f,r)
                    new_field = False
                    break

            if new_field:
                req_fields.append(create_req_field(rf))
                rfield_add_comp(req_fields[-1], r)

    return req_fields

#-------------------------------------------------------------------------------

def get_req_field_info(req_fields, r):

    for i in range(len(req_fields)):
        if req_fields[i]['dim'] == 1:
            if req_fields[i]['name'] == r:
                return i, -1, req_fields[i]['dim']

        else:
            if r in req_fields[i]['components']:
                return i, req_fields[i]['components'].index(r), req_fields[i]['dim']

    return None, None, None

#-------------------------------------------------------------------------------

def dump_req_fields(req_fields):

    print("===========================")
    for f in req_fields:
        print("Name: %s" % (f['name']))
        print("Dim: %d" % (f['dim']))
        print(f['components'])
        print("===========================")

#===============================================================================
# Mathematical expressions parser
#===============================================================================

class cs_math_parser:
    """
    Class for a mathematical formula parser.
    """

    #---------------------------------------------------------------------------

    def __init__(self):
        # Nothing to do

        return

    #---------------------------------------------------------------------------

    def find_c_comment_close(self, l, c_id, quotes):
        """
        Return index in given string of closing C-style comment,
        or -1 if not found after given c_id column.
        Escape characters or strings are handled.
        """

        w = len(l)

        while c_id < w:

            # Quoting (highest priority);
            # left inside expressions, but separators in quoted
            # regions are not considered as separators.

            if l[c_id] == "\\": # use 'e' to represent escape character
                if quotes[0] == "e":
                    quotes.pop(0)
                else:
                    quotes.insert(0, "e")
            elif l[c_id] == "'":
                if quotes[0] == "'":
                    quotes.pop(0)
                else:
                    quotes.insert(0, "'")
            elif l[c_id] == '"':
                if quotes[0] == '"':
                    quotes.pop(0)
                else:
                    quotes.insert(0, "'")

            elif quotes[0] == ' ': # found
                if l[c_id] == '*':
                    if c_id+1 < w:
                        if l[c_id+1] == '/':
                            return c_id

            c_id += 1

        return -1

    #---------------------------------------------------------------------------

    def separate_segments(self, lines):
        """
        Separate segments based on expected separators.
        This stage is not a parser, but simply splits lines into segments,
        separating comments (allowing for Python, C, and C++ style comments),
        and checking for quoting or escape characters.
        Returns a list of tuples containing the segments, and matching
        start line and column indexes in the original expression.
        """

        whitespace = (' ', '\t', '\n', '\l')
        separators = ('{', '}', ';')
        segments = []

        in_multiline_comment = False
        quotes = [' ']

        # Loop on lines

        l_id = 0
        for l in lines:

            w = len(l)

            # Loop on columns

            s_id = 0
            c_id = 0
            while c_id < w:

                # Quoting (highest priority);
                # left inside segments, but separators in quoted
                # regions are not considered as separators.

                if l[c_id] == "\\": # use 'e' to represent escape character
                    if quotes[0] == "e":
                        quotes.pop(0)
                    else:
                        quotes.insert(0, "e")
                elif l[c_id] == "'":
                    if quotes[0] == "'":
                        quotes.pop(0)
                    else:
                        quotes.insert(0, "'")
                elif l[c_id] == '"':
                    if quotes[0] == '"':
                        quotes.pop(0)
                    else:
                        quotes.insert(0, "'")

                if quotes[0] != ' ':
                    if quotes[0] == 'e' and l[c_id] in whitespace:
                        # escape character may be a line continuation character
                        # in this case; handle it like a separator
                        segments.append(('\\', l_id, c_id))
                        c_id += 1
                        s_id = c_id
                    continue

                # In multiline C-style comment
                # (transform to '#' for easier testing)

                elif in_multiline_comment:
                    j = self.find_c_comment_close(l, c_id, quotes)
                    if j >= 0: # on same line
                        j += 2
                        in_multiline_comment = False
                    else:
                        j = w
                    segments.append(('# ' + l[c_id:j].strip(), l_id, c_id))
                    c_id = j
                    s_id = c_id

                # Whitespace (handle here rather than using strip()
                # type functions to keep track of expression start columns

                elif l[c_id] in whitespace:
                    if s_id == c_id:
                        s_id += 1
                    c_id += 1
                    continue

                # Comments (allow C/C++ style, transform all to '#'
                # for easier testing)

                elif l[c_id] == '#':
                    e = l[s_id:c_id].strip()
                    if len(e):
                        segments.append((e, l_id, s_id))
                    segments.append((l[c_id:].strip(), l_id, c_id))
                    c_id = w
                    s_id = c_id

                elif l[c_id:c_id+2] in ('//', '/*'):
                    e = l[s_id:c_id].strip()
                    if len(e):
                        segments.append((e, l_id, s_id))
                    if l[c_id:c_id+2] == '//':
                        segments.append(('# ' + l[c_id+2:].strip(),
                                         l_id, c_id))
                        c_id = w
                        s_id = c_id
                    else:
                        j = self.find_c_comment_close(l, c_id+2, quotes)
                        if j >= 0: # on same line
                            segments.append(('# ' + l[c_id+2:j].strip(),
                                             l_id, c_id))
                            j += 2
                        else:
                            j = w
                            segments.append(('# ' + l[c_id+2:j].strip(),
                                             l_id, c_id))
                            in_multiline_comment = True
                        c_id = j
                        s_id = c_id

                else:

                    if l[c_id] in separators:
                        e = l[s_id:c_id].strip()
                        if len(e):
                            segments.append((e, l_id, s_id))
                        segments.append((l[c_id:c_id+1], l_id, c_id))
                        c_id += 1
                        s_id = c_id
                    else:
                        c_id += 1

            # End of loop on line:

            if s_id < c_id:
                e = l[s_id:c_id].strip()
                if len(e):
                    segments.append((e, l_id, s_id))

            l_id += 1

        return segments

    #---------------------------------------------------------------------------

    def parse_parentheses(self, line):

        istart = []
        d = {}

        for i, c in enumerate(line):
            if c == '(':
                istart.append(i)
            if c == ')':
                try:
                    d[istart.pop()] = i
                except IndexError:
                    print('There are too many closing parentheses')

        if istart:
            print('A closing parenthese is missing!')

        return d

    #---------------------------------------------------------------------------

    def get_start_lc(self, expr):
        """
        Return start line and column for a given expression
        """
        if isinstance(expr, (list,)):
            return self.get_start_lc(expr[0])

        else:
            return expr[1], expr[2]

    #---------------------------------------------------------------------------

    def recurse_expressions_syntax(self, expressions):
        """
        Recursively Update expressions
        """

        new_exp = []
        skip_to = 0

        for i, e in enumerate(expressions):

            if i < skip_to:
                continue

            if isinstance(e, (list,)):
                new_exp.append(self.recurse_expressions_syntax(e))

            else:

                # Translate "power" syntax
                if e[0] in ('^', '**') and i > 0:
                    valid = True
                    x = new_exp.pop()
                    y = None
                    try:
                        y = expressions[i+1]
                    except Exception:
                        valid = False
                        new_exp.append(x) # replace after pop() above
                    if valid:
                        sub_exp = []
                        li, ci = self.get_start_lc(x)
                        sub_exp.append(('(', li, ci))
                        sub_exp.append(x)
                        if y[0] in ('2', '3', '4'):
                            new_exp.append(('cs_math_pow'+y[0], li, ci))
                        else:
                            new_exp.append(('pow', li, ci))
                            sub_exp.append((',', li, ci))
                            li, ci = self.get_start_lc(x)
                            sub_exp.append(y)
                        sub_exp.append((')', li, ci))
                        new_exp.append((sub_exp, li, ci))
                        skip_to = i+2

                else:
                    new_exp.append(e)

        return new_exp

    #---------------------------------------------------------------------------

    def rename_math_functions(self, expressions):
        """
        Rename mathematical functions using the internal functions of
        code_saturne or standard math library.
        """
        _cs_math_internal_name = {'abs':'cs_math_fabs',
                                  'min':'cs_math_fmin',
                                  'max':'cs_math_fmax',
                                  'mod':'fmod',
                                  'square_norm':'cs_math_3_square_norm'}

        new_exp = []

        for e in expressions:

            if isinstance(e, (list, )):
                new_exp.append(self.rename_math_functions(e))

            else:
                if e[0] in _cs_math_internal_name.keys():
                    li, ci = self.get_start_lc(e)
                    en = _cs_math_internal_name[e[0]]
                    new_exp.append((en, li, ci))

                else:
                    new_exp.append(e)

        return new_exp

    #---------------------------------------------------------------------------

    def rebuild_text(self, expressions, comments,
                     level=0, s_line=0, s_col=0, t_prev=''):
        """
        Rebuild source code from expressions and comments.
        Reinsert comments at recorded lines in expressions.
        Comments are never inserted before an active token on a given
        line, event if this was the case in the original expression,
        both for simplification and because it is not recommeneded.
        """

        text = ''

        new_exp = []

        # operators adding spacing or not (minimalist prettyfier).
        spacing_operators = ('+', '-', '*', '%', '/', '=', '>', '<',
                             '==', '>=', '<=', 'if', 'then', 'else')
        no_spacing_operators = ('^', '**')

        for i, e in enumerate(expressions):

            li, ci = self.get_start_lc(e)

            # Restore comments from previous lines

            # column info does not need to be fully updated here,
            # as comments are always added at the end of a line,
            # and new lines reset column info.

            if comments:
                line_ref = comments[0][1]
                while comments:
                    if comments[0][1] >= li:
                        break
                    c = comments.pop(0)
                    while s_line < c[1]:
                        text += '\n'
                        s_line += 1
                        s_col = 0
                    if s_col > 0: # add space to nonempty column
                        text += ' '
                    if s_col < c[2]:
                        for j in range(c[2]):
                            text += ' '
                    text += '//' + c[0][1:]

            # Recursive handling of code

            if isinstance(e, list):
                sub_text, comments, e_line, e_col \
                    = self.rebuild_text(e, comments, level+1,
                                        s_line, s_col, t_prev)
                text += sub_text
                t_prev = sub_text[-1:]
                if e_line > s_line:
                    s_col = 0
                else:
                    s_col = e_col
                s_line = e_line
            elif isinstance(e[0], list):
                sub_text, comments, e_line, e_col \
                    = self.rebuild_text(e[0], comments, level+1,
                                        s_line, s_col, t_prev)
                text += sub_text
                t_prev = sub_text[-1:]
                if e_line > s_line:
                    s_col = 0
                else:
                    s_col = e_col
                s_line = e_line
            else:

                if s_line < li:
                    text += '\n'
                    s_col = 0
                    for i in range(ci):
                        text += ' '
                        s_col += 1

                else:      # Try to put spaces in recommended places
                    add_space = True
                    if t_prev in ('(', '[', '{', '', '% (int)'):
                        add_space = False
                    elif e[0] in (';', ':', ',', ')', ']', '}') \
                       or e[0] in no_spacing_operators:
                        add_space = False
                    elif e[0] in ('(', '['):
                        if t_prev not in spacing_operators:
                            add_space = False
                    if add_space:
                        text += ' '
                        s_col += 1

                # Add actual token

                text += e[0]
                t_prev = e[0]
                s_line = li
                s_col += len(e[0])

        # Restore comments after code

        if len(comments) > 0:
            if comments[0][1]:
                while comments[0][1] < e[1]:
                    c = comments.pop(0)
                    line_cur = c[1]
                    while line_cur < line_ref:
                        text += '\n'
                        line_cur += 1
                        for j in range(c[2]):
                            text += ' '
                            text += '//' + c[0][1:]
                        line_ref = line_cur
                    text += '\n\n'

        return text, comments, s_line, s_col

    #---------------------------------------------------------------------------

    def tokenize(self, segments):
        """
        Tokenize segments and separate comments.
        """

        whitespace = (' ', '\t', '\n', '\l', '\r')
        sep2 = ('<=', '>=', '!=', '==', '||', '&&', '+=', '-=', '*=', '/=', '**')
        sep1 = ('=', '(', ')', ';', ',', ':', '{', '}',
                '+', '-', '*', '/', '<', '>',  '^', '%', '!', '?')
        digits_p = ('.', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9')

        tokens = []
        comments = []

        for s in segments:

            prv = ' '
            s_id = 0
            s0 = s[0]
            if (s0[0] == '#'):
                comments.append(s)
                continue

            for i, c in enumerate(s0):
                if s_id > i:
                    continue
                elif c in whitespace:
                    if (not prv in whitespace) and (s_id < i):
                        tokens.append((s0[s_id:i], s[1], s[2]+s_id))
                    s_id = i+1
                elif s0[i:i+2] in sep2:
                    if (not prv in whitespace) and (s_id < i):
                        tokens.append((s0[s_id:i], s[1], s[2]+s_id))
                    tokens.append((s0[i:i+2], s[1], s[2]+i))
                    s_id = i+2
                elif c in sep1:
                    # special case: e+ or e- might not be a separator
                    is_exp = False
                    if c in ('+', '-'):
                        if s0[i-1:i+1] in ('e+', 'e-', 'E+', 'E-'):
                            if s0[i-2:i-1] in digits_p and s0[i+1:i+2] in digits_p:
                                is_exp = True
                    if not is_exp:
                        if (not prv in whitespace) and (s_id < i):
                            tokens.append((s0[s_id:i], s[1], s[2]+s_id))
                        tokens.append((s0[i:i+1], s[1], s[2]+i))
                        s_id = i+1
                prv = s0[i]
            r = s0[s_id:]
            if len(r) > 0:
                tokens.append((r, s[1], s[2]+s_id))

        return tokens, comments

    #---------------------------------------------------------------------------

    def build_expressions(self, exp_lines, tokens):
        """
        Organize expressions as lists of subexpressions based on levels
        """

        # Now we have a fully tokenized expression, we can buid a list of assignments

        opening_tokens = ('(', '{', '[')
        closing_tokens = (')', '}', ']')
        open_match = {')': '(',
                      '}': '{',
                      ']': '['}
        close_match = {'(': ')',
                      '{': '}',
                      '[': ']'}

        parent = []
        current = []

        rvalues = []
        expression = []

        previous = None
        level_open = []

        # Build levels based on parenthesis and braces (check consistency)

        match_error = None

        for t in tokens:

            if t[0] in opening_tokens:
                level_open.append(t)
                parent.append(current)
                current = []
                current.append(t)
            elif t[0] in closing_tokens:
                match_open = False
                if level_open:
                    t_open = level_open.pop()
                    if t_open[0] == open_match[t[0]]:
                        match_open = True
                        current.append(t)
                        sub = current
                        current = parent.pop()
                        current.append(sub)
                if not match_open:
                    match_error = (t[0], t[1], t[2], open_match[t[0]])
                    break
            else:
                current.append(t)

        if level_open:
            t = level_open.pop()
            match_error = (t[0], t[1], t[2], close_match[t[0]])

        if match_error:
            err_msg = []
            n_lines = 7
            j = match_error[1]  # error line
            for i in range(n_lines):
                if i + j + 1 - n_lines > -1:
                    err_msg.append(exp_lines[i + j + 1 - n_lines])
            c_line = ''
            for k in range(match_error[2]):  # error column
                c_line += ' '
            c_line += '^'
            err_msg.append(c_line)
            fmt = "Error: '{0}' at line {1} and column {2} does not have a matching '{3}'"
            err_msg.append(fmt.format(match_error[0], match_error[1],
                                      match_error[2], match_error[3]))
            err_msg.append('')
            for i in range(n_lines):
                if j+i < len(exp_lines):
                    err_msg.append(exp_lines[j-i])

            msg = ''
            for l in err_msg:
                msg += l + '\n'

            raise Exception(msg)

        return current

    #---------------------------------------------------------------------------

    def split_for_assignment(self, l):
        """
        Check for assignemnt (separating along = but not >=, <=)
        """
        lf = []
        c_idx = 0
        s_idx = 0
        e_idx = len(l)
        for c_idx in range(e_idx):
            if l[c_idx] == '=':
                if c_idx > 0 and l[c_idx-1] not in ('>', '<'):
                    lf.append(l[s_idx:c_idx])
                    s_idx = c_idx+1
        if s_idx <= e_idx:
            lf.append(l[s_idx:e_idx])

        return lf

    #---------------------------------------------------------------------------

    def parse_expression(self, expression, req, known_symbols,
                         func_type, glob_tokens, loop_tokens,
                         need_for_loop):
        """
        Parse an expression and return the corresponding C code, as well as
        the initialization block which needs to be used.
        """

        usr_defs = []
        usr_code = []

        nreq = len(req)
        req_fields = None
        if func_type == "vol":
            req_fields = split_req_components(req)

        # Parse the Mathematical expression and generate the C block code
        exp_lines = expression.split("\n")
        segments = self.separate_segments(exp_lines)
        tokens, comments = self.tokenize(segments)

        for t in tokens:
            # In the case of components, ensure the main field is added
            tklist = [t[0]]
            if (bool(re.search('\[[0-9]\]', t[0]))):
                tklist.append(re.sub('\[[0-9]\]', '', t[0]))

            for tk in tklist:
                if tk not in known_symbols:
                    # We use a double if and not if/else because some symbols
                    # may be present in both lists
                    if tk in glob_tokens.keys():
                        usr_defs.append(glob_tokens[tk] + '\n')
                        known_symbols.append(tk)
                    if tk in loop_tokens.keys():
                        usr_code.append(loop_tokens[tk] + '\n')
                        if tk not in known_symbols:
                            known_symbols.append(tk)

                        # For momentum source terms, check for velocity
                        if func_type == "src" and tk in ['u','v','w']:
                            if 'velocity' not in known_symbols:
                                if 'velocity' in glob_tokens:
                                    known_symbols.append('velocity')
                                    usr_defs.append(glob_tokens['velocity']+'\n')

        #-------------------------

        if len(usr_defs) > 0:
            usr_defs.append('\n')

        if len(usr_code) > 0:
            usr_code.append('\n')

        for t_i, t in enumerate(tokens):
            tk = t[0]
            # Check for assignments:
            if tk == "=" and t_i > 0:
                tk0 = tokens[t_i-1][0]
                if tk0 not in known_symbols:
                    usr_defs.append('cs_real_t %s = -1.;\n' % tk0)
                    known_symbols.append(tk0)


        req_to_replace = [elt for elt in req]
        for t_i, t in enumerate(tokens):
            tk = t[0]
            new_v = None
            if tk in req:
                if func_type == 'vol':
                    fid, fcomp, fdim = get_req_field_info(req_fields, tk)
                    if fid is None:
                        raise Exception("Uknown field: %s" %(tk))

                    if fcomp < 0:
                        new_v = new_v = 'fvals[%d][c_id]' % (fid)
                    else:
                        new_v = 'fvals[%d][c_id*%d + %d]' % (fid, fdim, fcomp)

                elif func_type == 'bnd':
                    ir = req.index(tk)
                    if need_for_loop:
                        new_v = 'retvals[%d * n_elts + e_id]' % (ir)
                    else:
                        new_v = 'retvals[%d]' % (ir)

                elif func_type in ['src', 'ini']:
                    if nreq > 1:
                        ir = req.index(tk)
                        new_v = 'retvals[%d * e_id + %d]' % (nreq, ir)
                    else:
                        new_v = 'retvals[e_id]'

                elif func_type == 'ibm':
                    new_v = '*ipenal'

                elif func_type == 'pca':
                    if nreq > 1:
                        ir = req.index(tk)
                        new_v = 'retvals[%d * e_id + %d]' % (nreq, ir)
                    else:
                        new_v = 'retvals[e_id]'

                if tk in req_to_replace:
                    req_to_replace.remove(tk)

            if new_v != None:
                tokens[t_i] = (new_v, t[1], t[2])

        if req_to_replace:
            msg="The following variables were not assigned in a formula:"
            for r in req_to_replace:
                msg+="\n"
                msg+=str(r)
            raise Exception(msg)
        #-------------------------

        tokens = self.rename_math_functions(tokens)
        tokens = self.build_expressions(exp_lines, tokens)
        tokens = self.recurse_expressions_syntax(tokens)

        #-------------------------

        # Rebuild lines
        new_text = self.rebuild_text(tokens, comments)
        for line in new_text[0].split('\n'):
            usr_code.append(line + '\n')
        if len(new_text[1]) > 0:
            for c in new_text[1]:
                usr_code.append('//' + c[0][1:] + '\n')

        #-------------------------

        return usr_code, usr_defs

#-------------------------------------------------------------------------------
