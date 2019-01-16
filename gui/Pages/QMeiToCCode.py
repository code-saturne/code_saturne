import os, shutil

from code_saturne.Pages.FluidCharacteristicsView import FluidCharacteristicsView


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

_file_header = \
"""/*----------------------------------------------------------------------------*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

"""

_file_footer = \
"""}

END_C_DECLS

"""

_vol_function_header = \
"""void
cs_meg_volume_function(cs_field_t              *f,
                       const cs_volume_zone_t  *z)
{

"""


class mei_to_c_interpreter:

    # -------------------------------
    def __init__(self, case):

        self.case = case

        # function name to file name dictionary
        self.funcs = {}

        self.code_to_write = ""
    # -------------------------------


    # -------------------------------
    def init_c_block(self,
                     expression,
                     required,
                     symbols,
                     scalars,
                     name,
                     zone_name='all_cells'):

        if name in self.funcs.keys():
            if self.funcs[name][zone] == zone_name:
                msg = "Formula for variable %s in zone %s was allready defined:\n %s" \
                        % (name, zone_name, self.funcs[name][exp])
                raise Exception(msg)

        self.funcs[name] = {'exp':expression,
                            'sym':symbols,
                            'sca':scalars,
                            'req':required,
                            'zone':zone_name}

    # -------------------------------

    # -------------------------------
    def write_c_block(self, name):

        # Check if function exists
        if name not in self.funcs.keys():
            return

        expression = self.funcs[name]['exp']
        symbols    = self.funcs[name]['sym']
        scalars    = self.funcs[name]['sca']
        required   = self.funcs[name]['req']
        zone       = self.funcs[name]['zone']

        # Get user definitions and code
        usr_defs = ''
        usr_code = ''
        usr_blck = ''


        tab   = "  "
        ntabs = 2

        known_symbols = []
        coords = ['x', 'y', 'z']
        need_coords = False

        internal_fields = []

        for s in symbols:
            if s[0] == 'dt':
                usr_defs += ntabs*tab + 'cs_real_t dt = cs_glob_time_step->dt;\n'
                known_symbols.append(s[0])
            elif s[0] == 'time':
                usr_defs += ntabs*tab + 'cs_real_t time = cs_glob_time_step->t_cur;\n'
                known_symbols.append(s[0])
            elif s[0] == 'iter':
                usr_defs += ntabs*tab + 'int iter = cs_glob_time_step->nt_cur;\n'
                known_symbols.append(s[0])
            elif s[0] in coords:
                ic = coords.index(s[0])
                lxyz = 'cs_real_t %s = xyz[c_id][%s];\n' % (s[0], str(ic))
                usr_code += (ntabs+1)*tab + lxyz
                known_symbols.appned(s[0])
                need_coords = True
            elif s not in scalars:
                if len(s[1].split('=')) > 1:
                    sval = s[1].split('=')[-1]
                    usr_defs += ntabs*tab + 'cs_real_t '+s[0]+' = '+str(sval)+';\n'
                    known_symbols.append(s[0])
                else:
                    internal_fields.append((s[0], s[1].lower()))

        if need_coords:
            user_defs = ntabs*tab \
                      + 'cs_real_3_t xyz = (cs_real_3_t *)cs_glob_mesh_quantities->cell_cen;' \
                      + '\n\n' \
                      + user_defs

        for f in internal_fields:
            l = 'cs_real_t *%s_vals = cs_field_by_name("%s")->val;\n' % (f[0], f[1])
            usr_defs += ntabs*tab + l
            known_symbols.append(f[0])
            usr_code += (ntabs+1)*tab + 'cs_real_t %s = %s_vals[c_id];\n' % (f[0], f[0])
        usr_code += '\n'

        for s in scalars:
            sn = s[0]
            l = 'cs_real_t *%s_vals = cs_field_by_name("%s")->val;\n' % (sn,sn)
            usr_defs += ntabs*tab + l
            known_symbols.append(sn)
            usr_code += (ntabs+1)*tab + 'cs_real_t %s = %s_vals[c_id];\n' % (sn, sn)
        usr_code += '\n'

        exp_str = expression
        for s in required:
            known_symbols.append(s[0]);

        exp_lines = exp_str.split("\n")

        ntabs += 1
        if_loop = False
        for l in exp_lines:
            if len(l) > 0:
                lf = l.split('=')
                if len(lf) > 1 and lf[0].rstrip() not in known_symbols:
                    known_symbols.append(lf[0])
                    l = 'cs_real_t '+l

                if l[:2] == 'if' and l[-1] != '{':
                    l = l + ' {'
                    if_loop = True

                if 'else' in l:
                    if l[0] != '}':
                        l = '} '+l
                    if l[-1] != '{':
                        l = l + ' {'

                if l[0] == '}':
                    ntabs -= 1
                usr_code += '\n'
                usr_code += ntabs*tab + l
                if l[-1] == '{':
                    ntabs += 1
            else:
                usr_code += '\n'

        if if_loop and usr_code[-1] != '}':
            usr_code += '\n' + (ntabs-1)*tab + '}\n'
        else:
            usr_code += '\n'

        usr_code = usr_code.replace(name, 'f->val[c_id]')

        # Write the block
        usr_blck = tab + 'if (strcmp(f->name, "%s") == 0 && strcmp(z->name, "%s") == 0) {\n' \
                % (name, zone)

        usr_blck += '\n' + usr_defs + '\n'

        usr_blck += 2*tab + 'for (cs_lnum_t e_id; e_id < z->n_cells; e_id++) {\n'
        usr_blck += 3*tab + 'cs_lnum_t c_id = z->cell_ids[e_id];'

        usr_blck += usr_code
        usr_blck += 2*tab + '}\n'
        usr_blck += tab + '}\n'

        return usr_blck
    # -------------------------------

    # -------------------------------
    def delete_file(self, c_file_name):

        # Copy function file if needed
        fpath = os.path.join(self.case['case_path'], 'SRC', c_file_name);
        if os.path.isfile(fpath):
            os.remove(fpath)
    # -------------------------------

    # -------------------------------
    def save_all_functions(self, root):

        self.root = root
        if self.case['package'].name == 'code_saturne':
            from code_saturne.Pages.FluidCharacteristicsView import FluidCharacteristicsView
            fcv = FluidCharacteristicsView(root, self.case)

            # Density
            if fcv.mdl.getPropertyMode('density') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaRhoComponents()
                self.init_c_block(exp, req, sym, sca, 'density')

            # Viscosity
            if fcv.mdl.getPropertyMode('molecular_viscosity') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaMuComponents()
                self.init_c_block(exp, req, sym, sca, 'molecular_viscosity')

            # Specific heat
            if fcv.mdl.getPropertyMode('specific_heat') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaCpComponents()
                self.init_c_block(exp, req, sym, sca, 'specific_heat')

            # Thermal conductivity
            if fcv.mdl.getPropertyMode('thermal_conductivity') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaAlComponents()
                self.init_c_block(exp, req, sym, sca, 'thermal_conductivity')

        # Delete previous existing file
        file2write = 'cs_meg_volume_function.c'
        self.delete_file(file2write)

        # Write the functions if needed
        if len(self.funcs.keys()) > 0:
            code_to_write = _file_header + _vol_function_header
            for key in self.funcs.keys():
                m1 = '/* User defined formula for variable %s over zone %s  */' \
                        % (key, self.funcs[key]['zone'])
                m2 = '/* '+'-'*(len(m1)-6)+' */\n'

                code_to_write += "  " + m2
                code_to_write += "  " + m1 + '\n'
                code_to_write += self.write_c_block(key)
                code_to_write += "  " + m2 + '\n'

            code_to_write += _file_footer

            fpath = os.path.join(self.case['case_path'], 'SRC', file2write)
            new_file = open(fpath, 'w')
            new_file.write(code_to_write)
            new_file.close()
    # -------------------------------

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
