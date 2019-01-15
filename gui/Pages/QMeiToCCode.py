import os, shutil

from code_saturne.Pages.FluidCharacteristicsView import FluidCharacteristicsView


#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

_cs_fdict = {'density':'cs_user_meg_density.c',
             'viscosity':'cs_user_meg_viscosity.c',
             'specific_heat':'cs_user_meg_cp.c',
             'thermal_cond':'cs_user_meg_thermal_conductivity.c'}

_ncfd_fdict = {}

class mei_to_c_interpreter:

    # -------------------------------
    def __init__(self, case):

        self.case = case

        # function name to file name dictionary
        self.fdict = {}
        if self.case['package'].name == 'code_saturne':
            self.fdict = _cs_fdict
        else:
            self.fdict = _ncfd_fdict

        self.funcs = {}
    # -------------------------------


    # -------------------------------
    def init_c_function(self,
                        expression,
                        required,
                        symbols,
                        scalars,
                        function_name):

        self.funcs[function_name] = {'exp':expression,
                                     'sym':symbols,
                                     'sca':scalars,
                                     'req':required}

    # -------------------------------

    # -------------------------------
    def write_c_function(self, function_name):

        # Check if function exists
        try:
            c_file_name = self.fdict[function_name]
        except:
            error_msg = "Function %s is uknown!" %(function_name)
            raise Exception(error_msg)

        if function_name not in self.funcs.keys():
            return

        expression = self.funcs[function_name]['exp']
        symbols    = self.funcs[function_name]['sym']
        scalars    = self.funcs[function_name]['sca']
        required   = self.funcs[function_name]['req']

        # Copy function file if needed
        fpath = os.path.join(self.case['case_path'], 'SRC', c_file_name);
        if not os.path.isfile(fpath):
            fexample = os.path.join(self.case['package'].dirs['pkgdatadir'][1],
                                    'user_meg',
                                    c_file_name)

            shutil.copy2(fexample, fpath)

        # Get user definitions and code
        usr_defs = ''
        usr_code = ''

        tab   = "  "
        ntabs = 1

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
            exp_str = exp_str.replace(s[0],'*'+s[0])
            known_symbols.append('*'+s[0]);
            usr_code += (ntabs+1)*tab + 'cs_real_t *'+s[0]+' = &values[c_id];\n'

        exp_lines = exp_str.split("\n")

        ntabs += 1
        if_loop = False
        for l in exp_lines:
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

        if if_loop and usr_code[-1] != '}':
            usr_code += '\n' + (ntabs-1)*tab + '}\n'
        else:
            usr_code += '\n'

        old_code = open(fpath, 'r').readlines()
        # Find positions to insert user defs and code
        usr_defs_idx1 = -1
        usr_defs_str1 = '/* START USER DEFS */'
        usr_defs_idx2 = -1
        usr_defs_str2 = '/* END USER DEFS */'
        usr_code_idx1 = -1
        usr_code_str1 = '/* START USER CODE */'
        usr_code_idx2 = -1
        usr_code_str2 = '/* END USER CODE */'

        for ii in range(len(old_code)):
            if usr_defs_str1 in old_code[ii]:
                usr_defs_idx1 = ii

            if usr_defs_str2 in old_code[ii]:
                usr_defs_idx2 = ii

            if usr_code_str1 in old_code[ii]:
                usr_code_idx1 = ii

            if usr_code_str2 in old_code[ii]:
                usr_code_idx2 = ii
                break

        new_code = open(fpath,'w')
        for ii in range(usr_defs_idx1+1):
            if 'return;' not in old_code[ii]:
                new_code.write(old_code[ii])

        new_code.write(usr_defs)

        for ii in range(usr_defs_idx2,usr_code_idx1+1):
            new_code.write(old_code[ii])

        new_code.write(usr_code)

        for ii in range(usr_code_idx2,len(old_code)):
            new_code.write(old_code[ii])
    # -------------------------------

    # -------------------------------
    def delete_c_function(self, function_name):

        # Check if function exists
        try:
            c_file_name = self.fdict[function_name]
        except:
            error_msg = "Function %s is uknown!" %(function_name)
            raise Exception(error_msg)

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
                self.init_c_function(exp, req, sym, sca, 'density')

            # Viscosity
            if fcv.mdl.getPropertyMode('molecular_viscosity') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaMuComponents()
                self.init_c_function(exp, req, sym, sca, 'viscosity')

            # Specific heat
            if fcv.mdl.getPropertyMode('specific_heat') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaCpComponents()
                self.init_c_function(exp, req, sym, sca, 'specific_heat')

            # Thermal conductivity
            if fcv.mdl.getPropertyMode('thermal_conductivity') == 'variable':
                exp, req, sca, sym, exa = fcv.getFormulaAlComponents()
                self.init_c_function(exp, req, sym, sca, 'thermal_cond')

        # Write the functions
        for key in self.fdict.keys():
            if key in self.funcs.keys():
                self.write_c_function(key)
            else:
                self.delete_c_function(key)
    # -------------------------------

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
