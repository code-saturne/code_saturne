#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_doc.py - Generate documentation from XML rule files


Generates Sphinx documentation for:
- TurbulenceRules.xml
- ThermalRules.xml
- TimeSteppingRules.xml
"""

from pathlib import Path
from datetime import datetime
import xml.etree.ElementTree as ET
import sys
import os

# Fix PYTHONPATH
script_dir = os.path.dirname(os.path.abspath(__file__))
site_packages = os.path.abspath(os.path.join(
    script_dir, '..', '..', '..', '..', '..',
    'arch', 'master.test.model', 'lib', 'python3.11', 'site-packages'))
sys.path.insert(0, site_packages)

DATA_DIR = os.path.abspath(os.path.join(script_dir, '..'))


# =============================================================================
# BASE GENERATOR
# =============================================================================

class BaseXMLDocGenerator:
    def __init__(self, xml_path, module_name):
        self.xml_path = xml_path
        self.module_name = module_name
        self.tree = ET.parse(xml_path)
        self.root = self.tree.getroot()

    def _write_file(self, path, content):
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"  -> {path.name}")

    def _get_typedef_values(self, name):
        defs = self.root.find("Definitions")
        if defs is None:
            return []
        typedef = defs.find(f"TypeDef[@name='{name}']")
        if typedef is None:
            return []
        return [v.text.strip() for v in typedef.findall("Value") if v.text]

    def generate(self, output_dir):
        raise NotImplementedError


# =============================================================================
# TURBULENCE GENERATOR
# =============================================================================

class TurbulenceDocGenerator(BaseXMLDocGenerator):

    def __init__(self, xml_path):
        super().__init__(xml_path, "Turbulence")
        self.valid_models = set(self._get_typedef_values("TurbulenceModel"))
        self.valid_variables = set(self._get_typedef_values("TurbulenceVariable"))
        self.valid_properties = set(self._get_typedef_values("TurbulenceProperty"))
        self._load_model_groups()
        self._load_model_requirements()

    def _load_model_groups(self):
        self.model_groups = {}
        groups_node = self.root.find("ModelGroups")
        if groups_node is None:
            return
        for group in groups_node.findall("Group"):
            name = group.get("Name")
            models = [m.text.strip() for m in group.findall("Model") if m.text]
            self.model_groups[name] = models

    def _load_model_requirements(self):
        self.model_requirements = {}
        for rule in self.root.findall(".//Rule[@Type='RequiredComponents']"):
            model = rule.get("Model")
            if model:
                vars_list = [{"name": v.get("Name")} for v in rule.findall("RequireVariable")]
                props = [{"name": p.get("Name")} for p in rule.findall("RequireProperty")]
                self.model_requirements[model] = {"variables": vars_list, "properties": props}

    def get_model_group(self, group_name):
        return self.model_groups.get(group_name, [])

    def is_model_in_group(self, model, group):
        return model in self.model_groups.get(group, [])

    def get_required_variables(self, model):
        return self.model_requirements.get(model, {}).get("variables", [])

    def get_required_properties(self, model):
        return self.model_requirements.get(model, {}).get("properties", [])

    def generate(self, output_dir="source/turbulence"):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        print("Generating TurbulenceRules.xml documentation...")
        self._gen_index(output_path)
        self._gen_overview(output_path)
        self._gen_models(output_path)
        self._gen_variables(output_path)
        self._gen_rules(output_path)
        self._gen_groups(output_path)
        self._gen_reference(output_path)
        print(f"  -> RST files in {output_path}\n")

    def _gen_index(self, p):
        content = f"""TurbulenceRules.xml -- Technical Reference
==========================================

Auto-generated documentation from ``TurbulenceRules.xml``

.. note::
   Generated on {datetime.now().strftime("%d/%m/%Y at %H:%M")}.

.. toctree::
   :maxdepth: 2
   :numbered:

   overview
   models
   variables
   rules
   groups
   reference

* :ref:`genindex`
* :ref:`search`
"""
        self._write_file(p / "index.rst", content)

    def _gen_overview(self, p):
        content = """Overview
========

``TurbulenceRules.xml`` is the single source of truth for turbulence
rules in Code_Saturne. It is read by:

* Le **kernel C++** via ``cs_turbulence_rules_manager``
* La **GUI Python** via ``TurbulenceRulesManager`` (ElementTree)

Main sections
--------------------

* ``Definitions`` : valid types (models, variables, properties)
* ``TurbulenceModels`` : numerical configuration for each model
* ``ModelGroups`` : RANS / LES / Rij / TwoEquations classification
* ``ValidationRules`` : GUI and kernel validation rules
* ``Defaults`` : default values
* ``Dependencies`` : compatibility with other physics modules
"""
        self._write_file(p / "overview.rst", content)

    def _gen_models(self, p):
        content = """Turbulence Models
=================

.. contents::
   :local:
   :depth: 2

"""
        for group_name in ["RANS", "LES", "RijModels", "TwoEquations"]:
            models = self.get_model_group(group_name)
            if models:
                content += f"\n{group_name} Models\n" + "-" * (len(group_name) + 7) + "\n\n"
                for m in sorted(models):
                    content += f"\n{m}\n" + "^" * len(m) + "\n\n"
                    vars_list = self.get_required_variables(m)
                    if vars_list:
                        content += "Required variables : " + ", ".join(f"``{v['name']}``" for v in vars_list) + "\n\n"
                    props = self.get_required_properties(m)
                    if props:
                        content += "Required properties : " + ", ".join(f"``{p['name']}``" for p in props) + "\n\n"
                    groups = [g for g in ["RANS", "LES", "RijModels", "TwoEquations"] if self.is_model_in_group(m, g)]
                    if groups:
                        content += f"Groups : {', '.join(groups)}\n\n"
        self._write_file(p / "models.rst", content)

    def _gen_variables(self, p):
        content = """Variables and Properties
=======================

Turbulence variables
---------------------

"""
        var_desc = {
            "k": "Turbulent kinetic energy (m2/s2)",
            "epsilon": "Dissipation rate (m2/s3)",
            "omega": "Specific dissipation rate (1/s)",
            "rij": "Reynolds stress tensor (m2/s2)",
            "nu_tilda": "Modified turbulent viscosity (m2/s)",
            "phi": "Elliptic variable v2f (m2/s2)",
            "alpha": "Elliptic variable EBRSM (-)",
        }
        for v in sorted(self.valid_variables):
            desc = var_desc.get(v, "Turbulence variable")
            content += f"* ``{v}`` : {desc}\n"
        content += "\n\nTurbulence properties\n----------------------\n\n"
        for prop in sorted(self.valid_properties):
            content += f"* ``{prop}``\n"
        self._write_file(p / "variables.rst", content)

    def _gen_rules(self, p):
        content = """Validation rules
====================

Kernel parameter ranges
---------------------------

.. list-table::
   :widths: 30 20 20 30
   :header-rows: 1

   * - Parameter
     - Min
     - Max
     - Description
   * - iwallf
     - 0
     - 8
     - Wall function type
   * - iclkep
     - 0
     - 2
     - k-epsilon clipping mode
   * - ikecou
     - 0
     - 2
     - k-epsilon coupling mode

LES constraints
---------------

* Wall functions forbidden
* Time scheme order 2 required
"""
        self._write_file(p / "rules.rst", content)

    def _gen_groups(self, p):
        content = """Model groups
==================

"""
        for group_name, models in self.model_groups.items():
            content += f"{group_name}\n" + "-" * len(group_name) + "\n\n"
            for m in sorted(models):
                content += f"* {m}\n"
            content += "\n"
        self._write_file(p / "groups.rst", content)

    def _load_wall_functions(self):
        """Lire les wall functions depuis TurbulenceRules.xml"""
        wf = {}
        for rule in self.root.findall(".//Rule[@Type='WallFunction']"):
            model = rule.get("Model")
            group = rule.get("Group")
            key = model if model else group
            allowed = [n.text for n in rule.findall("Allowed")]
            default = rule.findtext("Default", "-")
            wf[key] = {"allowed": allowed, "default": default}
        return wf

    def _get_wall_function(self, model_name, wf_rules):
        """Obtenir la wall function pour un modele (modele ou groupe)"""
        if model_name in wf_rules:
            return wf_rules[model_name]
        # Chercher dans les groupes
        for group in ["TwoEquations", "RijModels", "LES"]:
            if self.is_model_in_group(model_name, group) and group in wf_rules:
                return wf_rules[group]
        return {"allowed": ["-"], "default": "-"}

    def _gen_reference(self, p):
        wf_rules = self._load_wall_functions()

        content = """Reference table
====================

Models and their parameters
----------------------------

.. list-table::
   :widths: 25 15 25 10 25
   :header-rows: 1

   * - Model
     - Family
     - Required variables
     - Default WF
     - Allowed WF
"""
        for model in sorted(self.valid_models):
            famille = [g for g in ["RANS", "LES"] if self.is_model_in_group(model, g)]
            famille_str = ", ".join(famille) if famille else "-"
            vars_list = self.get_required_variables(model)
            vars_str = ", ".join(v['name'] for v in vars_list) if vars_list else "-"
            wf = self._get_wall_function(model, wf_rules)
            wf_default = wf["default"]
            wf_allowed = ", ".join(wf["allowed"]) if wf["allowed"] else "-"
            content += f"   * - {model}\n"
            content += f"     - {famille_str}\n"
            content += f"     - {vars_str}\n"
            content += f"     - {wf_default}\n"
            content += f"     - {wf_allowed}\n"

        content += """
Wall Functions — Value to description mapping
---------------------------------------------

.. list-table::
   :widths: 10 30 60
   :header-rows: 1

   * - Value
     - Name
     - Description
   * - 0
     - No wall function
     - Wall-resolved (y+ ~ 1). Required for v2f and some LES models.
   * - 2
     - Standard log law
     - Classical log law. Requires y+ ~ 30-300.
   * - 3
     - All y+ (two-layer)
     - Works for any mesh (low or high y+). Recommended by default.
   * - 4
     - Scalable wall function
     - Log law variant, robust for variable y+.
   * - 7
     - Two-scale model (EBRSM)
     - Specific to Rij-EBRSM model. Two-velocity-scale model.

.. note::
   The default value is automatically reset when the user
   changes the turbulence model in the GUI (wall function bug fix).
"""
        self._write_file(p / "reference.rst", content)


# =============================================================================
# THERMAL GENERATOR
# =============================================================================

class ThermalDocGenerator(BaseXMLDocGenerator):

    def __init__(self, xml_path):
        super().__init__(xml_path, "Thermal")
        self.thermal_variables = self._get_typedef_values("ThermalVariable")
        self.eos_values = self._get_typedef_values("EquationOfState")
        self.temp_scales = self._get_typedef_values("TemperatureScale")
        self.thermo_planes = self._get_typedef_values("ThermoPlane")

    def generate(self, output_dir="source/thermal"):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        print("Generating ThermalRules.xml documentation...")
        self._gen_index(output_path)
        self._gen_overview(output_path)
        self._gen_variables(output_path)
        self._gen_fields(output_path)
        self._gen_validation(output_path)
        self._gen_reference(output_path)
        print(f"  -> RST files in {output_path}\n")

    def _gen_index(self, p):
        content = f"""ThermalRules.xml -- Technical Reference
=======================================

Auto-generated documentation from ``ThermalRules.xml``

.. note::
   Generated on {datetime.now().strftime("%d/%m/%Y at %H:%M")}.

.. toctree::
   :maxdepth: 2
   :numbered:

   overview
   variables
   fields
   validation
   reference

* :ref:`genindex`
* :ref:`search`
"""
        self._write_file(p / "index.rst", content)

    def _gen_overview(self, p):
        content = """Overview
========

``ThermalRules.xml`` centralizes the thermal model rules of Code_Saturne.

Il est lu par :

* Le **kernel C++** via ``cs_thermal_rules_manager``
* La **GUI Python** via ``ThermalRulesManager`` (ElementTree)

Supported thermal variables
--------------------------------

"""
        for v in self.thermal_variables:
            content += f"* ``{v}``\n"
        content += "\n\nEquations of state\n----------------\n\n"
        for e in self.eos_values:
            content += f"* ``{e}``\n"
        content += "\n\nThermodynamic planes\n----------------------\n\n"
        for plane in self.thermo_planes:
            content += f"* ``{plane}``\n"
        self._write_file(p / "overview.rst", content)

    def _gen_variables(self, p):
        content = """Thermal variables
====================

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

   * - Variable
     - Diffusivity formula
     - Thermodynamic plane
     - is_temperature flag
"""
        diffusivity_map = {
            "temperature": ("lambda", "PT", "1"),
            "enthalpy": ("lambda/cp", "PH", "0"),
            "internal_energy": ("lambda/cp", "PH", "0"),
            "none": ("-", "-", "-"),
        }
        for v in self.thermal_variables:
            diff, plane, flag = diffusivity_map.get(v, ("-", "-", "-"))
            content += f"   * - ``{v}``\n     - {diff}\n     - {plane}\n     - {flag}\n"
        content += "\n\nTemperature scales\n-----------------------\n\n"
        scale_desc = {
            "none": "No active thermal model",
            "kelvin": "Temperature in Kelvin (K)",
            "celsius": "Temperature in Celsius (C)",
        }
        for s in self.temp_scales:
            desc = scale_desc.get(s, s)
            content += f"* ``{s}`` : {desc}\n"
        self._write_file(p / "variables.rst", content)

    def _gen_fields(self, p):
        content = """Fields created per thermal variable
====================================

Temperature
-----------

* ``temperature`` (variable, cellules)
* ``boundary_temperature`` (frontiere, faces de bord)

Enthalpy
---------

* ``enthalpy`` (variable, cellules)
* ``boundary_temperature`` (frontiere, faces de bord)
* ``thermal_diffusivity`` (propriete, cellules)

Internal energy
---------------

* ``total_energy`` (variable, cellules)
* ``temperature`` (propriete, cellules)
* ``boundary_temperature`` (frontiere, faces de bord)

Conditional fields
--------------------

**If** ``has_kinetic_st`` :

* ``kinetic_energy_thermal_st`` (propriete)
* ``rho_k_prev`` (interne)
* ``inner_face_velocity`` (faces internes, dim=3)
* ``boundary_face_velocity`` (faces de bord, dim=3)

**If** ``ieos == moist_air`` :

* ``yv`` (property, water vapor mass fraction)
"""
        self._write_file(p / "fields.rst", content)

    def _gen_validation(self, p):
        content = """Validation rules
====================

Parameter ``is_temperature``
----------------------------

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Min
     - Max
     - Description
   * - 0
     - 1
     - 0 = passive scalar/enthalpy, 1 = temperature

Parameter ``turbulent_schmidt``
-------------------------------

Must be strictly positive (> 0).

Lagrangian second order
-----------------------

The Lagrangian second-order scheme is forbidden if:

* Thermal variable = ``temperature`` AND scale = ``kelvin``
* Thermal variable = ``enthalpy``
"""
        self._write_file(p / "validation.rst", content)

    def _gen_reference(self, p):
        content = """Reference table
====================

.. list-table::
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Variable
     - Diffusivity
     - Thermo plane
     - is_temperature
     - Main field
   * - temperature
     - lambda
     - PT
     - 1
     - temperature
   * - enthalpy
     - lambda/cp
     - PH
     - 0
     - enthalpy
   * - internal_energy
     - lambda/cp
     - PH
     - 0
     - total_energy
   * - none
     - -
     - -
     - -
     - -
"""
        self._write_file(p / "reference.rst", content)


# =============================================================================
# TIME STEPPING GENERATOR
# =============================================================================

class TimeSteppingDocGenerator(BaseXMLDocGenerator):

    def __init__(self, xml_path):
        super().__init__(xml_path, "TimeStepping")
        self.time_schemes = self._get_typedef_values("TimeScheme")
        self.extrap_methods = self._get_typedef_values("ExtrapolationMethod")
        self.source_orders = self._get_typedef_values("SourceTermOrder")
        self._load_theta_schemes()

    def _load_theta_schemes(self):
        self.theta_schemes = {}
        for rule in self.root.findall(".//Rule[@Type='ThetaScheme']"):
            prop = rule.get("Property")
            desc_node = rule.find("Description")
            desc = desc_node.text.strip() if desc_node is not None else prop
            extrap = {e.get("Method"): e.get("Theta") for e in rule.findall("Extrapolation")}
            source = {s.get("Order"): s.get("Theta") for s in rule.findall("SourceTermOrder")}
            self.theta_schemes[prop] = {"description": desc, "extrapolation": extrap, "source": source}

    def generate(self, output_dir="source/timestepping"):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        print("Generating TimeSteppingRules.xml documentation...")
        self._gen_index(output_path)
        self._gen_overview(output_path)
        self._gen_theta(output_path)
        self._gen_limits(output_path)
        self._gen_validation(output_path)
        self._gen_reference(output_path)
        print(f"  -> RST files in {output_path}\n")

    def _gen_index(self, p):
        content = f"""TimeSteppingRules.xml -- Technical Reference
============================================

Auto-generated documentation from ``TimeSteppingRules.xml``

.. note::
   Generated on {datetime.now().strftime("%d/%m/%Y at %H:%M")}.

.. toctree::
   :maxdepth: 2
   :numbered:

   overview
   theta
   limits
   validation
   reference

* :ref:`genindex`
* :ref:`search`
"""
        self._write_file(p / "index.rst", content)

    def _gen_overview(self, p):
        content = """Overview
========

``TimeSteppingRules.xml`` centralizes the time discretization rules.

Il est lu par :

* Le **kernel C++** via ``cs_time_stepping_rules_manager``
* La **GUI Python** via ``TimeSteppingRulesManager`` (ElementTree)

Supported time schemes
----------------------------

"""
        scheme_desc = {
            "steady": "Steady state",
            "first_order": "First-order implicit Euler scheme",
            "second_order": "Second-order Crank-Nicolson scheme",
        }
        for s in self.time_schemes:
            desc = scheme_desc.get(s, s)
            content += f"* ``{s}`` : {desc}\n"
        content += "\n\nExtrapolation methods\n---------------------\n\n"
        extrap_desc = {
            "none": "No extrapolation (theta=0)",
            "linear": "Linear extrapolation (theta=0.5)",
            "second_order": "Second-order extrapolation (theta=1)",
        }
        for e in self.extrap_methods:
            desc = extrap_desc.get(e, e)
            content += f"* ``{e}`` : {desc}\n"
        self._write_file(p / "overview.rst", content)

    def _gen_theta(self, p):
        content = """Theta Schemes
=============

The theta scheme defines the time extrapolation level:

* **theta = 0** : explicite
* **theta = 0.5** : extrapolation en n+1/2
* **theta = 1** : extrapolation en n+1

Configurations per property
-----------------------------

"""
        for prop, info in self.theta_schemes.items():
            content += f"\n{prop}\n" + "^" * len(prop) + "\n\n"
            content += f"{info['description']}\n\n"
            if info["extrapolation"]:
                content += ".. list-table::\n   :widths: 40 30 30\n   :header-rows: 1\n\n"
                content += "   * - Methode\n     - Theta\n     - Description\n"
                extrap_desc = {
                    "none": "Explicit",
                    "linear": "Extrapolated at n+1/2",
                    "second_order": "Extrapolated at n+1",
                }
                for method, theta in info["extrapolation"].items():
                    desc = extrap_desc.get(method, method)
                    content += f"   * - {method}\n     - {theta}\n     - {desc}\n"
                content += "\n"
            if info["source"]:
                content += ".. list-table::\n   :widths: 40 30 30\n   :header-rows: 1\n\n"
                content += "   * - Ordre terme source\n     - Theta\n     - Description\n"
                source_desc = {
                    "explicit": "Explicit source terms",
                    "semi_implicit_1": "Semi-implicit theta=0.5",
                    "semi_implicit_2": "Semi-implicit theta=1",
                }
                for order, theta in info["source"].items():
                    desc = source_desc.get(order, order)
                    content += f"   * - {order}\n     - {theta}\n     - {desc}\n"
                content += "\n"
        self._write_file(p / "theta.rst", content)

    def _gen_limits(self, p):
        content = """Time step limits
=======================

.. list-table::
   :widths: 25 15 15 45
   :header-rows: 1

   * - Parameter
     - Default
     - Min
     - Description
   * - dt_ref
     - 0.1
     - > 0
     - Reference time step
   * - dt_min
     - 1e-5
     - >= 0
     - Minimum time step
   * - dt_max
     - 1e4
     - > 0
     - Maximum time step
   * - cfl_max
     - 1.0
     - > 0
     - Maximum Courant number
   * - courant_max
     - 1.0
     - > 0
     - Maximum Courant
   * - fourier_max
     - 10.0
     - > 0
     - Maximum Fourier number

Constraint: ``dt_min < dt_ref < dt_max``
"""
        self._write_file(p / "limits.rst", content)

    def _gen_validation(self, p):
        content = """Validation rules
====================

Parameter theta
---------------

* Must be in [0, 1]
* Do not modify manually: automatically initialized from XML
* Sentinel value: -999 (triggers auto-init)

Time step
---------

* ``dt_ref > 0``
* ``dt_min >= 0``
* ``dt_max > 0``
* ``dt_max > dt_min``
"""
        self._write_file(p / "validation.rst", content)

    def _gen_reference(self, p):
        content = """Reference table
====================

Extrapolation method -> theta
---------------------------------

.. list-table::
   :widths: 40 20 40
   :header-rows: 1

   * - Method
     - Theta
     - Properties
   * - none (iviext=0)
     - 0.0
     - viscosity, heat_capacity, scalar_diffusivity
   * - linear (iviext=1)
     - 0.5
     - viscosity, heat_capacity, scalar_diffusivity
   * - second_order (iviext=2)
     - 1.0
     - viscosity, heat_capacity, scalar_diffusivity

Source term order -> theta
----------------------------

.. list-table::
   :widths: 40 20 40
   :header-rows: 1

   * - Order
     - Theta
     - Source terms
   * - explicit (isno2t=0)
     - 0.0
     - NS, turbulence, scalaires
   * - semi_implicit_1 (isno2t=1)
     - 0.5
     - NS, turbulence, scalaires
   * - semi_implicit_2 (isno2t=2)
     - 1.0
     - NS, turbulence, scalaires
"""
        self._write_file(p / "reference.rst", content)



# =============================================================================
# BOUNDARY CONDITIONS GENERATOR
# =============================================================================

class BoundaryConditionsDocGenerator(BaseXMLDocGenerator):

    def __init__(self, xml_path):
        super().__init__(xml_path, "BoundaryConditions")
        self.valid_natures = self._get_typedef_values("BCNature")
        self.valid_velocity_choices = self._get_typedef_values("VelocityChoice")
        self.valid_velocity_directions = self._get_typedef_values("VelocityDirection")
        self.valid_turbulence_choices = self._get_typedef_values("TurbulenceChoice")
        self.valid_scalar_bc_types = self._get_typedef_values("ScalarBC")
        self._load_turb_inlet_rules()

    def _load_turb_inlet_rules(self):
        self.turb_inlet_rules = {}
        for rule in self.root.findall(".//Rule[@Type='TurbulenceInlet']"):
            model = rule.get("Model")
            if model:
                nbv = rule.find("NBoundaryValues")
                meg = rule.find("MEGFunction")
                vars_list = [v.text.strip() for v in rule.findall("Variable") if v.text]
                self.turb_inlet_rules[model] = {
                    "n_boundary_values": int(nbv.text) if nbv is not None else 1,
                    "meg_function": meg.text.strip() if meg is not None else "",
                    "variables": vars_list
                }

    def generate(self, output_dir="source/boundary_conditions"):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        print("Generating BoundaryConditionsRules.xml documentation...")
        self._gen_index(output_path)
        self._gen_overview(output_path)
        self._gen_natures(output_path)
        self._gen_turbulence_inlet(output_path)
        self._gen_scalars(output_path)
        self._gen_validation(output_path)
        self._gen_reference(output_path)
        print(f"  -> RST files in {output_path}\n")

    def _gen_index(self, p):
        content = f"""BoundaryConditionsRules.xml -- Technical Reference
====================================================

Auto-generated documentation from ``BoundaryConditionsRules.xml``

.. note::
   Generated on {datetime.now().strftime("%d/%m/%Y at %H:%M")}.

.. toctree::
   :maxdepth: 2
   :numbered:

   overview
   natures
   turbulence_inlet
   scalars
   validation
   reference

* :ref:`genindex`
* :ref:`search`
"""
        self._write_file(p / "index.rst", content)

    def _gen_overview(self, p):
        content = """Overview
========

``BoundaryConditionsRules.xml`` centralizes the boundary condition
rules of Code_Saturne. It is read by:

* Le **kernel C++** via ``cs_boundary_conditions_rules_manager``
* La **GUI Python** via ``BoundaryConditionsRulesManager`` (ElementTree)

The main contribution is the turbulence model -> inlet BC variables mapping,
which replaces 200 lines of hardcoded if/else if in
``cs_gui_boundary_conditions.cpp``.

Main sections
--------------------

* ``Definitions`` : valid types (natures, velocity choices, turbulence, scalars)
* ``ValidationRules`` : mapping turbulence -> BC inlet, contraintes
* ``Defaults`` : default values
* ``Mappings`` : mappings to C++ flags
"""
        self._write_file(p / "overview.rst", content)

    def _gen_natures(self, p):
        content = """Boundary condition types
========================

BC nature
---------

"""
        nature_desc = {
            "inlet": "Fluid inlet (imposed velocity)",
            "wall": "Wall (zero or slip velocity)",
            "outlet": "Fluid outlet (imposed pressure)",
            "symmetry": "Symmetry plane",
            "free_inlet_outlet": "Free inlet/outlet",
            "imposed_p_outlet": "Imposed pressure outlet",
            "free_surface": "Free surface",
            "groundwater": "Porous medium flow",
            "undefined": "Undefined type",
        }
        for n in self.valid_natures:
            desc = nature_desc.get(n, n)
            content += f"* ``{n}`` : {desc}\n"

        content += """
\nChoix de vitesse a l'entree
-----------------------------

"""
        vel_desc = {
            "norm": "Imposed velocity norm",
            "flow1": "Imposed mass flow rate",
            "flow2": "Imposed volume flow rate",
            "norm_formula": "Norm via MEG formula",
            "flow1_formula": "Mass flow rate via MEG formula",
            "flow2_formula": "Volume flow rate via MEG formula",
        }
        for v in self.valid_velocity_choices:
            desc = vel_desc.get(v, v)
            content += f"* ``{v}`` : {desc}\n"

        content += """
\nDirection de la vitesse
------------------------

"""
        dir_desc = {
            "normal": "Normal direction to the face",
            "coordinates": "Direction by Cartesian coordinates",
            "formula": "Direction via MEG formula",
        }
        for d in self.valid_velocity_directions:
            desc = dir_desc.get(d, d)
            content += f"* ``{d}`` : {desc}\n"

        self._write_file(p / "natures.rst", content)

    def _gen_turbulence_inlet(self, p):
        content = """Turbulence at inlet
===================

Model -> BC variables mapping
------------------------------

This table replaces 200 lines of hardcoded if/else if in
``cs_gui_boundary_conditions.cpp``. Generated from
``BoundaryConditionsRules.xml``.

.. list-table::
   :widths: 25 15 25 35
   :header-rows: 1

   * - Model
     - N values
     - MEG function
     - Variables
"""
        for model, rule in sorted(self.turb_inlet_rules.items()):
            vars_str = ", ".join(f"``{v}``" for v in rule["variables"])
            content += f"   * - {model}\n"
            content += f"     - {rule['n_boundary_values']}\n"
            content += f"     - {rule['meg_function']}\n"
            content += f"     - {vars_str}\n"

        content += """
\nMethodes de specification de turbulence
------------------------------------------

"""
        turb_desc = {
            "hydraulic_diameter": "Hydraulic diameter (smooth wall law)",
            "turbulent_intensity": "Turbulent intensity in percent",
            "formula": "User MEG formula",
        }
        for t in self.valid_turbulence_choices:
            desc = turb_desc.get(t, t)
            content += f"* ``{t}`` : {desc}\n"

        self._write_file(p / "turbulence_inlet.rst", content)

    def _gen_scalars(self, p):
        content = """Scalar boundary conditions
==========================

Scalar BC types
---------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Type
     - Description
"""
        scalar_desc = {
            "dirichlet": "Imposed value (Dirichlet)",
            "neumann": "Imposed flux (Neumann)",
            "dirichlet_formula": "Imposed value via MEG formula",
            "dirichlet_implicit": "Implicit Dirichlet (coupling)",
            "neumann_implicit": "Implicit Neumann (coupling)",
        }
        for s in self.valid_scalar_bc_types:
            desc = scalar_desc.get(s, s)
            content += f"   * - ``{s}``\n     - {desc}\n"

        self._write_file(p / "scalars.rst", content)

    def _gen_validation(self, p):
        content = """Validation rules
====================

Velocity
--------

* ``velocity_norm`` : must be strictly positive

Hydraulic diameter
------------------

* ``hydraulic_diameter`` : must be strictly positive

Turbulent intensity
-------------------

* ``turbulent_intensity`` : must be in [0, 100] (percentage)

BC nature
---------

Doit etre une des valeurs definies dans ``Definitions/TypeDef[@name='BCNature']``.
"""
        self._write_file(p / "validation.rst", content)

    def _gen_reference(self, p):
        content = """Reference table
====================

Nature -> C++ flag
------------------

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Nature
     - Flag C++
   * - inlet
     - CS_BOUNDARY_INLET
   * - wall
     - CS_BOUNDARY_WALL
   * - outlet
     - CS_BOUNDARY_OUTLET
   * - symmetry
     - CS_BOUNDARY_SYMMETRY
   * - free_inlet_outlet
     - CS_BOUNDARY_FREE_INLET_OUTLET
   * - imposed_p_outlet
     - CS_BOUNDARY_IMPOSED_P_OUTLET
   * - free_surface
     - CS_BOUNDARY_FREE_SURFACE
   * - groundwater
     - CS_BOUNDARY_GROUNDWATER

Defaults
--------

.. list-table::
   :widths: 40 60
   :header-rows: 1

   * - Parameter
     - Default value
   * - bc_nature
     - inlet
   * - velocity_choice
     - norm
   * - velocity_direction
     - normal
   * - turbulence_choice
     - hydraulic_diameter
   * - scalar_bc_type
     - dirichlet
"""
        self._write_file(p / "reference.rst", content)

# =============================================================================
# MAIN INDEX + CONF
# =============================================================================

def generate_main_index(output_dir="source"):
    p = Path(output_dir)
    p.mkdir(parents=True, exist_ok=True)
    content = f"""Code_Saturne -- Model Rules Documentation
=========================================

Auto-generated documentation from XML business rules files.

.. note::
   Generated on {datetime.now().strftime("%d/%m/%Y at %H:%M")}.

Architecture
------------

Business rules are centralized in XML files shared between
the C++ kernel and the Python GUI (Repository + Layers pattern).

.. code-block:: text

   data_model/
   |-- TurbulenceRules.xml          <- Turbulence models
   |-- ThermalRules.xml             <- Thermal model
   |-- TimeSteppingRules.xml        <- Time stepping schemes
   |-- BoundaryConditionsRules.xml  <- Boundary conditions
   |-- LagrangianRules.xml          <- Lagrangian module
   `-- ElectricalRules.xml          <- Electrical module

Documented modules
------------------

.. toctree::
   :maxdepth: 1

   turbulence/index
   thermal/index
   timestepping/index
   boundary_conditions/index
   lagrangian/index
   electrical/index

* :ref:`genindex`
* :ref:`search`
"""
    with open(p / "index.rst", 'w', encoding='utf-8') as f:
        f.write(content)
    print("  -> index.rst (main)")


def generate_main_conf(output_dir="source"):
    p = Path(output_dir)
    p.mkdir(parents=True, exist_ok=True)
    content = """# Configuration Sphinx -- Code_Saturne Model Rules
import os
import sys
from pathlib import Path

project = "Code_Saturne Model Rules"
release = "1.0.0"
language = "en"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
]

napoleon_google_docstring = True
napoleon_numpy_docstring = False
autodoc_member_order = "bysource"

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
}

exclude_patterns = []
"""
    with open(p / "conf.py", 'w', encoding='utf-8') as f:
        f.write(content)
    print("  -> conf.py")



# =============================================================================
# LAGRANGIAN GENERATOR
# =============================================================================
class LagrangianDocGenerator(BaseXMLDocGenerator):

    def __init__(self, xml_path):
        super().__init__(xml_path, "Lagrangian")

    def generate(self, output_dir="source/lagrangian"):
        p = Path(output_dir)
        p.mkdir(parents=True, exist_ok=True)
        print(f"Generating Lagrangian documentation in {output_dir}/")
        self._gen_index(p)
        self._gen_overview(p)
        self._gen_coupling_modes(p)
        self._gen_physical_models(p)
        self._gen_defaults(p)
        self._gen_reference(p)

    def _gen_index(self, p):
        content = """Lagrangian Module
=================

.. toctree::
   :maxdepth: 2

   overview
   coupling_modes
   physical_models
   defaults
   reference
"""
        self._write_file(p / "index.rst", content)

    def _gen_overview(self, p):
        content = """Overview
========

The Lagrangian module in Code_Saturne allows simulation of dispersed
phase flows (particles, droplets, bubbles) coupled with the continuous phase.

Rules are centralized in ``LagrangianRules.xml`` and parsed at runtime
by ``cs_lagr_rules_manager`` (C++ singleton).

Coupling Modes
--------------

Four coupling modes are available:

- **off**: Lagrangian module disabled
- **one_way**: particles follow the flow, no feedback
- **two_way**: mutual influence between phases
- **frozen**: continuous phase is frozen (restart only)

Physical Models
---------------

Particle physics can be enriched with:

- **thermal**: heat transfer
- **coal**: pulverised coal combustion
- **ctwr**: cooling tower
"""
        self._write_file(p / "overview.rst", content)

    def _gen_coupling_modes(self, p):
        defs = self.root.find("Definitions")
        typedef = defs.find("TypeDef[@name='LagrangianCouplingMode']") if defs is not None else None

        lines = ["Coupling Modes", "==============", ""]
        lines.append(".. list-table:: Lagrangian coupling modes")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 15 25 10 50")
        lines.append("")
        lines.append("   * - Name")
        lines.append("     - C++ Enum")
        lines.append("     - Value")
        lines.append("     - Description")

        if typedef is not None:
            for v in typedef.findall("Value"):
                name  = v.get("name", "")
                enum  = v.get("Enum", "")
                iv    = v.get("Int", "")
                desc  = v.get("Description", "")
                lines.append(f"   * - ``{name}``")
                lines.append(f"     - ``{enum}``")
                lines.append(f"     - {iv}")
                lines.append(f"     - {desc}")

        self._write_file(p / "coupling_modes.rst", "\n".join(lines) + "\n")

    def _gen_physical_models(self, p):
        defs = self.root.find("Definitions")
        typedef = defs.find("TypeDef[@name='ParticlesModel']") if defs is not None else None

        lines = ["Physical Models", "===============", ""]
        lines.append(".. list-table:: Particle physical models")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 15 25 10 50")
        lines.append("")
        lines.append("   * - Name")
        lines.append("     - C++ Enum")
        lines.append("     - Value")
        lines.append("     - Description")

        if typedef is not None:
            for v in typedef.findall("Value"):
                name  = v.get("name", "")
                enum  = v.get("Enum", "")
                iv    = v.get("Int", "")
                desc  = v.get("Description", "")
                lines.append(f"   * - ``{name}``")
                lines.append(f"     - ``{enum}``")
                lines.append(f"     - {iv}")
                lines.append(f"     - {desc}")

        self._write_file(p / "physical_models.rst", "\n".join(lines) + "\n")

    def _gen_defaults(self, p):
        dflts = self.root.find("Defaults")
        lines = ["Default Values", "==============", ""]
        lines.append(".. list-table:: Default parameter values")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 40 60")
        lines.append("")
        lines.append("   * - Parameter")
        lines.append("     - Default value")

        if dflts is not None:
            for d in dflts.findall("Default"):
                key = d.get("Key", "")
                val = d.get("Value", "")
                lines.append(f"   * - ``{key}``")
                lines.append(f"     - ``{val}``")

        self._write_file(p / "defaults.rst", "\n".join(lines) + "\n")

    def _gen_reference(self, p):
        content = """Reference
=========

Source files
------------

- ``data_model/LagrangianRules.xml`` — XML rules file
- ``src/model/cs_lagr_rules_manager.h`` — C++ manager header
- ``src/model/cs_lagr_rules_manager.cpp`` — C++ manager implementation
- ``python/code_saturne/model/LagrangianModel.py`` — Python GUI model

C++ API
-------

.. code-block:: cpp

   // Get singleton instance
   cs_lagr_rules_manager *mgr = cs_get_lagr_rules_manager();

   // Get coupling mode int value
   int mode = mgr->get_coupling_mode("two_way");  // returns 2

   // Get physical model int value
   int phys = mgr->get_phys_model("coal");  // returns 2

   // Get default values
   double t_thresh = mgr->get_default_double("threshold_temperature");  // 600.0
   int scheme     = mgr->get_default_int("scheme_order");               // 1
"""
        self._write_file(p / "reference.rst", content)


# =============================================================================
# ELECTRICAL GENERATOR
# =============================================================================
class ElectricalDocGenerator(BaseXMLDocGenerator):

    def __init__(self, xml_path):
        super().__init__(xml_path, "Electrical")

    def generate(self, output_dir="source/electrical"):
        p = Path(output_dir)
        p.mkdir(parents=True, exist_ok=True)
        print(f"Generating Electrical documentation in {output_dir}/")
        self._gen_index(p)
        self._gen_overview(p)
        self._gen_models(p)
        self._gen_properties(p)
        self._gen_validation(p)
        self._gen_reference(p)

    def _gen_index(self, p):
        content = """Electrical Module
=================

.. toctree::
   :maxdepth: 2

   overview
   models
   properties
   validation
   reference
"""
        self._write_file(p / "index.rst", content)

    def _gen_overview(self, p):
        content = """Overview
========

The Electrical module in Code_Saturne handles two types of models:

- **Joule effect**: resistive heating (AC/DC, three-phase, with or without transformer)
- **Electric arc**: high-temperature plasma arc simulation

Rules are centralized in ``ElectricalRules.xml`` and parsed at runtime
by ``cs_elec_rules_manager`` (C++ singleton).

Key parameters
--------------

- ``ieljou`` : Joule effect model integer code (-1, 1, 2, 3, 4)
- ``ielarc`` : Electric arc model integer code (-1, 2)
- ``ielcor`` : Scaling correction mode (0=off, 1=on)
"""
        self._write_file(p / "overview.rst", content)

    def _gen_models(self, p):
        defs = self.root.find("Definitions")

        lines = ["Models", "======", ""]

        # Joule models
        typedef = defs.find("TypeDef[@name='JouleModel']") if defs is not None else None
        lines.append("Joule Effect Models")
        lines.append("-------------------")
        lines.append("")
        lines.append(".. list-table::")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 30 10 60")
        lines.append("")
        lines.append("   * - Name")
        lines.append("     - ieljou")
        lines.append("     - Description")
        if typedef is not None:
            for v in typedef.findall("Value"):
                name = v.get("name", "")
                iv   = v.get("Int", "")
                desc = v.get("Description", "")
                lines.append(f"   * - ``{name}``")
                lines.append(f"     - {iv}")
                lines.append(f"     - {desc}")

        lines.append("")

        # Arc models
        typedef = defs.find("TypeDef[@name='ElectricArcModel']") if defs is not None else None
        lines.append("Electric Arc Models")
        lines.append("-------------------")
        lines.append("")
        lines.append(".. list-table::")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 30 10 60")
        lines.append("")
        lines.append("   * - Name")
        lines.append("     - ielarc")
        lines.append("     - Description")
        if typedef is not None:
            for v in typedef.findall("Value"):
                name = v.get("name", "")
                iv   = v.get("Int", "")
                desc = v.get("Description", "")
                lines.append(f"   * - ``{name}``")
                lines.append(f"     - {iv}")
                lines.append(f"     - {desc}")

        self._write_file(p / "models.rst", "\n".join(lines) + "\n")

    def _gen_properties(self, p):
        props = self.root.find("PhysicalProperties")
        lines = ["Physical Properties", "===================", ""]
        lines.append(".. list-table::")
        lines.append("   :header-rows: 1")
        lines.append("   :widths: 25 35 15 25")
        lines.append("")
        lines.append("   * - Name")
        lines.append("     - Label")
        lines.append("     - Unit")
        lines.append("     - Model")
        if props is not None:
            for p_node in props.findall("Property"):
                name  = p_node.get("name", "")
                label = p_node.get("Label", "")
                unit  = p_node.get("Unit", "")
                model = p_node.get("Model", "all")
                lines.append(f"   * - ``{name}``")
                lines.append(f"     - {label}")
                lines.append(f"     - {unit}")
                lines.append(f"     - {model}")
        self._write_file(p / "properties.rst", "\n".join(lines) + "\n")

    def _gen_validation(self, p):
        content = """Validation Rules
================

Source: ``cs_elec_model.cpp`` — ``_cs_electrical_model_verify()``

ieljou allowed values
---------------------

Only the following values are permitted for ``ieljou``:

- ``-1`` : Joule effect disabled
- ``1``  : AC/DC
- ``2``  : Three-phase
- ``3``  : AC/DC with transformer
- ``4``  : Three-phase with transformer

ielarc allowed values
---------------------

Only the following values are permitted for ``ielarc``:

- ``-1`` : Electric arc disabled
- ``2``  : Electric arc enabled

ielcor allowed values
---------------------

- ``0`` : Scaling disabled
- ``1`` : Scaling enabled

Constraints when scaling is active (ielcor=1)
---------------------------------------------

With electric arc (ielarc > 0):

- ``couimp`` must be strictly positive
- ``pot_diff`` must be strictly positive

With Joule effect (ieljou > 0):

- ``puisim`` must be strictly positive
- ``coejou`` must be strictly positive
- ``pot_diff`` must be strictly positive
"""
        self._write_file(p / "validation.rst", content)

    def _gen_reference(self, p):
        content = """Reference
=========

Source files
------------

- ``data_model/ElectricalRules.xml`` — XML rules file
- ``src/model/cs_elec_rules_manager.h`` — C++ manager header
- ``src/model/cs_elec_rules_manager.cpp`` — C++ manager implementation
- ``src/elec/cs_elec_model.h/.cpp`` — Electrical model kernel
- ``python/code_saturne/model/ElectricalModel.py`` — Python GUI model

C++ API
-------

.. code-block:: cpp

   // Get singleton instance
   cs_elec_rules_manager *mgr = cs_get_elec_rules_manager();

   // Get ieljou value from model name
   int ieljou = mgr->get_ieljou("AC/DC");  // returns 1

   // Get ielarc value
   int ielarc = mgr->get_ielarc("on");  // returns 2

   // Validate values
   bool ok = mgr->is_valid_ieljou(3);  // returns true
   bool ok2 = mgr->is_valid_ielarc(5); // returns false

   // Get default value
   double coejou = mgr->get_default_double("coejou");  // 1.0

   // Get properties for arc model
   auto props = mgr->get_properties("arc");
"""
        self._write_file(p / "reference.rst", content)

# =============================================================================
# MAIN
# =============================================================================

def main():
    xml_turb = os.path.join(DATA_DIR, "TurbulenceRules.xml")
    xml_therm = os.path.join(DATA_DIR, "ThermalRules.xml")
    xml_time = os.path.join(DATA_DIR, "TimeSteppingRules.xml")
    xml_bc = os.path.join(DATA_DIR, "BoundaryConditionsRules.xml")
    xml_lagr = os.path.join(DATA_DIR, "LagrangianRules.xml")
    xml_elec = os.path.join(DATA_DIR, "ElectricalRules.xml")

    for xml_path, name in [(xml_turb, "TurbulenceRules.xml"),
                            (xml_therm, "ThermalRules.xml"),
                            (xml_time, "TimeSteppingRules.xml"),
                            (xml_bc, "BoundaryConditionsRules.xml"),
                            (xml_lagr, "LagrangianRules.xml"),
                            (xml_elec, "ElectricalRules.xml")]:
        if not os.path.exists(xml_path):
            print(f"Error: {name} not found at {xml_path}")
            sys.exit(1)

    print("=" * 60)
    print("  CODE_SATURNE -- MODEL RULES DOCUMENTATION GENERATOR")
    print("=" * 60)
    print()

    generate_main_index("source")
    generate_main_conf("source")
    print()

    TurbulenceDocGenerator(xml_turb).generate("source/turbulence")
    ThermalDocGenerator(xml_therm).generate("source/thermal")
    TimeSteppingDocGenerator(xml_time).generate("source/timestepping")
    BoundaryConditionsDocGenerator(xml_bc).generate("source/boundary_conditions")
    LagrangianDocGenerator(xml_lagr).generate("source/lagrangian")
    ElectricalDocGenerator(xml_elec).generate("source/electrical")

    print("=" * 60)
    print("RST documentation generated.")
    print()
    print("To compile: make html")
    print("Open: build/html/index.html")
    print("=" * 60)


if __name__ == "__main__":
    main()
