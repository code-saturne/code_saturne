"""
Generate Sphinx toctree structure.
Must be run after generate_rst.py:
    python3 <srcdir>/doc/sphinx/generate_toctree.py <srcdir> <builddir>
"""
import os
import sys

if len(sys.argv) < 3:
    raise SystemExit(
        "Usage: python3 generate_toctree.py <srcdir> <builddir>\n"
        "Example: python3 /path/to/src/doc/sphinx/generate_toctree.py "
        "/path/to/src /path/to/build"
    )

srcdir   = sys.argv[1]
builddir = sys.argv[2]

sphinx_srcdir   = os.path.join(srcdir,   "doc", "sphinx")
sphinx_builddir = os.path.join(builddir, "doc", "sphinx")
rst_dir = os.path.join(sphinx_builddir, "api")
os.makedirs(sphinx_builddir, exist_ok=True)

structure = [
    ("General Information",        "general_information", ["cs_ug_quick_start","cs_ug_case_structure","cs_ug_run_computation","cs_ug_user_settings"]),
    ("General Principles",         "general_principles",  ["cs_ug_base_setup","cs_ug_udf","cs_ug_file_formats"]),
    ("Setting Up a Computation",   "setting_up",          ["cs_ug_mesh_prepare","cs_ug_mesh_select_c","advanced_ale","advanced_atmospheric","advanced_compressible","advanced_conjugate_heat_transfer","advanced_coupling","advanced_electric_arcs","advanced_particle_tracking","advanced_radiative_thermal","advanced_solidification","advanced_turbomachinery","cs_ug_cdo_gwf","low_level_boundary_condition_definitions","cs_ug_output"]),
    ("Running Computations",       "running_computations",["cs_ug_parallel","cs_ug_troubleshooting","cs_ug_studymanager","cs_ug_smgr_slurm"]),
    ("BPG",                        "bpg",                 []),
    ("Installation",               "installation",        []),
    ("Developer Guide",            "developer_guide",     ["cs_dg_developer_guide","cs_dg_programming_languages","cs_dg_coding_style","cs_dg_common_constructs","cs_dg_node_level_parallelism","cs_dg_build_system","cs_dg_writing_theory","cs_dg_further_reading","cs_dg_debugging","cs_dg_profiling"]),
    ("User Examples",              "user_examples",       ["cs_user_examples"]),
    ("Source Code Documentation",  "source_code",         []),
    ("Bug Report",                 "bug_report",          []),
]

section_files = []
for title, safe, pages in structure:
    section_rst = os.path.join(sphinx_builddir, f"section_{safe}.rst")
    with open(section_rst, "w") as f:
        f.write(f"{title}\n{'=' * len(title)}\n\n")
        if safe == "bpg":
            f.write("Best Practice Guidelines for Code_Saturne simulations.\n\n")
            f.write(".. note::\n\n   This section is under construction.\n\n")
        elif safe == "installation":
            f.write("Installation instructions for Code_Saturne.\n\n")
            f.write(".. note::\n\n   See the ``INSTALL.md`` file in the source repository.\n\n")
        elif safe == "bug_report":
            f.write("How to report a bug:\n\n")
            f.write("- \`EDF GitLab <https://code.visualworks.fr/code_saturne/code_saturne>\`_\n")
            f.write("- \`Public GitLab <https://github.com/code-saturne/code_saturne>\`_\n\n")
        elif safe == "source_code":
            f.write("Auto-generated API documentation from source code.\n\n")
        if pages:
            f.write(".. toctree::\n   :maxdepth: 2\n\n")
            for page in pages:
                rst_path = os.path.join(rst_dir, page + ".rst")
                if os.path.exists(rst_path):
                    f.write(f"   api/{page}\n")
                else:
                    print(f"  Missing: {page}")
    section_files.append(f"section_{safe}")
    print(f"Section created: {title}")

with open(os.path.join(sphinx_builddir, "index.rst"), "w") as f:
    f.write("Code_Saturne Documentation\n==========================\n\n")
    f.write("Code_Saturne is a general purpose CFD tool developed by EDF R&D.\n\n")
    f.write(".. toctree::\n   :maxdepth: 2\n   :hidden:\n\n")
    for s in section_files:
        f.write(f"   {s}\n")
    f.write("   section_model\n   section_api\n")

already_in_sections = [p for _, _, pages in structure for p in pages]
with open(os.path.join(sphinx_builddir, "section_api.rst"), "w") as f:
    f.write("API Reference\n=============\n\n")
    f.write("Full C/C++ API auto-generated from source code.\n\n")
    f.write(".. toctree::\n   :maxdepth: 1\n\n")
    if os.path.exists(rst_dir):
        for rst_file in sorted(os.listdir(rst_dir)):
            if rst_file == "index.rst":
                continue
            name = rst_file.replace(".rst", "")
            if name not in already_in_sections:
                f.write(f"   api/{name}\n")

print("\nMain index created!")
