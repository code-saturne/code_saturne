Contributing to code_saturne
============================

This page provides a high level overview of steps and recommendations for
contributions to code_saturne.

Licence agreement
-----------------

EDF distributes code_saturne under a dual licence model:
- GPL licence (v2 or above)
  - With a more permissive LGPL licence for the PLE library.
- Specific licence agreements allow compatibility with non-free software based
  on code_saturne, such as the neptune_cfd multiphase solver, which is owned
  by a consortium in which EDF is not the sole owner.

To ensure contributions are compatible with this dual-licence model,
code provided by authors not working for EDF or in the context
of a collaborative project already including a similar agreement
should sign a *Contributor Licencing Agreement*, for which a
[model is provided here](docs/doxygen/developer_guide/EDF_Open_Source_CLA.pdf).
This ensures that both EDF and the contributing party each have full rights to
use and distribute the contributed code under the licence of their choice.

If no similar licence agreement is possible, contributions provided
only under a GPL-compatible open source licence could be used by code_saturne as
external libraries, provided the dependency remains optional. This is an
option for libraries providing significant features, not minor contributions.

Minor changes such as simple bug fixes or typo corrections which do not imply
copyright aspects can be provided direclty.

Coding practice
---------------

Following the coding recommendations provided in the
[developer guide](docs/doxygen/developer_guide) is essential, as ensuring
consistency with the existing code base will reduce the workload for code
integration, and significantly increase the chances that code will actually
be integrated.

Contacting the code_saturne Development Team
---------------------------------------------

The development team may be contacted through several means, including
the [user's forum](https://www.code-saturne.org/forum/), GitHub issues or merge
requests, and the saturne-support@edf.fr contact address. The forum is
the preferred means of initial contact, though detailed discussions may be taken
offline.

We strongly urge contributors to discuss planned developments at an early
stage with the development team, to allow providing suggestions. As parts of
the code_saturne code base evolve rapidly, regular discussions with contributors
allow sharing some insights and planned changes, and ensuring contributions
may be integrated smoothly in the code.
