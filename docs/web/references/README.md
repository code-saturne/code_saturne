
Instructions for updating references on web site:

- Update references.bib
- Generate web page using "make html"
  - This requires having sphinx and the sphinxcontrib.bibtex extension
  - Installable using pip install, though this seems to be incomplete in Scibian 9.
    - Using a separate machine might be simpler...
- Using a web browser, copy-paste the results.

A direct edition is also possible, and other conversion tools may be useful,
but the sphinx-based method is very efficient and generates clean output.

Running this (with the appropriate tools) on the web-site itself could
be interesting, as it would allow using the actual (sphinx-generated)
page, with better looking output (and requiring the appropriate
sphinx install only at one place).