SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

clean:
	@$(SPHINXBUILD) -M clean "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS)

.PHONY: html clean
