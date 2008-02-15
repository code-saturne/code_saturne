#############################
#
# Makefile pour Documentation
#
#############################


# MACROS ENVIRONNEMENT
# ====================

SHELL           =/bin/sh
TERM            =xterm
HOME            =./.
FIG             =$(HOME)/graphics
TEX             =$(HOME)/latex


# DEFINITION DES OBJETS A CREER
#==============================

DOC             =user


# DEFINITIONS DES CIBLES
#=======================

all doc: pdf install


fig:
	@cd $(FIG) ;\
	 $(MAKE)


pdf: fig
	@cd $(TEX) ;\
	 $(MAKE) NOTE=$(DOC) pdf


clean:
	@cd $(TEX) ;\
	 $(MAKE) clean
	@cd $(FIG) ;\
	 $(MAKE) clean


distclean: clean
	@\rm -f $(FIG)/*~
	@\rm -f $(TEX)/*~
	@\rm -f $(DOC).pdf


install: 
	@\mv $(TEX)/$(DOC).pdf .


help:
	@echo
	@echo "Usage : make"
	@echo "   or : make <target>"
	@echo
	@echo " with :"
	@echo "    <target> : doc       : Compile the Code_Saturne user guide"
	@echo "                           in PDF format"
	@echo "               all       : Same as doc"
	@echo "               clean     : Remove temporary files"
	@echo "               distclean : Clean directories completely"
	@echo "               help      : This message"
	@echo
