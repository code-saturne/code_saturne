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

DOC             =theory


# DEFINITIONS DES CIBLES
#=======================

all doc: pdf install


fig:
	@cd $(FIG) ;\
	 $(MAKE)


pdf:
	@cd $(TEX) ;\
	 $(MAKE) NOTE=$(DOC) pdf


clean:
	@cd $(TEX) ;\
	 $(MAKE) clean


distclean:
	@make clean
	@cd $(FIG) ;\
	 $(MAKE) clean
	@\rm -f $(FIG)/*~
	@\rm -f $(TEX)/*~
	@\rm -f $(DOC).pdf


install: 
	@\mv $(TEX)/$(DOC).pdf .


help info aide:
	@echo
	@echo "Usage : make"
	@echo "   or : make <target>              (in  $(HOME)/)"
	@echo
	@echo " avec :"
	@echo "    <cible> : doc       : Compilation de la documentation"
	@echo "                          utilisateur de Code_Saturne"
	@echo "                          au format PDF"
	@echo "              all       : idem doc"
	@echo "              clean     : Suppression des fichiers temporaires"
	@echo "              distclean : Nettoyage complet des répertoires"
	@echo "              help      : ce message"
	@echo "              info      : idem help"
	@echo "              aide      : idem help"
	@echo
