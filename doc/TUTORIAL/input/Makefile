FIGFILE = $(wildcard *.fig)
FIGLIST = ${patsubst %.fig, %.pdf, $(FIGFILE) }

EPSFILE = $(wildcard *.eps)
EPSLIST = ${patsubst %.eps, %.pdf, $(EPSFILE) }

PSFILE = $(wildcard *.ps)
PSLIST = ${patsubst %.ps, %.pdf, $(PSFILE) }

all: $(FIGLIST) $(EPSLIST) $(PSLIST)

%.pdf: %.fig
	fig2dev -L pdf $^ > $@

%.pdf: %.eps
	epstopdf $^

%.pdf: %.ps
	pstopdf $^

clean:
	\rm -f $(FIGLIST) $(EPSLIST) $(PSLIST) *.bak
