DOC := mslt-eqns
PDF := $(DOC).pdf

TEX_CC := xelatex

TMP_EXTS := aux log
TMP := $(foreach EXT,$(TMP_EXTS),$(wildcard *.${EXT}))


all: $(PDF)

%.pdf: %.tex
	@$(TEX_CC) $*
	@$(TEX_CC) $*

clean:
	@rm -f $(TMP)

distclean: clean
	@rm -f $(PDF)

.PHONY: all clean distclean
