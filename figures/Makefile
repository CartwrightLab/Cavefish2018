RSCRIPT=Rscript --vanilla

default: all

all: fig-molino.pdf molino.mp4

.PHONY: all default

.DELETE_ON_ERROR:

.SECONDARY:

%.pdf: %.R
	$(RSCRIPT) $<

fig-molino.pdf molino.mp4: fig-molino.R
