RSCRIPT=Rscript --vanilla
SIM=../script/k-allele-sc2m-sim.R
ARGS=k=12 h=0.5 u=1e-6 Q=0.1 q0=0.1 G1=170000 G2=20000 N=30522

M=0e-6 2e-6 4e-6 6e-6 8e-6 10e-6
S=1e-5 1e-4 1e-3 1e-2 1e-1 1

RUNS=$(foreach m,$(M),$(foreach s,$(S),m=$(m)_s=$(s)))
SIMCSVS=$(addsuffix .csv,$(RUNS))

REPS:=$(shell seq -w 0 99)

.DELETE_ON_ERROR:

.SECONDARY:

default: all

all: results.csv

.PHONY: all default

header.csv:
	$(RSCRIPT) $(SIM) > $@

results.csv: header.csv $(SIMCSVS)
	cat $^ > $@

%.run.csv:
	test -d $(dir $@) || mkdir $(dir $@)
	$(RSCRIPT) $(SIM) $(ARGS) $(subst _, ,$(patsubst %/,%,$(dir $@))) > $@

$(SIMCSVS) : %.csv : $(foreach r,$(REPS),%/$(r).run.csv)
	cat $^ > $@
