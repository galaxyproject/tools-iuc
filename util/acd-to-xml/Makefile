EMBOSS_ACD_LOCATION=~/.linuxbrew/Cellar/emboss/6.6.0/share/EMBOSS/acd/
EMBOSS_DOC_LOCATION=${EMBOSS_ACD_LOCATION}../doc/programs/html/
ACDS:= $(wildcard ${EMBOSS_ACD_LOCATION}*.acd)
OUTPUT_ACDS:=$(addprefix output/, $(addsuffix .xml, $(notdir $(ACDS))))
OUTPUT_GXML:=$(addprefix output/, $(addsuffix .gxml, $(notdir $(ACDS))))

.PHONY: all clean full_clean

all: $(OUTPUT_ACDS) $(OUTPUT_GXML)
	@echo "Done"

output:
	mkdir output

output/%.xml: acd output
	-cat $(EMBOSS_ACD_LOCATION)$(notdir $(basename $@)).acd | ./acd | xmllint --pretty 1 - > $@

output/%.gxml: output/%.xml
	python emboss.py $(EMBOSS_DOC_LOCATION)$(notdir $(basename $@)).html $(basename $@).xml | xmllint --pretty 1 - > $@

lex.yy.c:
	flex acd.l

y.tab.c:
	yacc -Wall -d acd.y

acd: lex.yy.c y.tab.c
	g++ -Wall -Wno-unused-function -g lex.yy.c y.tab.c -o acd -lfl -D_GNU_SOURCE
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=30426


# In case we need it?
groups.json: output
	python convert-groups.py ${EMBOSS_ACD_LOCATION}groups.standard | json_pp > output/groups.json

clean:
	@rm -fr lex.yy.c acd y.tab.c y.tab.h

full_clean: clean
	@rm -fr output
