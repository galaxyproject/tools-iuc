EMBOSS_ACD_LOCATION=~/.linuxbrew/Cellar/emboss/6.6.0/share/EMBOSS/acd
ACDS:= $(wildcard ${EMBOSS_ACD_LOCATION}/*.acd)

.PHONY: all clean

all: $(ACDS:=.xml)
	@echo "Done"

output:
	mkdir output

%.xml: acd output
	-cat $(basename $@) | ./acd | xmllint --pretty 1 - > output/$(notdir $@)

lex.yy.c:
	flex acd.l

y.tab.c:
	yacc -Wall -d acd.y

acd: lex.yy.c y.tab.c
	g++ -Wall -Wno-unused-function -g lex.yy.c y.tab.c -o acd -lfl -D_GNU_SOURCE
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=30426


# In case we need it?
groups.json: output
	python convert-groups.py ${EMBOSS_ACD_LOCATION}/groups.standard | json_pp > output/groups.json

clean:
	@rm -fr output lex.yy.c acd y.tab.c y.tab.h 
