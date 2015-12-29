TEMPDIR=/tmp/jbrowse-guarunteed-non-conflicting-temp-dir
TARGETS:=test-data/bam/trackList.json.comparable test-data/blastxml/trackList.json.comparable test-data/bw/trackList.json.comparable test-data/gencode/trackList.json.comparable test-data/gff3/trackList.json.comparable test-data/vcf/trackList.json.comparable

all: $(TARGETS)
clean:
	rm $(TARGETS)

test-data/vcf/trackList.json.comparable:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/vcf/test.xml
	cat ${TEMPDIR}/data/trackList.json  | sed 's/,\s*$$//g' | sort > test-data/vcf/trackList.json.comparable
	rm -rf ${TEMPDIR}

test-data/blastxml/trackList.json.comparable:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/blastxml/test.xml
	cat ${TEMPDIR}/data/trackList.json  | sed 's/,\s*$$//g' | sort > test-data/blastxml/trackList.json.comparable
	rm -rf ${TEMPDIR}

test-data/bw/trackList.json.comparable:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/bw/test.xml
	cat ${TEMPDIR}/data/trackList.json  | sed 's/,\s*$$//g' | sort > test-data/bw/trackList.json.comparable
	rm -rf ${TEMPDIR}

test-data/bam/trackList.json.comparable:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/bam/test.xml
	cat ${TEMPDIR}/data/trackList.json  | sed 's/,\s*$$//g' | sort > test-data/bam/trackList.json.comparable
	rm -rf ${TEMPDIR}

test-data/gff3/trackList.json.comparable:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/gff3/test.xml
	cat ${TEMPDIR}/data/trackList.json  | sed 's/,\s*$$//g' | sort > test-data/gff3/trackList.json.comparable
	rm -rf ${TEMPDIR}

test-data/gencode/trackList.json.comparable:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/gencode/test.xml
	cat ${TEMPDIR}/data/trackList.json  | sed 's/,\s*$$//g' | sort > test-data/gencode/trackList.json.comparable
	rm -rf ${TEMPDIR}
