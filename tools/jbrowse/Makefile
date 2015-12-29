TEMPDIR=/tmp/jbrowse-guarunteed-non-conflicting-temp-dir

all: test-data/bam/trackList.json test-data/blastxml/trackList.json test-data/bw/trackList.json test-data/gencode/trackList.json test-data/gff3/trackList.json test-data/vcf/trackList.json 

test-data/vcf/trackList.json:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/vcf/test.xml
	cp ${TEMPDIR}/data/trackList.json test-data/vcf/
	rm -rf ${TEMPDIR}

test-data/blastxml/trackList.json:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/blastxml/test.xml
	cp ${TEMPDIR}/data/trackList.json test-data/blastxml/
	rm -rf ${TEMPDIR}

test-data/bw/trackList.json:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/bw/test.xml
	cp ${TEMPDIR}/data/trackList.json test-data/bw/
	rm -rf ${TEMPDIR}

test-data/bam/trackList.json:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/bam/test.xml
	cp ${TEMPDIR}/data/trackList.json test-data/bam/
	rm -rf ${TEMPDIR}

test-data/gff3/trackList.json:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/gff3/test.xml
	cp ${TEMPDIR}/data/trackList.json test-data/gff3/
	rm -rf ${TEMPDIR}

test-data/gencode/trackList.json:
	python jbrowse.py --jbrowse JBrowse-1.12.0 --outdir ${TEMPDIR} --standalone test-data/gencode/test.xml
	cp ${TEMPDIR}/data/trackList.json test-data/gencode/
	rm -rf ${TEMPDIR}
