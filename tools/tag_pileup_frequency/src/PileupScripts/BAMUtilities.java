package PileupScripts;

import java.io.File;
import java.io.IOException;

import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;

public class BAMUtilities {
        public static double calculateStandardizationRatio(File BAM) throws IOException {
                SamReader inputSam = SamReaderFactory.makeDefault().open(BAM);
                AbstractBAMFileIndex bai = (AbstractBAMFileIndex) inputSam.indexing().getIndex();
                double counter = 0;
                double totalAligned = 0;
                double totalGenome = 0;
                
                for (int x = 0; x < bai.getNumberOfReferences(); x++) {
                        SAMSequenceRecord seq = inputSam.getFileHeader().getSequence(x);
                        totalAligned += inputSam.indexing().getIndex().getMetaData(x).getAlignedRecordCount();
                        totalGenome += seq.getSequenceLength();
                }
                CloseableIterator<SAMRecord> iter = inputSam.iterator();
                while (iter.hasNext()) {
                        SAMRecord sr = iter.next();
                        if(sr.getReadPairedFlag()) {
                                if(sr.getSecondOfPairFlag()) { counter++; } //count read 2 to remove from aligned reads
                        }
                }
                inputSam.close();
                bai.close();
                iter.close();
                                
                totalAligned -= counter;
                if(totalAligned > 0) { return (totalAligned / totalGenome); }
                else { return 1; }
        }
        
}
