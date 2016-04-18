import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Date;

public class BAMtoscIDX {
        private static File BAM = null;
        private static File OUTPUT = null;
        
        private static int STRAND = 0;
        private static String READ = "READ1";
                
        private static int PAIR = 1;
        private static String MIN = "NaN";
        private static int MIN_INSERT = -9999;
        private static String MAX = "NaN";
        private static int MAX_INSERT = -9999;
        
        private static SamReader inputSam = null;
        private static PrintStream OUT = null;
        
        private static ArrayList<Integer> BP;
        private static ArrayList<Integer> F_OCC;
        private static ArrayList<Integer> R_OCC;
                
        private static int CHROMSTOP = -999;
        
        public static void main(String[] args) {
                loadConfig(args);
                
                System.out.println("\n" + getTimeStamp());
                
                //Open Output File
                try { OUT = new PrintStream(OUTPUT); }
                catch (FileNotFoundException e) { e.printStackTrace(); }
                
                //Check to Make Sure BAI-index file exists
                File f = new File(BAM.getAbsolutePath() + ".bai");
                if(f.exists() && !f.isDirectory()) {
                        //Print Header
                        OUT.println("#" + getTimeStamp() + ";" + BAM.getName() + ";" + READ);
                        OUT.println("chrom\tindex\tforward\treverse\tvalue");
                        
                        processREADS(); //Begin processing reads in BAM file
                        
                        OUT.close();
                } else { OUT.println("BAI Index File does not exist for: " + BAM.getName() + "\n"); }
                System.out.println(getTimeStamp());
        }
                
        public static void addTag(SAMRecord sr) {
                //Get the start of the record 
                int recordStart = sr.getUnclippedStart();//.getAlignmentStart();
                //Accounts for reverse tag reporting 3' end of tag and converting BED to IDX/GFF format
                if(sr.getReadNegativeStrandFlag()) { recordStart = sr.getUnclippedEnd(); }//.getAlignmentEnd(); }                                        
                                
                //Make sure we only add tags that have valid starts
                if(recordStart > 0 && recordStart <= CHROMSTOP) {
                        if(BP.contains(new Integer(recordStart))) {
                                int index = BP.indexOf(new Integer(recordStart));
                                if(sr.getReadNegativeStrandFlag()) {
                                        R_OCC.set(index, new Integer(R_OCC.get(index).intValue() + 1));
                                } else {
                                        F_OCC.set(index, new Integer(F_OCC.get(index).intValue() + 1));
                                }
                        } else {
                                //Sometimes the start coordinate will be out of order due to (-) strand correction
                                //Need to efficiently identify where to place it relative to the other bps
                                int index = BP.size() - 1;
                                if(index >= 0) {
                                        while(index >= 0 && recordStart < BP.get(index).intValue()) {
                                                index--;
                                        }
                                }
                                if(index < BP.size() - 1) {
                                        BP.add(index + 1, new Integer(recordStart));
                                        if(sr.getReadNegativeStrandFlag()) {
                                                R_OCC.add(index + 1, new Integer(1));
                                                F_OCC.add(index + 1, new Integer(0));
                                        } else {
                                                F_OCC.add(index + 1, new Integer(1));
                                                R_OCC.add(index + 1, new Integer(0));
                                        }
                                } else {
                                        BP.add(new Integer(recordStart));
                                        if(sr.getReadNegativeStrandFlag()) {
                                                R_OCC.add(new Integer(1));
                                                F_OCC.add(new Integer(0));
                                        } else {
                                                F_OCC.add(new Integer(1));
                                                R_OCC.add(new Integer(0));
                                        }
                                }
                        }
                }
        }
        
        public static void dumpExcess(String chrom) {
                int trim = 9000;
                while(trim > 0) {
                        int sum = F_OCC.get(0).intValue() + R_OCC.get(0).intValue();
                        OUT.println(chrom + "\t" + BP.get(0).intValue() + "\t" + F_OCC.get(0).intValue() + "\t" + R_OCC.get(0).intValue() + "\t" + sum);
                        BP.remove(0);
                        F_OCC.remove(0);
                        R_OCC.remove(0);
                        trim--;
                }
        }
        
        public static void processREADS() {
                inputSam = SamReaderFactory.makeDefault().open(BAM);//factory.open(BAM);
                AbstractBAMFileIndex bai = (AbstractBAMFileIndex) inputSam.indexing().getIndex();
                                        
                for(int numchrom = 0; numchrom < bai.getNumberOfReferences(); numchrom++) {
                        SAMSequenceRecord seq = inputSam.getFileHeader().getSequence(numchrom);
                        System.out.println("Processing: " + seq.getSequenceName());

                        CHROMSTOP = seq.getSequenceLength();
                        BP = new ArrayList<Integer>();
                        F_OCC = new ArrayList<Integer>();
                        R_OCC = new ArrayList<Integer>();
                        
                        CloseableIterator<SAMRecord> iter = inputSam.query(seq.getSequenceName(), 0, seq.getSequenceLength(), false);
                        while (iter.hasNext()) {
                                //Create the record object 
                                SAMRecord sr = iter.next();
                                
                                if(STRAND == 2) { //Output combined READ 1 && READ 2
                                        if(PAIR == 0) { addTag(sr); } //Output read if proper mate-pairing is NOT required
                                        else if(sr.getReadPairedFlag()) { //otherwise, check for PE flag
                                                if(sr.getProperPairFlag()) { addTag(sr); } //output read if proper mate-pair is detected
                                        }
                                } else if(STRAND == 0) { //Output READ 1
                                        if(sr.getReadPairedFlag()) { //Check if PAIRED-END
                                                if(((sr.getProperPairFlag() && PAIR == 1) || PAIR == 0) && sr.getFirstOfPairFlag()) { //mate must be mapped if PAIR requirement, must be read1
                                                        boolean flag1 = (Math.abs(sr.getInferredInsertSize()) >= MIN_INSERT && MIN_INSERT != -9999) || MIN_INSERT == -9999; //check if insert size >= min if in use
                                                        boolean flag2 = (Math.abs(sr.getInferredInsertSize()) <= MAX_INSERT && MAX_INSERT != -9999) || MAX_INSERT == -9999; //check if insert size <= max if in use
                                                        if(flag1 && flag2) { addTag(sr); } //add tag if both flags true
                                                }
                                        } else if(PAIR == 0) { addTag(sr); } //Output if not paired-end, by default it is Read1, and mate-pair not required
                                } else if(STRAND == 1) { //Output READ 2
                                        if(sr.getReadPairedFlag()) { ////Must be PAIRED-END for valid Read 2
                                                if(((sr.getProperPairFlag() && PAIR == 1) || PAIR == 0) && !sr.getFirstOfPairFlag()) { //mate must be mapped if PAIR requirement, must be read2
                                                        boolean flag1 = (Math.abs(sr.getInferredInsertSize()) >= MIN_INSERT && MIN_INSERT != -9999) || MIN_INSERT == -9999; //check if insert size >= min if in use
                                                        boolean flag2 = (Math.abs(sr.getInferredInsertSize()) <= MAX_INSERT && MAX_INSERT != -9999) || MAX_INSERT == -9999; //check if insert size <= max if in use
                                                        if(flag1 && flag2) { addTag(sr); } //add tag if both flags true
                                                }
                                        }
                                }
                                
                                //Dump ArrayLists to OUT if they get too big in order to save RAM and therefore time
                                if(BP.size() > 10000) {        dumpExcess(seq.getSequenceName()); }
                                
                        }
                        iter.close();
                        for(int z = 0; z < BP.size(); z++) {
                                int sum = F_OCC.get(z).intValue() + R_OCC.get(z).intValue();
                                OUT.println(seq.getSequenceName() + "\t" + BP.get(z).intValue() + "\t" + F_OCC.get(z).intValue() + "\t" + R_OCC.get(z).intValue() + "\t" + sum);                
                        }
                }
                bai.close();
        }
        
        public static void loadConfig(String[] command){
                for (int i = 0; i < command.length; i++) {
                        switch (command[i].charAt((1))) {
                                case 'b':
                                        BAM = new File(command[i + 1]);
                                        i++;
                                        break;
                                case 'o':
                                        OUTPUT = new File(command[i + 1]);
                                        i++;
                                        break;
                                case 'p':
                                        PAIR = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'r':
                                        STRAND = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'm':
                                        MIN = command[i + 1];
                                        i++;
                                        break;
                                case 'M':
                                        MAX = command[i + 1];
                                        i++;
                                        break;
                        }
                }
                if(BAM == null) {
                        System.out.println("Invalid BAM File!!!\n");
                        printUsage();
                        System.exit(0);
                }
                if(PAIR < 0 || PAIR > 1) {
                        System.out.println("Invalid Mate-Pair requirement!!!\n");
                        printUsage();
                        System.exit(0);                
                }
                if(STRAND > 2 || STRAND < 0) {
                        System.out.println("Invalid Strand Output!!!\n");
                        printUsage();
                        System.exit(0);
                }
                if(STRAND == 0) { READ = "READ1"; }
                else if(STRAND == 1) { READ = "READ2"; }
                else if(STRAND == 2) { READ = "COMBINED"; }
                
                
                if(OUTPUT == null) {
                        OUTPUT = new File(BAM.getName().split("\\.")[0] + "_" + READ + ".scidx");
                }
                                
                System.out.println("-----------------------------------------\nCommand Line Arguments:");
                System.out.println("BAM file: " + BAM);
                System.out.println("Output: " + OUTPUT);
                
                System.out.print("Require proper Mate-pair: ");
                if(PAIR == 0) { System.out.println("no"); }
                else { System.out.println("yes"); }
                
                System.out.println("Output Read: " + READ);
                
                if(MIN.matches("\\d+")) { 
                        MIN_INSERT = Integer.parseInt(MIN);
                        System.out.println("Minimum insert size required to output: " + MIN_INSERT);
                } else { System.out.println("Minimum insert size required to output: NaN"); }
                if(MAX.matches("\\d+")) { 
                        MAX_INSERT = Integer.parseInt(MAX);
                        System.out.println("Maximum insert size required to output: " + MAX_INSERT);
                } else { System.out.println("Maximum insert size required to output: NaN"); }
        }
        
        public static void printUsage() {
                System.out.println("Usage: java -jar BAMtoscIDX.jar -b [BAMFile] -o [OutputFile] [Options]");
                System.out.println("-----------------------------------------");
                System.out.println("Required: BAM file must be sorted and BAI index must be in same folder as BAM file.");
                System.out.println("\nRequired Parameter:");
                System.out.println("BAM File:\t\t-b\t\tBAM file");
                
                System.out.println("\nOptional Parameters:");
                System.out.println("Output File Name:\t-o\t\tOutput file");
                                
                System.out.println("\nOptional Paired-end Parameters:");
                System.out.println("Require proper mate-pairing:\t\t-p\t\t0,1 (1 - yes[Default], 0 - no)");
                System.out.println("Read to Output:\t\t\t\t-r\t\t0,1,2 (Read1 = 0[Default], Read2 = 1, Combined = 2)");
                System.out.println("Minimum Insert Size to Output:\t\t-m\t\tInteger or NaN[Default]");
                System.out.println("Maximum Insert Size to Output:\t\t-M\t\tInteger or NaN[Default]");
                System.out.println("NOTE: Filtering by insert size parameters will NOT filter out single-end Read 1\nunless proper mate-pairing is required!!!");

        }
        
        private static String getTimeStamp() {
                Date date= new Date();
                String time = new Timestamp(date.getTime()).toString();
                return time;
        }
}
