import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Scanner;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import PileupObjects.BEDCoord;
import PileupObjects.PileupParameters;
import PileupScripts.BAMUtilities;
import PileupScripts.PileupExtract;
import PileupScripts.TransformArray;
import htsjdk.samtools.AbstractBAMFileIndex;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class TagPileup {
        //Global Transformations
        private static boolean TAGSEQUAL = false; //set tags to be equal
        private static int SHIFT = 0; //shift from 5' to 3'
        private static int BINSIZE = 1; //bin tags, useful for mammalian or larger windows
        
        //Composite Transformations
        private static int SMOOTHWINDOW = 0;
        
        //Read parameters
        private static int STRAND = 0; //0 - strand sep, 1 - strand combined
        private static int READ = 0; //0 - read 1, 1 - read 2, 2 - combined
        private static boolean PE_STATUS = false; //require read 1 and 2 to be properly paired
        
        //Run parameters
        private static int CPU = 1;
        private static boolean OUTPUT_HEATMAP = true;
        private static boolean OUTPUT_COMP = true;
        
        private static File BAM = null;
        private static File BAI = null;
        private static File BED = null;
        private static Vector<BEDCoord> INPUT = null;
        
        private static String HPATH = null;
        private static String CNAME = null;
        private static PileupParameters PARAM = null;
        private static PrintStream COMPOSITE = null;
        private static PrintStream OUT_S1 = null;
        private static PrintStream OUT_S2 = null;
        
        public static void main(String[] args) throws IOException {
                loadConfig(args);

                System.out.println("\n" + getTimeStamp() + "\nLoading BED Coordinates..");
                INPUT = loadCoord(BED);
                System.out.println("BED Coordinates Loaded\n" + getTimeStamp());
                
                //Validate BED coordinates exist in BAM file
                validateBED();
                
                System.out.println("\n" + getTimeStamp() + "\nParsing Tags...");
                run();
                System.out.println("Tags Parsed\n" + getTimeStamp());
        }
        
        //Get size of largest array for composite generation
        public static int getMaxBEDSize(Vector<BEDCoord> sites) {
                int maxSize = 0;
                for(int x = 0; x < sites.size(); x++) {
                        int SIZE = sites.get(x).getFStrand().length;
                        if(SIZE > maxSize) {
                                maxSize = SIZE;
                        }
                }
                return maxSize;
        }
        
        //Validate BED coordinates exist within BAM file and satisfy BED format
        public static void validateBED() throws IOException {
                //Get chromosome IDs
                SamReader inputSam = SamReaderFactory.makeDefault().open(BAM);
                AbstractBAMFileIndex bai = (AbstractBAMFileIndex) inputSam.indexing().getIndex();
                ArrayList<String> chrom = new ArrayList<String>();
                for(int x = 0; x < bai.getNumberOfReferences(); x++) {
                        chrom.add(inputSam.getFileHeader().getSequence(x).getSequenceName());
                }
                inputSam.close();
                bai.close();
                
                ArrayList<Integer> indexFail = new ArrayList<Integer>();
                int FAILcount = 0;
                for(int x = 0; x < INPUT.size(); x++) {
                        //check for existence of chromosome in BAM file
                        if(!chrom.contains(INPUT.get(x).getChrom())) {
                                INPUT.get(x).setStatus(false);
                                indexFail.add(new Integer(x));
                                FAILcount++;
                        }
                        //check that the start is smaller than the stop
                        if(INPUT.get(x).getStart() > INPUT.get(x).getStop()) {
                                INPUT.get(x).setStatus(false);
                                indexFail.add(new Integer(x));
                                FAILcount++;
                        }
                }
                if(FAILcount == INPUT.size()) {
                        System.err.println("No BED Coordinates exist within BAM file!!!");
                        printUsage();
                        System.exit(1);
                }
                //Remove failed indexes to more efficiently use CPUs
                for(int x = indexFail.size() - 1; x >= 0; x--) {
                        INPUT.remove(indexFail.get(x).intValue());
                }
        }
        
        public static void run() throws IOException {
                if(PARAM.getOutputCompositeStatus()) {
                        try { 
                                if(CNAME != null) { COMPOSITE = new PrintStream(CNAME); }
                                else { COMPOSITE = new PrintStream("composite.out"); }
                        } catch (FileNotFoundException e) {        e.printStackTrace(); }
                }
                
                //Code to standardize tags sequenced to genome size (1 tag / 1 bp)
                if(PARAM.getStandard()) { PARAM.setRatio(BAMUtilities.calculateStandardizationRatio(BAM)); }
                                                        
                        if(PARAM.getOutputType() != 0) {
                                if(STRAND == 0) {
                                        try {
                                                if(HPATH != null) {
                                                        OUT_S1 = new PrintStream(HPATH + File.separator + generateFileName(BED.getName(), BAM.getName(), 0));
                                                        OUT_S2 = new PrintStream(HPATH + File.separator + generateFileName(BED.getName(), BAM.getName(), 1));
                                                } else {
                                                        OUT_S1 = new PrintStream(generateFileName(BED.getName(), BAM.getName(), 0));
                                                        OUT_S2 = new PrintStream(generateFileName(BED.getName(), BAM.getName(), 1));
                                                }                                                
                                        } catch (FileNotFoundException e) {        e.printStackTrace(); }
                                } else {
                                        try {
                                                if(HPATH != null) { OUT_S1 = new PrintStream(HPATH + File.separator + generateFileName(BED.getName(), BAM.getName(), 2)); }
                                                else { OUT_S1 = new PrintStream(generateFileName(BED.getName(), BAM.getName(), 2)); }
                                        } catch (FileNotFoundException e) {        e.printStackTrace(); }
                                }
                        }
                                        
                        //Split up job and send out to threads to process                                
                        ExecutorService parseMaster = Executors.newFixedThreadPool(CPU);
                        if(INPUT.size() < CPU) CPU = INPUT.size();
                        int subset = 0;
                        int currentindex = 0;
                        for(int x = 0; x < CPU; x++) {
                                currentindex += subset;
                                if(CPU == 1) subset = INPUT.size();
                                else if(INPUT.size() % CPU == 0) subset = INPUT.size() / CPU;
                                else {
                                        int remainder = INPUT.size() % CPU;
                                        if(x < remainder ) subset = (int)(((double)INPUT.size() / (double)CPU) + 1);
                                        else subset = (int)(((double)INPUT.size() / (double)CPU));
                                }
                                PileupExtract extract = new PileupExtract(PARAM, BAM, INPUT, currentindex, subset);
                                parseMaster.execute(extract);
                        }
                        parseMaster.shutdown();
                        while (!parseMaster.isTerminated()) {
                        }
                                        
                        //double[] AVG_S1 = new double[INPUT.get(0).getFStrand().length];
                        double[] AVG_S1 = new double[getMaxBEDSize(INPUT)];
                        double[] AVG_S2 = null;
                        if(STRAND == 0) AVG_S2 = new double[AVG_S1.length];
                        double[] DOMAIN = new double[AVG_S1.length];
        
                        //Account for the shifted oversized window produced by binning and smoothing
                        int OUTSTART = 0;
                        if(PARAM.getTrans() == 1) { OUTSTART = PARAM.getSmooth(); }
                        else if(PARAM.getTrans() == 2) { OUTSTART = (PARAM.getStdSize() * PARAM.getStdNum()); }
                                                                        
                        if(PARAM.getOutputType() == 2) {
                                if(OUT_S1 != null) OUT_S1.print("YORF\tNAME");
                                if(OUT_S2 != null) OUT_S2.print("YORF\tNAME");
                                                                                        
                                for(int i = OUTSTART; i < AVG_S1.length - OUTSTART; i++) {
                                        int index = i - OUTSTART;
                                        if(OUT_S1 != null) OUT_S1.print("\t" + index);
                                        if(OUT_S2 != null) OUT_S2.print("\t" + index);
                                }
                                if(OUT_S1 != null) OUT_S1.println();
                                if(OUT_S2 != null) OUT_S2.println();
                        }
                                        
                        //Output individual sites
                        for(int i = 0; i < INPUT.size(); i++) {
                                double[] tempF = INPUT.get(i).getFStrand();
                                double[] tempR = INPUT.get(i).getRStrand();
                                if(OUT_S1 != null) OUT_S1.print(INPUT.get(i).getName());
                                if(OUT_S2 != null) OUT_S2.print(INPUT.get(i).getName());
                                
                                if(PARAM.getOutputType() == 2) {
                                        if(OUT_S1 != null) OUT_S1.print("\t" + INPUT.get(i).getName());
                                        if(OUT_S2 != null) OUT_S2.print("\t" + INPUT.get(i).getName());
                                }
                                
                                for(int j = 0; j < tempF.length; j++) {
                                        if(j >= OUTSTART && j < tempF.length - OUTSTART) {
                                                if(OUT_S1 != null) OUT_S1.print("\t" + tempF[j]);
                                                if(OUT_S2 != null) OUT_S2.print("\t" + tempR[j]);
                                        }
                                        AVG_S1[j] += tempF[j];
                                        if(AVG_S2 != null) AVG_S2[j] += tempR[j];
                                }
                                if(OUT_S1 != null) OUT_S1.println();
                                if(OUT_S2 != null) OUT_S2.println();
                        }
        
                        //Calculate average and domain here
                        int temp = (int) (((double)AVG_S1.length / 2.0) + 0.5);
                        for(int i = 0; i < AVG_S1.length; i++) {
                                DOMAIN[i] = (double)((temp - (AVG_S1.length - i)) * PARAM.getBin()) + 1;
                                AVG_S1[i] /= INPUT.size();
                                if(AVG_S2 != null) AVG_S2[i] /= INPUT.size();
                        }
                        
                        //Transform average given transformation parameters
                        if(PARAM.getTrans() == 1) { 
                                AVG_S1 = TransformArray.smoothTran(AVG_S1, PARAM.getSmooth());
                                if(AVG_S2 != null) AVG_S2 = TransformArray.smoothTran(AVG_S2, PARAM.getSmooth());
                        } else if(PARAM.getTrans() == 2) {
                                AVG_S1 = TransformArray.gaussTran(AVG_S1, PARAM.getStdSize(), PARAM.getStdNum());
                                if(AVG_S2 != null) AVG_S2 = TransformArray.gaussTran(AVG_S2, PARAM.getStdSize(), PARAM.getStdNum());
                        }
                                        
                        //Trim average here and output to statistics pane
                        double[] AVG_S1_trim = new double[AVG_S1.length - (OUTSTART * 2)];
                        double[] AVG_S2_trim = null;
                        if(STRAND == 0) AVG_S2_trim = new double[AVG_S1_trim.length];
                        double[] DOMAIN_trim = new double[AVG_S1_trim.length];
                        for(int i = OUTSTART; i < AVG_S1.length - OUTSTART; i++) {
                                if(AVG_S2 != null) {
                                        AVG_S2_trim[i - OUTSTART] = AVG_S2[i];
                                }
                                AVG_S1_trim[i - OUTSTART] = AVG_S1[i];
                                DOMAIN_trim[i - OUTSTART] = DOMAIN[i];
                        }
                        AVG_S1 = AVG_S1_trim;
                        AVG_S2 = AVG_S2_trim;
                        DOMAIN = DOMAIN_trim;
                                        
                        //Output composite data to tab-delimited file
                        if(COMPOSITE != null) {
                                for(int a = 0; a < DOMAIN.length; a++) {
                                        COMPOSITE.print("\t" + DOMAIN[a]);
                                }
                                COMPOSITE.println();
                                if(STRAND == 0) {
                                        COMPOSITE.print(generateFileName(BED.getName(), BAM.getName(), 0));
                                        for(int a = 0; a < AVG_S1.length; a++) {
                                                COMPOSITE.print("\t" + AVG_S1[a]);
                                        }
                                        COMPOSITE.println();
                                        COMPOSITE.print(generateFileName(BED.getName(), BAM.getName(), 1));
                                        for(int a = 0; a < AVG_S2.length; a++) {
                                                COMPOSITE.print("\t" + AVG_S2[a]);
                                        }
                                        COMPOSITE.println();
                                } else {
                                        COMPOSITE.print(generateFileName(BED.getName(), BAM.getName(), 2));
                                        for(int a = 0; a < AVG_S1.length; a++) {
                                                COMPOSITE.print("\t" + AVG_S1[a]);
                                        }
                                        COMPOSITE.println();
                                }                        
                        }
                                                                                
                        if(OUT_S1 != null && PARAM.getOutputType() == 2) { OUT_S1.close(); }
                        if(OUT_S2 != null && PARAM.getOutputType() == 2) { OUT_S2.close(); }
        }
        
        public static String generateFileName(String bed, String bam, int strandnum) {
                String[] bedname = bed.split("\\.");
                String[] bamname = bam.split("\\.");
                
                String strand = "sense";
                if(strandnum == 1) strand = "anti";
                else if(strandnum == 2) strand = "combined";
                String read = "read1";
                if(PARAM.getRead() == 1) read = "read2";
                else if(PARAM.getRead() == 2) read = "readc";
                
                String filename = bedname[0] + "_" + bamname[0] + "_" + read + "_" + strand + ".tabular";
                return filename;
        }
        
    public static Vector<BEDCoord> loadCoord(File INPUT) throws FileNotFoundException {
                Scanner scan = new Scanner(INPUT);
                Vector<BEDCoord> c = new Vector<BEDCoord>();
                while (scan.hasNextLine()) {
                        String[] temp = scan.nextLine().split("\t");
                        if(temp.length > 2) { 
                                if(!temp[0].contains("track") && !temp[0].contains("#")) {
                                        String name = "";
                                        if(temp.length > 3) { name = temp[3]; }
                                        else { name = temp[0] + "_" + temp[1] + "_" + temp[2]; }
                                        if(Integer.parseInt(temp[1]) >= 0) {
                                                if(temp.length > 4) { 
                                                        if(temp[5].equals("+")) { c.add(new BEDCoord(temp[0], Integer.parseInt(temp[1]), Integer.parseInt(temp[2]), "+", name)); }
                                                        else { c.add(new BEDCoord(temp[0], Integer.parseInt(temp[1]), Integer.parseInt(temp[2]), "-", name)); }
                                                } else { c.add(new BEDCoord(temp[0], Integer.parseInt(temp[1]), Integer.parseInt(temp[2]), "+", name)); }

                                        } else {
                                                System.out.println("Invalid Coordinate in File!!!\n" + Arrays.toString(temp));
                                        }
                                }
                        }
                }
                scan.close();
                return c;
    }

        public static void loadConfig(String[] command){
                for (int i = 0; i < command.length; i++) {
                        switch (command[i].charAt((1))) {
                                case 'b':
                                        BAM = new File(command[i + 1]);
                                        i++;
                                        break;
                                case 'i':
                                        BAI = new File(command[i + 1]);
                                        i++;
                                        break;
                                case 'c':
                                        BED = new File(command[i + 1]);
                                        i++;
                                        break;
                                case 's':
                                        SHIFT = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'n':
                                        BINSIZE = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'e':
                                        TAGSEQUAL = Boolean.parseBoolean(command[i + 1]); 
                                        i++;
                                        break;
                                case 'r':
                                        READ = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'a':
                                        STRAND = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'p':
                                        PE_STATUS = Boolean.parseBoolean(command[i + 1]);
                                        i++;
                                        break;
                                case 't':
                                        CPU = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'w':
                                        SMOOTHWINDOW = Integer.parseInt(command[i + 1]);
                                        i++;
                                        break;
                                case 'h':
                                        OUTPUT_HEATMAP = Boolean.parseBoolean(command[i + 1]);
                                        i++;
                                        break;
                                case 'm':
                                        OUTPUT_COMP = Boolean.parseBoolean(command[i + 1]);
                                        i++;
                                        break;
                                case 'o':
                                        HPATH = command[i + 1];
                                        i++;
                                        break;
                                case 'x':
                                        CNAME = command[i + 1];
                                        i++;
                                        break;
                        }
                }
                if(BAM == null) {
                        System.err.println("No BAM File loaded!!!\n");
                        printUsage();
                        System.exit(1);
                } else if(BAI == null) {
                        System.err.println("No BAI Index File loaded!!!\n");
                        printUsage();
                        System.exit(1);
                }else if(BED == null) {
                        System.err.println("No BED File loaded!!!\n");
                        printUsage();
                        System.exit(1);
                }
                
                if(BINSIZE < 1) {
                        System.err.println("Invalid Bin Size!!! Must be larger than 0 bp");
                        printUsage();
                        System.exit(1);
                } else if(READ < 0 || READ > 2) {
                        System.err.println("Invalid Read to Examine!!! Select 0-2");
                        printUsage();
                        System.exit(1);
                } else if(CPU < 1) {
                        System.err.println("Invalid Number of CPU's!!! Must use at least 1");
                        printUsage();
                        System.exit(1);
                } else if(SMOOTHWINDOW < 0) {
                        System.err.println("Invalid Smoothing Window Size!!! Must be larger than 0 bp");
                        printUsage();
                        System.exit(1);
                } 

                //Load up parameters for the pileup into single object
                PARAM = new PileupParameters();
                               
            //Initialize global transformation parameters
            PARAM.setShift(SHIFT); //SHIFT can be negative
            PARAM.setBin(BINSIZE);
            PARAM.setStandard(TAGSEQUAL);
            
            //Initialize read parameters
            PARAM.setRead(READ);
            PARAM.setPErequire(PE_STATUS);
            if(STRAND == 0) {
                    PARAM.setStrand(0);
                    PARAM.setSenseColor(Color.BLUE);
                    PARAM.setAntiColor(Color.RED);
            } else {
                    PARAM.setStrand(1);
                    PARAM.setCombinedColor(new Color(0, 100, 0));
            }
            
            //Initialize run parameters
            PARAM.setCPU(CPU);
            if(SMOOTHWINDOW > 0) {
                    PARAM.setTrans(1);
                    PARAM.setSmooth(SMOOTHWINDOW);
            } else { PARAM.setTrans(0); }
            
            if(!OUTPUT_HEATMAP) { PARAM.setOutputType(0); }
            else { PARAM.setOutputType(2); }
            if(!OUTPUT_HEATMAP && !OUTPUT_COMP) { PARAM.setOutput(null); }
            else { PARAM.setOutput(new File(System.getProperty("user.dir"))); }
            
            PARAM.setOutputCompositeStatus(OUTPUT_COMP);
            
                System.out.println("-----------------------------------------\nCommand Line Arguments:");
                System.out.println("BAM file: " + BAM);
                System.out.println("BAI file: " + BAI);
                System.out.println("BED file: " + BED);
                
                System.out.println("\nGlobal transformations:");
                System.out.println("5' to 3' Tag Shift (bp):\t" + SHIFT);
                System.out.println("Bin Size (bp):\t\t\t" + BINSIZE);
                System.out.println("Set Tags to be equal:\t\t" + TAGSEQUAL);
                System.out.println("\nRead parameters:");
                System.out.println("Reads to Examine:\t\t" + READ);
                System.out.println("Require Proper PE:\t\t" + PE_STATUS);
                System.out.println("Combine Strands:\t\t" + STRAND);
                System.out.println("\nRun parameters:");
                System.out.println("CPUs to use:\t\t\t" + CPU);
                System.out.println("Composite smoothing window (bp):" + SMOOTHWINDOW);
                System.out.println("Output heatmap:\t\t\t" + OUTPUT_HEATMAP);
                System.out.println("Output path for heatmap:\t" + HPATH);
                System.out.println("Output composite:\t\t" + OUTPUT_COMP);
                System.out.println("Output filename for composite:\t" + CNAME);

        }
        
        public static void printUsage() {
                System.err.println("Usage: java -jar TagPileup.jar -b [BAMFile] -c [BEDFile] [Options]");
                System.err.println("-----------------------------------------");
                System.err.println("\nRequired Parameter:");
                System.err.println("BAM File:\t\t\t-b\tBAM file");
                System.err.println("BAI File:\t\t\t-i\tBAI file");
                System.err.println("BED File:\t\t\t-c\tBED file");
                System.err.println("\nOptional Parameters:");
                System.err.println("\nGlobal transformations:");
                System.err.println("5' to 3' Tag Shift (bp):\t-s\t0 (default)");
                System.err.println("Bin Size (bp):\t\t\t-n\t1 (default)");
                System.err.println("Set Tags to be equal:\t\t-e\tfalse (default), true");
                System.err.println("\nRead parameters:");
                System.err.println("Reads to Examine:\t\t-r\t0 - Read 1 (default), 1 - Read 2, 2 - Combined");
                System.err.println("Require Proper PE:\t\t-p\tfalse (default), true");
                System.err.println("Combine Strands:\t\t-a\t0 - Strand seperate (default), 1 - Strand combined");
                System.err.println("\nRun parameters:");
                System.err.println("CPUs to use:\t\t\t-t\t1 (default)");
                System.err.println("Composite smoothing window (bp):-w\t0 (default)");
                System.err.println("Output heatmap:\t\t\t-h\ttrue (default), false");
                System.err.println("Output composite:\t\t-m\ttrue (default), false");
                System.err.println("Output path for heatmap:\t-o\tDefault local path");
                System.err.println("Output filename for composite:\t-x\tcomposite.out (default)");
        }
        
        private static String getTimeStamp() {
                Date date= new Date();
                String time = new Timestamp(date.getTime()).toString();
                return time;
        }
}
