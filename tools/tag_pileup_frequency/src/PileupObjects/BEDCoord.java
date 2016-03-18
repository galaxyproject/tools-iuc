package PileupObjects;

import java.util.Comparator;

public class BEDCoord {
        private String CHROM = "";
        private int START = 0;
        private int STOP = 0;
        private String DIR = "+";
        
        private String NAME = ".";
        private double SCORE = 0;

        private double[] Fstrand = null;
        private double[] Rstrand = null;
        
        private boolean VALIDCHROM = true;
        
        public BEDCoord() {
                
        }
        
        public BEDCoord(String c, int sta, int sto, String di, String na) {
                CHROM = c;
                START = sta;
                STOP = sto;
                DIR = di;
                NAME = na;
        }
        
        public BEDCoord(String na, double sco) {
                NAME = na;
                SCORE = sco;
        }
        
        public boolean getStatus() {
                return VALIDCHROM;
        }
        
        public void setStatus(boolean stat) {
                VALIDCHROM = stat;
        }
        
        public double[] getFStrand() {
                return Fstrand;
        }
        
        public double[] getRStrand() {
                return Rstrand;
        }
        
        public void setFstrand(double[] f) {
                Fstrand = f;
        }
        
        public void setRstrand(double[] r) {
                Rstrand = r;
        }
        
        public String getChrom() {
                return CHROM;
        }
        
        public void setChrom(String chr) {
                CHROM = chr;
        }
        
        public int getStart() {
                return START;
        }
        
        public void setStart(int sta) {
                START = sta;
        }
        
        public int getStop() {
                return STOP;
        }
        
        public void setStop(int sto) {
                STOP = sto;
        }
        
        public String getDir() {
                return DIR;
        }
        
        public void setDir(String di) {
                DIR = di;
        }
        
        public String getName() {
                return NAME;
        }
        
        public void setName(String na) {
                NAME = na;
        }
        
        public double getScore() {
                return SCORE;
        }
        
        public void setScore(double sco) {
                SCORE = sco;
        }
        
        public String toString() {
                String line = CHROM + "\t" + START + "\t" + STOP + "\t" + NAME + "\t" + SCORE + "\t" + DIR;
                return line;
        }
        
        public static Comparator<BEDCoord> PeakChromComparator = new Comparator<BEDCoord>() {
                public int compare(BEDCoord node1, BEDCoord node2) {
                        return node1.getChrom().compareTo(node2.getChrom());
        }
        };
                
        public static Comparator<BEDCoord> PeakPositionComparator = new Comparator<BEDCoord>() {
                public int compare(BEDCoord node1, BEDCoord node2) {
                        int PeakStart1 = node1.getStart();
                        int PeakStart2 = node2.getStart();
                        if (PeakStart1 > PeakStart2) return 1;
                        else if (PeakStart1 < PeakStart2) return -1;
                        else return 0;
        }
        };
        
        public static Comparator<BEDCoord> ScoreComparator = new Comparator<BEDCoord>() {
                public int compare(BEDCoord node1, BEDCoord node2) {
                        double PeakStart1 = node1.getScore();
                        double PeakStart2 = node2.getScore();
                        if (PeakStart1 < PeakStart2) return 1;
                        else if (PeakStart1 > PeakStart2) return -1;
                        else return 0;
        }
        };
}
