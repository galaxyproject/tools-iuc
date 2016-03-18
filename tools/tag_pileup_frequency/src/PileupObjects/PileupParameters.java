package PileupObjects;

import java.awt.Color;
import java.io.File;
import java.util.ArrayList;

public class PileupParameters {
        private File OUTPUT = null;
        private String COMPOSITE = null;
        private int READ = 0;
        private int STRAND = 0;
        private int TRANS = 0;
        private int SHIFT = 0;
        private int BIN = 1;
        private int SMOOTH = 0;
        private int STDSIZE = 0;
        private int STDNUM = 0;
        private int CPU = 1;
        private int OUTTYPE = 0;
        private boolean STANDARD = false;
        private boolean outputCOMPOSITE = false;
        private boolean requirePE = false;
        private double STANDRATIO = 1;
        
        private Color Sense = null;
        private Color Anti = null;
        private Color Combined = null;
        
        public PileupParameters() {
                
        }

        public boolean getPErequire() {
                return requirePE;
        }
        
        public void setPErequire(boolean status) {
                requirePE = status;
        }
        
        public String getCompositeFile() {
                return COMPOSITE;
        }
        
        public void setCompositeFile(String comp) {
                COMPOSITE = comp;
        }
        
        public boolean getOutputCompositeStatus() {
                return outputCOMPOSITE;
        }
        
        public void setOutputCompositeStatus(boolean out) {
                outputCOMPOSITE = out;
        }
        
        public ArrayList<Color> getColors() {
                ArrayList<Color> ALL = new ArrayList<Color>();
                if(Sense != null) ALL.add(Sense);
                if(Anti != null) ALL.add(Anti);
                if(Combined != null) ALL.add(Combined);
                return ALL;
        }
        
        public void setSenseColor(Color s) {
                Sense = s;
        }
        
        public Color getSenseColor() {
                return Sense;
        }
        
        public void setAntiColor(Color a) {
                Anti = a;
        }
        
        public Color getAntiColor() {
                return Anti;
        }
        
        public void setCombinedColor(Color c) {
                Combined = c;
        }
        
        public Color getCombinedColor() {
                return Combined;
        }
        
        public void setRatio(double rat) {
                STANDRATIO = rat;
        }
        
        public double getRatio() {
                return STANDRATIO;
        }
        
        public void setStandard(boolean stand) {
                STANDARD = stand;
        }
        
        public boolean getStandard() {
                return STANDARD;
        }

        public void setOutputType(int newtype) {
                OUTTYPE = newtype;
        }
        
        public int getOutputType() {
                return OUTTYPE;
        }
        
        public File getOutput() {
                return OUTPUT;
        }

        public void setOutput(File oUTPUT) {
                OUTPUT = oUTPUT;
        }

        public int getRead() {
                return READ;
        }

        public void setRead(int rEAD) {
                READ = rEAD;
        }

        public int getStrand() {
                return STRAND;
        }

        public void setStrand(int sTRAND) {
                STRAND = sTRAND;
        }

        public int getTrans() {
                return TRANS;
        }

        public void setTrans(int tRANS) {
                TRANS = tRANS;
        }

        public int getShift() {
                return SHIFT;
        }

        public void setShift(int sHIFT) {
                SHIFT = sHIFT;
        }

        public int getBin() {
                return BIN;
        }

        public void setBin(int bIN) {
                BIN = bIN;
        }

        public int getSmooth() {
                return SMOOTH;
        }

        public void setSmooth(int sMOOTH) {
                SMOOTH = sMOOTH;
        }

        public int getStdSize() {
                return STDSIZE;
        }

        public void setStdSize(int sTDSIZE) {
                STDSIZE = sTDSIZE;
        }

        public int getStdNum() {
                return STDNUM;
        }

        public void setStdNum(int sTDNUM) {
                STDNUM = sTDNUM;
        }

        public int getCPU() {
                return CPU;
        }

        public void setCPU(int cPU) {
                CPU = cPU;
        }
        
}
