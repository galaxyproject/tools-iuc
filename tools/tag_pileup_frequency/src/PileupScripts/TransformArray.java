package PileupScripts;

public class TransformArray {
        
        public static double[] reverseTran(double[] orig) {
                if(orig != null) {
                        double[] reverse = new double[orig.length];                
                        for(int x = orig.length - 1; x >= 0; x--) {
                                reverse[orig.length - 1 - x] = orig[x];
                        }
                        return reverse;
                }
                return null;
        }
        
        public static double[] smoothTran(double[] orig, int win) {
                int window = win / 2;
                double[] Sarray = new double[orig.length];
                for(int x = 0; x < orig.length; x++) {
                        double score = 0;
                        double weight = 0;
                        for(int y = x - window; y <= x + window; y++) {
                                if(y < 0) y = -1;
                                else if(y < orig.length) {
                                        score += orig[y];
                                        weight++;
                                } else y = x + window + 1;
                        }
                        if(weight != 0) Sarray[x] = score / weight;
                }
                return Sarray;
        }
        
        public static double[] gaussTran(double[] orig, int size, int num) {
                double[] Garray = new double[orig.length];
                int window = size * num;
                double SDSize = (double)size;
                for(int x = 0; x < orig.length; x++) {
                 double score = 0;
                 double weight = 0;
                 for(int y = x - window; y <= x + window; y++) {
                        if(y < 0) y = -1;
                        else if(y < orig.length) {
                                double HEIGHT = Math.exp(-1 * Math.pow((y - x), 2) / (2 * Math.pow(SDSize, 2)));
                                score += (HEIGHT * orig[y]);
                                weight += HEIGHT;
                        } else y = x + window + 1;
                 }
                 if(weight != 0) Garray[x] = score / weight;
                }
                return Garray;
         }
}
