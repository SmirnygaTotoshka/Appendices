import ij.IJ;
import ij.ImagePlus;
import ij.plugin.Converter;
import java.awt.image.DataBufferByte;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import javax.imageio.ImageIO;
import Jama.Matrix;
import Jama.QRDecomposition;

public class Polynomial {
    public static Matrix beta, coefficients;                // the polynomial regression coefficients
    public static double sse;                 // sum of squares due to error
    public static double sst, y, y1, y2;// total sum of squares
    public static int[] threshold, thresh;
    public static double[] coefficients2Proizv;


    public static Matrix PolynomialRegression(int[] x, double[] y, int degree) {
        // @param  x the values of the predictor variable
        //@param  y the corresponding values of the response variable
        // @param  degree the degree of the polynomial to fit
        // @param  variableName the name of the predictor variable
        // @throws IllegalArgumentException if the lengths of the two arrays are not equal
        int n = x.length;
        QRDecomposition qr = null;
        Matrix matrixX = null;
        // in case Vandermonde matrix does not have full rank, reduce degree until it does
        while (true) {
            // build Vandermonde matrix
            double[][] vandermonde = new double[n][degree + 1];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j <= degree; j++) {
                    vandermonde[i][j] = Math.pow(x[i], j);
                }
            }
            matrixX = new Matrix(vandermonde);

            // find least squares solution
            qr = new QRDecomposition(matrixX);
            if (qr.isFullRank()) break;

            // decrease degree and try again
            degree--;
        }
        // create matrix from vector
        Matrix matrixY = new Matrix(y, n);

        // linear regression coefficients
        beta = qr.solve(matrixY);

        // mean of y[] values
        //   double sum = 0.0;
        //  for (int i = 0; i < n; i++)
        //     sum += y[i];
        //  double mean = sum / n;

        // total variation to be accounted for
        //  for (int i = 0; i < n; i++) {
        //      double dev = y[i] - mean;
        //     sst += dev * dev;
        //   }

        // variation not accounted for
        //   Matrix residuals = matrixX.times(beta).minus(matrixY);
        //  sse = residuals.norm2() * residuals.norm2();
        return beta;
    }

    public static double beta(int j) {
        // Returns the {@code j}th regression coefficient.
        // to make -0.0 print as 0.0
        if (Math.abs(beta.get(j, 0)) < 1E-4) return 0.0;
        return beta.get(j, 0);
    }

    //public double predict(double x) {
    //Returns the expected response {@code y} given the value of the predictor variable {@code x}
    // horner's method
    //  double y = 0.0;
    //for (int j = degree; j >= 0; j--)
    //   y = beta(j) + (x * y);
    //return y;
    //}

    public static double R2() {
        if (sst == 0.0) return 1.0;   // constant function
        return 1.0 - sse / sst;
    }

    public static int[] thresholdCoord(double[][] b){
        int i = 0;
        coefficients2Proizv = new double[15];
        for (int j=2; j<b.length; j++) {
            coefficients2Proizv[i] = b[j][0]*j*(j-1);
       //     System.out.print( " " +b[j][0]);
       //     System.out.print( " " +coefficients2Proizv[i] + b[j][0]*j*(j-1));
            i++;
        }
        thresh = new int[256];
        double okruglenie = Math.pow(10,10);
        for (int j =0; j<257; j++) {
            y = coefficients2Proizv[0] + coefficients2Proizv[1]*j;
                    //+coefficients2Proizv[2]*(j^2)+coefficients2Proizv[3]*(j^3)+coefficients2Proizv[4]*(j^4)+coefficients2Proizv[5]*(j^5)+coefficients2Proizv[6]*(j^6)+coefficients2Proizv[7]*(j^7)+coefficients2Proizv[8]*(j^8)+coefficients2Proizv[9]*(j^9)+coefficients2Proizv[10]*(j^10)+coefficients2Proizv[11]*(j^11)+coefficients2Proizv[12]*(j^12)+coefficients2Proizv[13]*(j^13);
            y = Math.ceil(y*okruglenie)/okruglenie;
            if (y==0) {
             //   y1 = coefficients2Proizv[0] + coefficients2Proizv[1]*(j-1)+coefficients2Proizv[2]*((j-1)^2)+coefficients2Proizv[3]*((j-1)^3)+coefficients2Proizv[4]*((j-1)^4)+coefficients2Proizv[5]*((j-1)^5)+coefficients2Proizv[6]*((j-1)^6)+coefficients2Proizv[7]*((j-1)^7)+coefficients2Proizv[8]*((j-1)^8)+coefficients2Proizv[9]*((j-1)^9)+coefficients2Proizv[10]*((j-1)^10)+coefficients2Proizv[11]*((j-1)^11)+coefficients2Proizv[12]*((j-1)^12)+coefficients2Proizv[13]*((j-1)^13);
             //   y2 = coefficients2Proizv[0] + coefficients2Proizv[1]*(j+1)+coefficients2Proizv[2]*((j+1)^2)+coefficients2Proizv[3]*((j+1)^3)+coefficients2Proizv[4]*((j+1)^4)+coefficients2Proizv[5]*((j+1)^5)+coefficients2Proizv[6]*((j+1)^6)+coefficients2Proizv[7]*((j+1)^7)+coefficients2Proizv[8]*((j+1)^8)+coefficients2Proizv[9]*((j+1)^9)+coefficients2Proizv[10]*((j+1)^10)+coefficients2Proizv[11]*((j+1)^11)+coefficients2Proizv[12]*((j+1)^12)+coefficients2Proizv[13]*((j+1)^13);
             //   if (y1>0 && y2<0) {thresh[i]=j; i++;}
             //   if (y1<0 && y2>0) {thresh[i]=j; i++;}
                {thresh[i]=j; i++;};
            }
        }
        return thresh;
    }

    public static void main(String[] args) throws IOException {
        IJ.open("C:\\Users\\maria\\SynologyDrive\\Appendix\\Tiff\\0012_Wholeslide_Default_Extended_tif_1_1024_76800.tif");
        ImagePlus imp1 = IJ.getImage();
        ImagePlus dup = imp1.duplicate();
        Converter converter = new Converter();
        converter.run("8-bit");
        IJ.save("C:\\Users\\maria\\SynologyDrive\\Appendix\\Tiff\\0012_Wholeslide_Default_Extended_tif_1_1024_76800_8bit.tif");

        // getting massive of bytes from image
        BufferedImage image = ImageIO.read(new File("C:\\Users\\maria\\SynologyDrive\\appendix\\Tiff\\0012_Wholeslide_Default_Extended_tif_1_1024_76800_8bit.tif"));
        byte[] pixels = ((DataBufferByte) image.getRaster().getDataBuffer()).getData();

        // count byte's frequencies
        Map<Integer, Double> counter = new HashMap<>();
        for (int x : pixels) {
            double newValue = counter.getOrDefault(x, 0.0) + 1;
            counter.put(x, newValue);
        }
        //System.out.println(counter);
       // counter.forEach((pix, count) -> {
        //    count = count / 1024 / 1024;
        //    counter.put(pix, count);
      //  });

        // sorting
        TreeMap<Integer, Double> sortedMap = new TreeMap<>(counter);
        //System.out.println(sortedMap);

        int[] key = new int[sortedMap.size()];
        double[] object = new double[sortedMap.size()];
        int n;
        n = 0;
        // getting arrays for function: key - pixels, object - weighted mean
        for (Map.Entry entry : sortedMap.entrySet()) {
            key[n] = (int) entry.getKey()+128;
            object[n] = (double) entry.getValue();
            n++;
        }
    //    for (int i = 0; i < key.length; i++) {
    //        System.out.print( " " + key[i]);
    //    }

        // function
        double chislitel;
        chislitel = 0;
        double znamenatel;
        znamenatel = 0;
        double[] function = new double[key.length];
        for (int i = 0; i < key.length; i++) {
            chislitel = chislitel + key[i] * object[i];
            znamenatel = znamenatel + object[i];
            function[i] = chislitel / znamenatel;
        }

      //  for (int i = 0; i < function.length; i++) {
      //      System.out.print( " " + function[i]);
      //  }

        // the least squares method
        coefficients = PolynomialRegression(key, function, 15);
        double[][] coefficients_array = coefficients.getArray();
   //     for (int i = 0; i < coefficients_array.length; i++) {
     //       for(int j = 0; j < coefficients_array[i].length; j++) {
       //         System.out.print( " " + coefficients_array[i][j] );
        //    }
      //  }

        //вывод уравнения на экран
       // StringBuilder s = new StringBuilder();
       // int j = 15;
        //ignoring leading zero coefficients
       // while (j >= 0 && Math.abs(beta(j)) < 1E-5)
       //     j--;
        // create remaining terms
      //  while (j >= 0) {
       //     if (j == 0) s.append(String.format("%.2f ", beta(j)));
        //    else if (j == 1) s.append(String.format("%.2f %s + ", beta(j), "x"));
        //    else s.append(String.format("%.2f %s^%d + ", beta(j), "x", j));
        //    j--;
      //  }
      //  s = s.append("  (R^2 = " + String.format("%.3f", R2()) + ")");
        //replace "+ -2n" with "- 2n"
      //  System.out.printf(s.toString().replace("+ -", "- "));

        //координаты для порогов
        threshold = new int[256];
        threshold = thresholdCoord(coefficients_array);
     //   for (int i = 0; i < threshold.length; i++) {
     //       System.out.print( " " + threshold[i]);
     //   }
    }
}
