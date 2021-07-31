/**
 * Appendices is program based on the ImageJ. version 1.0.
 * Authors: Anton S. Smirnov, Maria P. Myshkina, Natalya A. Ermakova - Department of Bioinformatics, PRNRMU
 * Appendices is modificated SlideJ plugin for analysis whole-slide tiff snapshots of appendices.
 * It count number of lymphocytes in muscular tunic. It needs for automatic diagnostic of appendicitis.
 * Rewrite macros written by Myshkina and Ermakova.
 * */
import com.opencsv.CSVReader;
import com.opencsv.CSVWriter;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.ChannelSplitter;
import ij.plugin.Converter;
import ij.plugin.Thresholder;
import ij.plugin.filter.Analyzer;
import ij.plugin.filter.LutApplier;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.RankFilters;
import ij.process.ColorProcessor;
import loci.formats.FormatException;
import loci.formats.ImageReader;
import loci.formats.gui.BufferedImageReader;
import me.tongfei.progressbar.ProgressBar;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.GregorianCalendar;
import java.util.List;

public class Appendices_ {
    private static int serie = 1;
    private static int crop = 1024;
    private static int overlap = 0;
    private static final String[] saveModes = new String[]{"None","Muscle","All"};
    private static String saveMode = saveModes[0];
    private static FileWriter log;
    private static String inputPath;
    private static String outputPath;
    private static String finalTable;
    private static String tableLympho;
    private static String tableArea;
    private static CSVWriter writerLympho;
    private static CSVWriter writerArea;

    public static int count_muscle(String path) {
        IJ.open(path);
        ImagePlus imp1 = IJ.getImage();
        ImagePlus dup = imp1.duplicate();
        dup.getProcessor().setMinAndMax(100,103);
        LutApplier applier = new LutApplier();
        applier.setup(null,dup);
        applier.run(dup.getProcessor());
        Image magenta = ReplaceRedWithMagenta((ColorProcessor) dup.getProcessor());
        dup.setImage(magenta);
        int blue = 3;
        ImageStack blueChannel = ChannelSplitter.getChannel(dup,blue);
        ImagePlus impBlue = new ImagePlus("4",blueChannel.getProcessor(1));
        RankFilters rankFilters = new RankFilters();
        rankFilters.rank(impBlue.getProcessor(),4.0,RankFilters.MEDIAN);
        IJ.saveAs(impBlue,"Tiff",outputPath + File.separator + "blue.tif");
        IJ.open(outputPath + File.separator + "blue.tif");
        Thresholder thresholder = new Thresholder();
        thresholder.run("mask");
        int options = ParticleAnalyzer.INCLUDE_HOLES; //bitwise options 1024 + 256+ 64. see flags constant in particleAnalyzer
        ResultsTable muscle_res = new ResultsTable();
        ParticleAnalyzer particleAnalyzer = new ParticleAnalyzer(options, Analyzer.getMeasurements(), muscle_res, 100, 3000, 0.1, 0.5);
        boolean success = particleAnalyzer.analyze(impBlue);
        IJ.exit();
        imp1.close();
        dup.close();
        impBlue.close();
        delete(outputPath + File.separator + "blue.tif");
        if (muscle_res.size() == 0)
            return 0;
        if (success) {
            double[] round = muscle_res.getColumnAsDoubles(muscle_res.getColumnIndex("Round"));
            double[] angle = muscle_res.getColumnAsDoubles(muscle_res.getColumnIndex("Angle"));
            double[] circ = muscle_res.getColumnAsDoubles(muscle_res.getColumnIndex("Circ."));
            double[] feretD = muscle_res.getColumnAsDoubles(muscle_res.getColumnIndex("Feret"));
            double[] minferet = muscle_res.getColumnAsDoubles(muscle_res.getColumnIndex("MinFeret"));
            int num = 0;
            for (int i = 0; i < round.length; i++) {
                try{
                    double feret = feretD[i] / minferet[i];
                    if (round[i] < 0.96 && angle[i] < 201 && circ[i] < 0.61 && feret > 1.45)
                        num++;
                }
                catch (NullPointerException e){
                    //pass
                }
            }
            return num;
        }
        else
            return -1;
    }
    private static void delete(String path){
        try {
            if (!deleteTemp(path, 0))
                log.write("Couldn`t delete file " + outputPath + File.separator + "blue.tif" + " of the " + path + "\n");
        }
        catch (IOException e) {
            //pass
        }
    }
    private static boolean deleteTemp(String path,int attempt){
        try {
            if (attempt > 5)
                return false;
            Files.deleteIfExists(Paths.get(path));
            return true;
        } catch (IOException e) {
            try {
                Thread.sleep(2000);
                attempt++;
                return deleteTemp(path, attempt);
            } catch (InterruptedException ex) {
                ex.printStackTrace();
                return false;
            }
        }
    }
    public static int count_lymphocytes(String path){
        IJ.open(path);
        ImagePlus imp1 = IJ.getImage();
        ImagePlus dup = imp1.duplicate();
        dup.getProcessor().setMinAndMax(53,70);
        LutApplier applier = new LutApplier();
        applier.setup(null,dup);
        applier.run(dup.getProcessor());
        Image magenta = ReplaceRedWithMagenta((ColorProcessor) dup.getProcessor());
        dup.setImage(magenta);
        IJ.saveAs(dup,"Tiff",outputPath + File.separator + "magenta.tif");
        IJ.open(outputPath + File.separator + "magenta.tif");//todo - change name for identification region
        Converter converter = new Converter();
        converter.run("8-bit");
        dup = IJ.getImage();
        RankFilters rankFilters = new RankFilters();
        rankFilters.rank(dup.getProcessor(),4.0,RankFilters.MEDIAN);
        Thresholder thresholder = new Thresholder();
        thresholder.run("");
        int options = ParticleAnalyzer.INCLUDE_HOLES;
        ResultsTable lympho_res = new ResultsTable();
        ParticleAnalyzer particleAnalyzer = new ParticleAnalyzer(options, Analyzer.getMeasurements(), lympho_res, 250, 1200, 0.5, 1);
        boolean success = particleAnalyzer.analyze(IJ.getImage());
        if (success) {
            try {
                lympho_res.saveAs(outputPath + File.separator + "lympho.csv");
                CSVReader reader = new CSVReader(new FileReader(outputPath + File.separator + "lympho.csv"));
                List<String[]> lines = reader.readAll();
                if (lines.size() == 0)
                    return -2;
                for (int i = 1; i < lines.size(); i++) {
                    writerLympho.writeNext(lines.get(i));
                }
                writerLympho.flush();
                reader.close();
                delete(outputPath + File.separator + "magenta.tif");
                delete(outputPath + File.separator + "lympho.csv");
                return lines.size() - 1;//minus header
            } catch (IOException e) {
                e.printStackTrace();
                return -3;
            }
        }
        else
            return -1;
    }
    public static double count_area(String path){
        IJ.open(path);
        Converter converter = new Converter();
        converter.run("8-bit");
        Thresholder thresholder = new Thresholder();
        thresholder.run("mask");
        int options = ParticleAnalyzer.INCLUDE_HOLES;
        ResultsTable area_res = new ResultsTable();
        ParticleAnalyzer particleAnalyzer = new ParticleAnalyzer(options, Analyzer.getMeasurements(), area_res, 0.0D , 1.0D / 0.0);
        boolean success = particleAnalyzer.analyze(IJ.getImage());
        if (success) {
            try {
                area_res.saveAs(outputPath + File.separator + "area.csv");
                CSVReader reader = new CSVReader(new FileReader(outputPath + File.separator + "area.csv"));
                List<String[]> lines = reader.readAll();
                int area_column = 0;
                double total_area = 0;
                if (lines.size() == 0)
                    return -2;
                for (int i = 0; i < lines.get(0).length; i++) {
                    if (lines.get(0)[i].contentEquals("Area")) {
                        area_column = i;
                        break;
                    }
                }
                for (int i = 1; i < lines.size(); i++) {
                    total_area += Double.parseDouble(lines.get(i)[area_column]);
                    writerArea.writeNext(lines.get(i));
                }
                writerArea.flush();
                reader.close();
                delete(outputPath + File.separator + "area.csv");
                return total_area;//minus header
            } catch (IOException e) {
                e.printStackTrace();
                return -3;
            }
        }
        else
            return -1;
    }
    public static Image ReplaceRedWithMagenta(ColorProcessor ip) {
        int w = ip.getWidth(), h = ip.getHeight();
        int[] pixels = (int[]) ip.getPixels();
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                int value = pixels[i + j * w];
                int red = (value >> 16) & 0xff;
                int green = (value >> 8) & 0xff;
                int blue = value & 0xff;
                if (false && blue > 16)
                    continue;
                pixels[i + j * w] = (red << 16) | (green << 8) | red;
            }
        }
        ip.setPixels(pixels);
        return ip.createImage();
    }

    /**
     * CMD arguments
     * 1 - -i --input input path
     * 2 - -o --output output path
     * 3 - -s --serie serie default 1
     * 4 - -c --crop crop  default 1024
     * 5 - -r --overlay overlay default 0
     * 6 - -h --help help
     * 7 - -m --saveMode [None,Muscle,All]
     * */
    private static void printHelp(){
        System.out.println("Appendices is program based on the ImageJ. version 1.0.\n" +
                "Authors: Anton S. Smirnov, Maria P. Myshkina, Natalya A. Ermakova\t Department of Bioinformatics, PRNRMU\n"+
                "Appendices is modificated SlideJ plugin for analysis whole-slide tiff snapshots of appendices.\n"+
                "It count number of lymphocytes in muscular tunic. It needs for automatic diagnostic of appendicitis.\n"+
                "Options:\n"+
                "-i --input = path to input directory, where whole-slide snapshots is keeped. Folder should contain only tiff snapshots\n"+
                "-o --output = path to output directory. It will contain log file, 3 tables with statistics and (option) tiles\n" +
                "-s --serie = number of layers, default equal to 1(image is 2D)\n"+
                "-c --crop = width of square tile, default equal to1024\n"+
                "-r --overlap =  overlap size between tiles, default equal to 0\n"+
                "-h --help = print help and exit\n"+
                "-m --saveMode = one of this variants, default equal to None\n"+
                "\tNone - don`t save tiles\n"+
                "\tMuscle - save tiles,where muscular tunic is deteted\n"+
                "\tAll - save all tiles\n");
        System.exit(0);
    }
    private static void parseArgs(String[] args){
        if (args.length == 0)
            printHelp();
        else {
            for (int i = 0; i < args.length; i += 2) {
                if (args[i].contentEquals("-h") || args[i].contentEquals("--help"))
                    printHelp();
                else if (args[i].contentEquals("-i") || args[i].contentEquals("--input"))
                    inputPath = args[i + 1];
                else if (args[i].contentEquals("-o") || args[i].contentEquals("--output"))
                    outputPath = args[i + 1];
                else if (args[i].contentEquals("-s") || args[i].contentEquals("--serie")) {
                    try {
                        serie = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        System.out.println("Incorrect parameter \"serie\". Can`t parse it from:" + args[i + 1]);
                        System.exit(1);
                    }
                }
                else if (args[i].contentEquals("-c") || args[i].contentEquals("--crop")) {
                    try {
                        crop = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        System.out.println("Incorrect parameter \"crop\". Can`t parse it from:" + args[i + 1]);
                        System.exit(1);
                    }
                }
                else if (args[i].contentEquals("-r") || args[i].contentEquals("--overlap")) {
                    try {
                        overlap = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        System.out.println("Incorrect parameter \"overlap\". Can`t parse it from:" + args[i + 1]);
                        System.exit(1);
                    }
                }
                else if (args[i].contentEquals("-m") || args[i].contentEquals("--saveMode")) {
                    boolean f = false;
                    for (int j = 0; j < saveModes.length; j++) {
                        if (args[i + 1].contentEquals(saveModes[j])) {
                            saveMode = saveModes[j];
                            f = true;
                        }
                    }
                    if (!f) {
                        System.out.println("Incorrect parameter \"saveMode\". Possible values = " + saveModes);
                        System.exit(1);
                    }
                }
                else {
                    System.out.println("Incorrect parameter " + args[i]);
                    System.exit(1);
                }
            }
        }
    }
    private static void setupStatics(){
        Prefs.blackBackground = false;
        Thresholder.setMethod("Default");
        Thresholder.setBackground("Light");
        Analyzer.precision = 3;
        //Analyzer.setMeasurements(945889);
        Analyzer.setMeasurements(Measurements.AREA + Measurements.CENTROID + Measurements.CENTER_OF_MASS + Measurements.PERIMETER +
                Measurements.RECT + Measurements.ELLIPSE + Measurements.SHAPE_DESCRIPTORS + Measurements.FERET + Measurements.SKEWNESS +
                Measurements.KURTOSIS + Measurements.AREA_FRACTION);
        Analyzer.setRedirectImage(null);
    }
    //TODO - несколько файлов
    public static void main(String[] args) throws Exception {
        parseArgs(args);
        setupStatics();
        File[] inputFiles = new File(inputPath).listFiles();
        int num_files = inputFiles.length;

        for (int i = 0; i < num_files; i++) {
            ImageReader imageReader = new ImageReader(); //Byte reader
            String inputName = inputFiles[i].getName();
            String imgPath = inputPath + File.separator + inputName;
            finalTable = outputPath + File.separator + inputName + "_finalTable.csv";
            tableLympho = outputPath + File.separator + inputName + "_tableLympho.csv";
            tableArea = outputPath + File.separator + inputName + "_tableArea.csv";
            log = new FileWriter(outputPath + File.separator + inputName + "_MeowLog.txt", true);
            GregorianCalendar now = new GregorianCalendar();
            SimpleDateFormat dateformat = new SimpleDateFormat("dd/MM/yyyy - HH:mm:ss");

            log.write("Start Meow Plugin: " + dateformat.format(now.getTime()) + "\n\r");
            log.flush();

            CSVWriter writer = new CSVWriter(new FileWriter(finalTable, true), ';');
            writer.writeNext(new String[]{"ID", "Lymphocytes", "Area", "UnitLymphocytes"});

            writerLympho = new CSVWriter(new FileWriter(tableLympho, true), ';');
            writerArea = new CSVWriter(new FileWriter(tableArea, true), ';');

            long startTime = System.currentTimeMillis();//Time recorder

            log.write(dateformat.format(startTime) + "\t" + imgPath + "\n\r");
            try {
                imageReader.setId(imgPath);
                BufferedImageReader buffImageReaderTest = new BufferedImageReader();
                buffImageReaderTest.setId(imgPath);
                int s = (serie - 1);
                buffImageReaderTest.setSeries(s);
                int buffYTest = buffImageReaderTest.getSizeY();
                int buffXTest = buffImageReaderTest.getSizeX();
                int tilex = crop;
                int x_coor = 0;
                int byX = (int) Math.ceil((double) buffXTest / crop);
                int byY = (int) Math.ceil((double) buffYTest / crop);
                int total = byX * byY;
                String taskName = "Processing " + inputName;
                ProgressBar pb = new ProgressBar(taskName, total);
                pb.setExtraMessage((i+1) + "files from " + num_files);
                for (int j = 0; j == 0; x_coor = x_coor + crop - overlap) {
                    int tiley = crop;
                    int y_coor = 0;
                    for (int z = 0; z == 0; y_coor = y_coor + crop - overlap) {
                        if ((x_coor + crop) >= buffXTest) {
                            tilex = buffXTest - x_coor;
                            j = 1;
                        }

                        if ((y_coor + crop) >= buffYTest) {
                            tiley = buffYTest - y_coor;
                            z = 1;
                        }
                        log.write("Processing... X = " + x_coor + " from " + buffXTest + " Y = " + y_coor + " from " + buffYTest + "\n\r");

                        BufferedImage rgbImage = buffImageReaderTest.openImage(0, x_coor, y_coor, tilex, tiley);
                        /* Tiles are stored in TIFF format with a file name that reflects their position
                           on the overall digital slide according to the following template:
                           <OriginalFileName.ext>__<series>_<Xorigin>_<Yorigin>.tif */
                        String title = inputName + "__" + serie + "_" + x_coor + "_" + y_coor + ".tif";
                        String path = outputPath + File.separator + title;
                        ImagePlus imp = new ImagePlus(title, rgbImage);
                        IJ.run(imp, "RGB Color", "");
                        IJ.saveAs(imp, "Tiff", path);
                        int num_muscle_cell = count_muscle(path);
                        if (num_muscle_cell == -1) {
                            log.write("Particle Analyzer couldnt analyse and find muscle cells in " + title + "\n");
                        } else if (num_muscle_cell > 25) {
                            int num_lympho_cell = count_lymphocytes(path);
                            switch (num_lympho_cell) {
                                case -1:
                                    log.write("Particle Analyzer couldnt analyse and find lymphocytes in " + title + "\n");
                                    break;
                                case -2:
                                    log.write("Don`t find result table " + title + "\n");
                                    break;
                                case -3:
                                    log.write("IOException when open result table " + title + "\n");
                                    break;
                                default:
                                    double area = count_area(path);
                                    if (area == -1) {
                                        log.write("Particle Analyzer couldnt analyse and calculate area in " + title + "\n");
                                    } else if (area == -2) {
                                        log.write("Don`t find area result table " + title + "\n");
                                    } else if (area == -3) {
                                        log.write("IOException when open area result table " + title + "\n");
                                    } else {
                                        writer.writeNext(new String[]{title, Integer.toString(num_lympho_cell), Double.toString(area), Double.toString((double) num_lympho_cell / area)});
                                        log.write("title = " + title + " has " + num_muscle_cell + " muscle cells and " + num_lympho_cell + " lymphocytes.\n");
                                    }
                                    break;
                            }
                        } else {
                            log.write("title = " + title + " has not detected enough muscle cells.\n\r");
                        }
                        imp.close();
                        if (saveMode.contentEquals(saveModes[0]) || (saveMode.contentEquals(saveModes[1]) && num_muscle_cell < 25)) {
                            delete(path);
                        }
                        pb.step();
                        log.flush();
                    }
                }
                long endTime = System.currentTimeMillis();
                long totalTime = endTime - startTime;
                pb.close();
                buffImageReaderTest.close();
                log.write("File " + inputName + ", series number " + serie + ", tile size " + crop + ", overlap " + overlap + ", time " + totalTime + " ms. " + "\n\r");
                log.flush();
            }
            catch (FormatException exc) {
                System.out.println("Format exception" + exc.getMessage());
                exc.printStackTrace();
            }
            catch (IOException exc) {
                System.out.println("I/O exception" + exc.getMessage());
                exc.printStackTrace();
            }
            finally {
                log.flush();
                log.close();
                writer.flush();
                writer.close();
                writerArea.flush();
                writerArea.close();
                writerLympho.flush();
                writerLympho.close();
                imageReader.close();
            }
        }
    }
}