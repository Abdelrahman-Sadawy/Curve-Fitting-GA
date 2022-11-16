import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;

public class Main {
    static int degree;
    final static double lower = -10.0;
    final static double upper = 10.0;
    static final double pc = 0.5;

    static final double pm = 0.01;
    static final int generationNum = 100;
    static final int b = 1;
    static final int T = generationNum - 1;

    static double bestValue;

    public static ArrayList<ArrayList<Double>> populationInit(int degree) {
        Random random = new Random();
        int populationSize = random.nextInt(500 - 100) + 100;
        ArrayList<ArrayList<Double>> chromosomes = new ArrayList<>();
        for (int i = 0; i < populationSize; i++) {
            chromosomes.add(new ArrayList<>());
            for (int j = 0; j <= degree; j++) {
                chromosomes.get(i).add(random.nextDouble(upper - lower) + lower);
            }
        }
        return chromosomes;
    }

    public static ArrayList<Double> fitnessCalc(
            double[] x, double[] y, ArrayList<ArrayList<Double>> chromosomes, int pointsNum, int degree) {
        ArrayList<Double> fitness = new ArrayList<>();

        for (ArrayList<Double> chromosome : chromosomes) {
            double totalError = 0.0;

            for (int i = 0; i < pointsNum; i++) {
                double error = 0.0;
                for (int j = 0; j <= degree; j++) {
                    error += chromosome.get(j) * Math.pow(x[i], j);
                }
                error -= y[i];
                error = Math.pow(error, 2);
                totalError += error;
            }
            totalError /= pointsNum;
            fitness.add(1.0 / totalError);
        }
        return fitness;
    }

    public static int bestFitness(ArrayList<Double> fitness, int ind1, int ind2) {
        if (fitness.get(ind1) > fitness.get(ind2))
            return ind1;

        return ind2;
    }

    public static ArrayList<Integer> tournamentSelection(ArrayList<Double> fitness, int chromosomesSize) {
        Random random = new Random();

        int selectionNum;
        ArrayList<Integer> parentsInd = new ArrayList<>();

        for (int i = 0; i < chromosomesSize; i++) {
            int rand1 = random.nextInt(chromosomesSize);
            int rand2 = random.nextInt(chromosomesSize);

            selectionNum = bestFitness(fitness, rand1, rand2);
            parentsInd.add(selectionNum);
        }
        return parentsInd;
    }

    public static ArrayList<ArrayList<Double>> crossover(ArrayList<Integer> parentsInd, ArrayList<ArrayList<Double>> chromosomes) {
        Random random = new Random();
        ArrayList<ArrayList<Double>> offSprings = new ArrayList<>();
        int size = chromosomes.size();

        for (int i = 0; i < chromosomes.size(); i++) {
            offSprings.add(new ArrayList<>());
        }

        if (chromosomes.size() % 2 != 0) {

            for (int j = 0; j < chromosomes.get(1).size(); j++) {
                offSprings.get(chromosomes.size() - 1).add(chromosomes.get(parentsInd.get(chromosomes.size() - 1)).get(j));
            }

            size = chromosomes.size() - 1;
        }

        for (int i = 0; i < size; i += 2) {
            double randomProbability = random.nextDouble(1);

            int r1 = random.nextInt((chromosomes.get(1).size() - 1) - 1) + 1;
            int r2 = random.nextInt((chromosomes.get(1).size() - 1) - 1) + 1;

            //System.out.println("RP " + randomProbability + " R1 " + r1 + " R2 " + r2);

            int crossoverPoint1, crossoverPoint2;

            if (r1 < r2) {
                crossoverPoint1 = r1;
                crossoverPoint2 = r2;
            } else {
                crossoverPoint1 = r2;
                crossoverPoint2 = r1;
            }

            if (randomProbability <= pc) {
                for (int j = 0; j < crossoverPoint1; j++) {
                    offSprings.get(i).add(chromosomes.get(parentsInd.get(i)).get(j));
                    offSprings.get(i + 1).add(chromosomes.get(parentsInd.get(i + 1)).get(j));
                }
                for (int j = crossoverPoint1; j < crossoverPoint2; j++) {
                    offSprings.get(i).add(chromosomes.get(parentsInd.get(i + 1)).get(j));
                    offSprings.get(i + 1).add(chromosomes.get(parentsInd.get(i)).get(j));
                }

                for (int k = crossoverPoint2; k < chromosomes.get(1).size(); k++) {
                    offSprings.get(i).add(chromosomes.get(parentsInd.get(i)).get(k));
                    offSprings.get(i + 1).add(chromosomes.get(parentsInd.get(i + 1)).get(k));

                }
            } else {
                for (int j = 0; j < chromosomes.get(1).size(); j++) {
                    offSprings.get(i).add(chromosomes.get(parentsInd.get(i)).get(j));
                    offSprings.get(i + 1).add(chromosomes.get(parentsInd.get(i + 1)).get(j));
                }
            }
        }


        return offSprings;
    }

    public static ArrayList<ArrayList<Double>> mutation(ArrayList<ArrayList<Double>> offSprings, int t) {
        Random random = new Random();

        double deltaLower, deltaUpper;
        double y;
        double deltaTY;
        double pow1, pow2;


        for (ArrayList<Double> offSpring : offSprings) {
            for (int j = 0; j < offSprings.get(1).size(); j++) {
                double probabilityOfMutation = random.nextDouble(1);
                deltaLower = offSpring.get(j) - lower;
                deltaUpper = upper - offSpring.get(j);

                float r1 = random.nextFloat(1);
                float r = random.nextFloat(1);

                if (probabilityOfMutation <= pm) {
                    if (r1 <= 0.5)
                        y = deltaLower;
                    else
                        y = deltaUpper;

                    pow1 = Math.pow((1.0 - (1.0 * t / T)), b);
                    pow2 = Math.pow(r, pow1);
                    deltaTY = y * (1.0 - pow2);

                    if (y == deltaLower) {
                        offSpring.set(j, offSpring.get(j) - deltaTY);
                    } else {
                        offSpring.set(j, offSpring.get(j) + deltaTY);
                    }


                }
            }
        }
        return offSprings;
    }

    public static ArrayList<ArrayList<Double>> replacement(ArrayList<ArrayList<Double>> oldGeneration,
                                                           ArrayList<ArrayList<Double>> offSprings, ArrayList<Double> parentFitness, ArrayList<Double> offSpringsFitness) {

        double parentMax, offSpringMax;
        ArrayList<ArrayList<Double>> newGeneration = new ArrayList<>();
        for (int i = 0; i < oldGeneration.size(); i++) {
            //newGeneration.add(new ArrayList<>());

            parentMax = Collections.max(parentFitness);
            offSpringMax = Collections.max(offSpringsFitness);


            if (parentMax > offSpringMax) {
                newGeneration.add(oldGeneration.get(parentFitness.indexOf(parentMax)));
                parentFitness.set(parentFitness.indexOf(parentMax), -Double.MAX_VALUE);
            } else {
                newGeneration.add(offSprings.get(offSpringsFitness.indexOf(offSpringMax)));
                offSpringsFitness.set(offSpringsFitness.indexOf(offSpringMax), -Double.MAX_VALUE);
            }

        }
        return newGeneration;

    }

    public static void selectBestChromosome(ArrayList<Double> fitness, ArrayList<ArrayList<Double>> chromosomes, int dataSetIndex) throws IOException, IOException {
        ArrayList<Double> bestIndividual;
        bestValue = Collections.max(fitness);
        bestIndividual = chromosomes.get(fitness.indexOf(bestValue));
        double meanSquareError = 1.0 / bestValue;
        FileWriter fileWriter = new FileWriter("output.txt", true);
        fileWriter.write("Dataset Index: " + dataSetIndex + "\nCoefficients of the polynomial function: \n" + bestIndividual + "\nMean square error: " + meanSquareError + "\n");
        fileWriter.close();

    }
    public static void main(String[] args) throws IOException {
        File file = new File("curve_fitting_input.txt");
        Scanner sc = new Scanner(file);
        int testCaseNum = sc.nextInt();
        int pointsNum;

        double[] x;
        double[] y;

        ArrayList<ArrayList<Double>> chromosomes;
        ArrayList<ArrayList<Double>> offSprings;

        ArrayList<Double> fitness;
        ArrayList<Integer> parentsInd;
        ArrayList<Double> offSpringFitness;


        for (int i = 0; i < testCaseNum; i++) {
            pointsNum = sc.nextInt();
            degree = sc.nextInt();
            x = new double[pointsNum];
            y = new double[pointsNum];

            for (int j = 0; j < pointsNum; j++) {
                x[j] = sc.nextDouble();
                y[j] = sc.nextDouble();
            }

            chromosomes = populationInit(degree);
            fitness = fitnessCalc(x, y, chromosomes, pointsNum, degree);

            for(int k = 0; k < generationNum - 1; k++){
               parentsInd = tournamentSelection(fitness, chromosomes.size());

                offSprings = crossover(parentsInd, chromosomes);
                mutation(offSprings, k);
                offSpringFitness = fitnessCalc(x, y, offSprings, pointsNum, degree);

                chromosomes = replacement(chromosomes, offSprings, fitness, offSpringFitness);

               fitness = fitnessCalc(x, y, chromosomes, pointsNum, degree);

            }
            selectBestChromosome(fitness, chromosomes, i);
        }

    }
}