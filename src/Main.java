import java.io.File;
import java.io.FileNotFoundException;
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

    public static ArrayList<ArrayList<Double>> populationInit(int degree) {
        Random random = new Random();
        int populationSize = random.nextInt(50 - 5) + 5;
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
            double totalError = 0;

            for (int i = 0; i < pointsNum; i++) {
                double error = 0;
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
        for (int i = 0; i < parentFitness.size(); i++) {
            //newGeneration.add(new ArrayList<>());

            parentMax = Collections.max(parentFitness);
            offSpringMax = Collections.max(offSpringsFitness);


            if (parentMax > offSpringMax) {
                newGeneration.add(oldGeneration.get(parentFitness.indexOf(parentMax)));
                parentFitness.remove(parentMax);
            } else {
                newGeneration.add(offSprings.get(offSpringsFitness.indexOf(offSpringMax)));
                offSpringsFitness.remove(offSpringMax);
            }

        }
        return newGeneration;

    }


    public static void main(String[] args) throws FileNotFoundException {
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

        /* chromosomes.add(new ArrayList<>());
        chromosomes.get(0).add(1.95);
        chromosomes.get(0).add(8.16);
        chromosomes.get(0).add(-2.0);

        chromosomes.add(new ArrayList<>());
        chromosomes.get(1).add(4.26);
        chromosomes.get(1).add(-7.4);
        chromosomes.get(1).add(-2.5);

        double[] x1 = new double[]{1.0, 2.0, 3.0, 4.0};
        double[] y2 = new double[]{5.0, 8.0, 13.0, 20.0};

        System.out.println("Fitness ");
        System.out.println(fitness = fitnessCalc(x1, y2, chromosomes, 4, 2));
        System.out.println("parentInd ");
        System.out.println(tournamentSelection(fitness, chromosomes.size()));
        System.out.println("Mutation ");
        System.out.println(mutation(chromosomes, 97));*/

        for (int i = 0; i < testCaseNum; i++) {
            pointsNum = sc.nextInt();
            degree = sc.nextInt();
            x = new double[pointsNum];
            y = new double[pointsNum];

            for (int j = 0; j < pointsNum; j++) {
                x[j] = sc.nextDouble();
                y[j] = sc.nextDouble();
            }

            System.out.println("Population");
            System.out.println(chromosomes = populationInit(degree));
            System.out.println("fitness ");
            System.out.println(fitness = fitnessCalc(x, y, chromosomes, pointsNum, degree));
            System.out.println("Parents ");
            System.out.println(parentsInd = tournamentSelection(fitness, chromosomes.size()));

            offSprings = crossover(parentsInd, chromosomes);
            System.out.println("Mutation");
            System.out.println(mutation(offSprings, 97));
            offSpringFitness = fitnessCalc(x, y, offSprings, pointsNum, degree);

            System.out.println("Replacement: ");
            System.out.println(replacement(chromosomes, offSprings, fitness, offSpringFitness));

        }

    }
}