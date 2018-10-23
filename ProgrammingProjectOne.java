package programmingprojectone;

// Author: Abraham Yepremian
// Class: CS 3010.01

import java.io.File;
import java.util.LinkedList;
import java.util.Scanner;

public class ProgrammingProjectOne {
    
    private static Scanner scanner;
    private static double matrix[][];
    // constant final value boolean to automatically run the extra credit or not
    // this is useful when I want to see how equation size impacts execution time
    private final static boolean AUTO_RUN_EXTRA_CREDIT = false;
    // constant final value for number of equations to generate for extra credit
    private final static int N_FOR_EXTRA_CREDIT = 25;
    // algorithm to run for extra credit if run immediately
    // 0 will be Gauss Elimination with Scaled Partial Pivoting, 1 will be Jacobi, and 2 will be Gauss-Siedel
    private final static int ALGORITHM_TO_RUN = 0;
    // default desired error for Jacobi or Gauss-Siedel. This is for extra credit when run immediately.
    // set to half a percent
    private final static double EC_DESIRED_ERROR = 0.005;
        
    private static void printMatrix(double[][] matrix) {
        for(int i =0; i <matrix.length ; i++) {
            System.out.print("| ");
            for(int j= 0 ; j<matrix[i].length ; j++) {
                System.out.print(matrix[i][j] +  " | ");
            }
            System.out.println();
        }
    }
    
    private static void launchMenu() {
        System.out.println("\nInteractive Menu");
        System.out.println("1. Print Matrix");
        System.out.println("2. The Scaled Partial Pivoting Method for Gaussian Elimination.");
        System.out.println("3. Jacobi Iterative Method");
        System.out.println("4. Gauss-Seidel Method");
        System.out.println("5. Exit");
        
        System.out.print("Choose option (1-5): ");
        int userChoice = scanner.nextInt();
        
        switch(userChoice) {
            case 1: 
                System.out.println("\nHere is your matrix.");
                printMatrix(matrix);
                launchMenu();
                break;
            case 2:
                System.out.println("\nGaussian Elimination with Scaled Partial Pivoting");
                printGaussAnswers(runGaussianElimination(matrix));
                launchMenu();
                break;
            case 3:
                System.out.println("Jacobi Iterative Method");
                System.out.println("A smaller desired error will produce more iterations.");
                System.out.print("Input your desired error (between 0 and 1): ");
                double desiredError = scanner.nextFloat();
                runJacobi(matrix, desiredError);
                launchMenu();
                break;
            case 4:
                System.out.println("Gauss Siedel Method");
                System.out.println("A smaller desired error will produce more iterations.");
                System.out.print("Input your desired error (between 0 and 1): ");
                double desiredErrorSeidel = scanner.nextFloat();
                runGaussSiedel(matrix, desiredErrorSeidel);
                launchMenu();
                break;
            case 5:
                System.out.println("Thank you for using my program.");
                System.exit(0);
                break;
            default:
                System.out.println("\nChoice not recognized, please try again.");
                launchMenu();
                break;
        }
    }
    
    private static double[][] pivot(double[][] matrix, int p, int q) {   
        for (int i = 0; i < matrix.length; i++) {
            double x = matrix[i][q] / matrix[p][q];
            for (int j = 0; j <= (matrix.length ); j++) {
                if (i != p && j != q) {
                        matrix[i][j] -= x * matrix[p][j];
                }
            }
        }

        for (int i = 0; i < matrix.length; i++) {
            if (i != p) {
                matrix[i][q] = 0.0;
            }
        }

        for (int j = 0; j <= (matrix.length); j++) {
            if (j != q) {
                matrix[p][j] /= matrix[p][q];
            }
        }

        matrix[p][q] = 1.0;
        return matrix;
    }
		
    private static double[][] runGaussianElimination(double[][] matrix) {
        // find maximum
        double[][] maxMatrix = new double[matrix.length][1]; 
        for(int i =0; i < matrix.length; i++) {
            maxMatrix[i][0] = Math.abs(matrix[i][0]); 
            for(int j= 0; j < matrix[i].length - 1 ; j++) {
                    if (Math.abs(matrix[i][j]) > maxMatrix[i][0]) {
                            maxMatrix[i][0] =  Math.abs(matrix[i][j]); 
                    }
            }
        }

        double pivotMaximum = -1000; 
        int rowSwapIndex = 0 ;
        boolean swapRow = false;  
        double[][] pivotMatrix = new double[matrix.length][1]; 
        for(int i =0; i < matrix.length; i++) {
            for(int j = 0; j < matrix.length; j++ ) {
                    pivotMatrix[j][0] = matrix[j][i] / maxMatrix[j][0];  
                    if (pivotMatrix[j][0] > pivotMaximum) {
                            pivotMaximum = pivotMatrix[j][0]; 
                            rowSwapIndex = j; 
                            swapRow = true; 
                    }
            }
            
            if(swapRow) {
                // swap row
                double[] tempValue = matrix[i];
                matrix[i] = matrix[rowSwapIndex];
                matrix[i] = tempValue; 
            }
            
            swapRow = false;
            // calls pivot method on matrix
            matrix = pivot(matrix, i, i);
            System.out.println("Intermediary Matrix");
            printMatrix(matrix);
        }
        return matrix;
    }
    
    private static void printGaussAnswers(double[][] finalMatrix) {
        System.out.println("\nHere are the final values from the matrix");
        for(int i = 0; i<finalMatrix.length; i++) {
            System.out.println("x" + (i+1) + "=" + finalMatrix[i][finalMatrix[i].length - 1]);
        }
    }
    
    private static void generateExtraCreditMatrix() {
        int n = N_FOR_EXTRA_CREDIT; 
        int columns = n + 1 ;
        matrix = new double[n][columns]; 
        for(int i =0; i < n ; i++) {
            for(int j= 0 ; j<matrix[i].length ; j++) {
                matrix[i][j] = (int )(Math.random() * 50 + 1);
            }
        }
    }
    
    private static void runJacobi(double[][] matrix, double desiredError) {
        int kmax = 50;
        double[] xArray = new double[matrix.length];
        double[] yArray = new double[matrix.length];
        int i, j, k, n = matrix.length, h = 0;
        double diag, sum;
        for (k = 0; k < kmax; k++) {
            for (i = 0; i < n; i++) {
                sum = matrix[i][n];
                diag = matrix[i][i];
                for (j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= matrix[i][j] * yArray[j];
                    }
                }
                xArray[i] = sum/diag;
            }

            System.out.println("\nIteration " + (k + 1));
            for (int g = 0; g < xArray.length; g++) {
                System.out.println("x" + (g+1) + "=" + yArray[g] + " ");
            }

            for (int u = 0; u < matrix.length; u++) {
                if (Math.abs(xArray[u] - yArray[u]/xArray[u]) <= desiredError) {
                    System.out.println("Desired error has been reached");
                    return;
                }
            }
            
            for (int m = 0; m < matrix.length; m++) {
                yArray[m] = xArray[m];
            }
        }
        System.out.println("Maximum iterations have been reached.");
    }
    
    private static void printColumnVector(double[] columnVector) {
        System.out.print("X = {");
        for(int i = 0; i < columnVector.length; i++) {
            if(i == columnVector.length - 1 ) {
                System.out.print(columnVector[i] + "}\n");
            } else {
                System.out.print(columnVector[i] + ", ");
            }
        }
    }

    private static void runGaussSiedel(double[][] matrix, double desiredError) {
        int kmax = 50;
        double[] xArray = new double[matrix.length];
        double[] yArray = new double[matrix.length];
        int i, j, k, n = matrix.length, h = 0;
        double diag, sum;
        for (k = 0; k < kmax; k++) {
            for (i = 0; i < n; i++) {
                sum = matrix[i][n];
                diag = matrix[i][i];
                for (j = 0; j < n; j++) {
                    if (j != i) {
                        sum -= matrix[i][j] * xArray[j];
                    }
                }
                xArray[i] = sum/diag;
            }

            System.out.println("\nIteration " + (k + 1));
            for (int g = 0; g < xArray.length; g++) {
                System.out.println("x" + (g+1) + "=" + yArray[g] + " ");
            }

            for (int u = 0; u < matrix.length; u++) {
                if (Math.abs(xArray[u] - yArray[u]/xArray[u]) <= desiredError) {
                    System.out.println("Desired error has been reached.");
                    return;
                }
            }

            for (int m = 0; m < matrix.length; m++) {
                yArray[m] = xArray[m];
            }
        }
        System.out.println("Maximum iterations have been reached.");
    }

    public static void main(String[] args) {
        if(AUTO_RUN_EXTRA_CREDIT) {
            long startTime = System.currentTimeMillis();
            System.out.println("Run Extra Credit Immediately");
            System.out.println("Generating matrix of " + N_FOR_EXTRA_CREDIT + " equations...");
            generateExtraCreditMatrix();
            System.out.println("\nMatrix has been set.");
            printMatrix(matrix);
            
            switch(ALGORITHM_TO_RUN) {
                case 0:
                    System.out.println("\nGaussian Elimination with Scaled Partial Pivoting");
                    printGaussAnswers(runGaussianElimination(matrix));
                    break;
                case 1:
                    System.out.println("\nJacobi Method");
                    runJacobi(matrix, EC_DESIRED_ERROR);
                    break;
                case 2:
                    System.out.println("\nGauss-Siedel Method");
                    runGaussSiedel(matrix, EC_DESIRED_ERROR);
                    break;
                default:
                    System.out.println("\nAlgorithm to run is unrecognized");
                    break;
            }
            System.out.println("Time taken: " + (System.currentTimeMillis() - startTime) + " millisecond(s).");
        } else {
            System.out.println("Programming Project 1 CS 301");
            System.out.println("You will have the option to use 3 algorithms on a selected matrix.");
            System.out.println("Choose 1 to manually input your matrix with equations \nChoose 2 to have the program read a file for your matrix \nChoose 3 for a randomly generate extra credit matrix");
            System.out.print("Choose input (1-3): ");
            scanner = new Scanner(System.in);
            int userInput = scanner.nextInt();
            while(userInput != 1 && userInput != 2 && userInput !=3) {
                System.out.println("Invalid input, please try again.");
                System.out.print("Choose input (1-3): ");
                userInput = scanner.nextInt();
            }
            if(userInput == 1) {
                System.out.print("How many linear equations to solve? (n <= 10): ");
                int n = scanner.nextInt();
                int columns = n + 1;
                matrix = new double[n][columns]; 
                while(n > 10 || n < 1) {
                    System.out.print("Invalid number of equations, try again (n <= 10): ");
                    n = scanner.nextInt();
                    columns = n + 1;
                }
                for(int i = 0; i < n; i++) {
                    System.out.println("For equation #" + (i+1));
                    for(int j = 0; j < columns; j++) {
                        if(j<n) {
                            System.out.println("Enter the coffiecent for variable #" + (j+1) + ": ");
                            matrix[i][j] = scanner.nextFloat(); 
                        } else {
                            System.out.println("Enter the constant term: ");
                            matrix[i][j] = scanner.nextFloat(); 
                        }
                    } 
                }
                System.out.println("\nMatrix has been set.");
                printMatrix(matrix);
                launchMenu();
            } else if(userInput == 2) {
                System.out.println("You chose to import your matrix by file name (ex. input.txt)");
                System.out.print("Enter the full file name to import from: ");
                // initialize linked list for storing lines from file
                LinkedList<String> fileLines = new LinkedList<>();
                try {
                    // reset scanner from taking integer to being able to take a string
                    scanner.nextLine();
                    String fileName = scanner.nextLine();
                    System.out.println("Importing " + fileName);
                    Scanner fileScanner = new Scanner(new File(fileName));
                    while(fileScanner.hasNext()) {
                        fileLines.add(fileScanner.nextLine());
                    }
                    fileScanner.close();
                } catch(Exception e) {
                    System.out.println("There was an exception when trying to load the file.");
                    System.out.println("Exception was: " + e);
                }
                // file was handled at this point
                // now we will extract the matrix from the fileLines stored from the file
                int i =0; 
                String arrayLine[]; 
                int n = fileLines.size(); 
                int columns = fileLines.size() + 1 ;
                matrix = new double[n][columns]; 
                while(i < fileLines.size()) {
                    for(int a = 0 ; a < matrix.length; a++) {
                        arrayLine = fileLines.get(i).split(" ");
                        for(int b = 0; b < matrix[a].length; b++) {
                            matrix[a][b] = Double.parseDouble(arrayLine[b]); 
                        }
                        i++;
                    }
                }
                System.out.println("\nMatrix has been set.");
                printMatrix(matrix);
                launchMenu();
            } else if(userInput == 3) {
                System.out.println("Extra Credit Matrix");
                System.out.println("Generating matrix of " + N_FOR_EXTRA_CREDIT + " equations...");
                generateExtraCreditMatrix();
                System.out.println("\nMatrix has been set.");
                printMatrix(matrix);
                launchMenu();
            } else {
                System.out.println("Input error.");
            }
        }
    }
    
}
