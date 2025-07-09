import java.util.InputMismatchException;
import java.lang.IllegalArgumentException;
import java.util.Scanner;

public class AllInOneCalculator2 {

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

         while (true){
            System.out.println("Calculator Menu:");
            System.out.println("1. Arithmetic calculation");
            System.out.println("2. Logarithmic Equation Solver");
            System.out.println("3. Linear Equation Solver");
            System.out.println("4. Polynomial Equation Solver");
            System.out.println("5. Matrix Operations");
            System.out.println("6. Trigonometry Calculator");
            System.out.println("7. Exit");
            System.out.print("Enter your choice (1-7): ");

            int choice;
            try {
                choice = scanner.nextInt();
            } catch (InputMismatchException e) {
                System.out.println("InputMismatchException");
                scanner.nextLine();
                continue;
            }

            if (choice == 7) {
                System.out.println("Exiting the calculator. Goodbye!");
                break;
            }

            switch (choice) {
                case 1:
                    runArithmeticCalculator(scanner);
                    break;
                case 2:
                    runLogarithmicEquationSolver(scanner);
                    break;
                case 3:
                    runLinearEquationSolver(scanner);
                    break;
                case 4:
                    runPolynomialEquationSolver(scanner);
                    break;
                case 5:
                    runMatrixOperations(scanner);
                    break;
                case 6:
                    runTrigonometryCalculator(scanner);
                    break;
                default:
                    System.out.println("Invalid choice. Please enter a number between 1 and 7.");
            }
        }

        scanner.close();
    }

    private static void runArithmeticCalculator(Scanner scanner) {
        while (true) {
            System.out.println("Scientific Calculator Menu:");
            System.out.println("1. Addition");
            System.out.println("2. Subtraction");
            System.out.println("3. Multiplication");
            System.out.println("4. Division");
            System.out.println("5. Square Root");
            System.out.println("6. Exponentiation");
            System.out.println("7. Back to Main Menu");
            System.out.print("Enter your choice (1-7): ");

            int choice;
            try {
                choice = scanner.nextInt();
            } catch (InputMismatchException e) {
                System.out.println("InputMismatchException");
                scanner.nextLine(); // Consume the invalid input
                continue; // Restart the loop
            }

            if (choice == 7) {
                System.out.println("Returning to the main menu.");
                break;
            }

            double result=0.0;

            switch (choice) {
                case 1:
                    result+= addition(scanner);
                    break;
                case 2:
                    result+= subtraction(scanner);
                    break;
                case 3:
                    result+= multiplication(scanner);
                    break;
                case 4:
                    result+= division(scanner);
                    break;
                case 5:
                    result+= squareRoot(scanner);
                    break;
                case 6:
                    result+= exponentiation(scanner);
                    break;
                default:
                    System.out.println("Invalid choice. Please enter a number between 1 and 7.");
            }

            System.out.println("Result: " + result);
        }
        
    }

    private static double addition(Scanner scanner) {
        
        System.out.print("Enter the first number: ");
        double num1 = scanner.nextDouble();
        System.out.print("Enter the second number: ");
        double num2 = scanner.nextDouble();
        return num1 + num2;
    
    }

    private static double subtraction(Scanner scanner) {
        
        System.out.print("Enter the first number: ");
        double num1 = scanner.nextDouble();
        System.out.print("Enter the second number: ");
        double num2 = scanner.nextDouble();
        return num1 - num2;

    }

    private static double multiplication(Scanner scanner) {
        
        System.out.print("Enter the first number: ");
        double num1 = scanner.nextDouble();
        System.out.print("Enter the second number: ");
        double num2 = scanner.nextDouble();
        return num1 * num2;

    }

    private static double division(Scanner scanner) {
        System.out.print("Enter the numerator: ");
        double numerator = scanner.nextDouble();
        System.out.print("Enter the denominator: ");
        double denominator = scanner.nextDouble();

        if (denominator == 0) {
            System.out.println("Error: Cannot divide by zero.");
            return Double.NaN; // Not a Number
        }

        return numerator / denominator;
    }

    private static double squareRoot(Scanner scanner) {
        
        System.out.print("Enter a number: ");
        double num = scanner.nextDouble();

        if (num < 0) {
            System.out.println("Error: Cannot calculate square root of a negative number.");
            return Double.NaN;
        }

        return Math.sqrt(num);

    }

    private static double exponentiation(Scanner scanner) {
        
            System.out.print("Enter the base: ");
            double base = scanner.nextDouble();
            if (base < 0) {
                System.out.println("Error: Cannot calculate exponentiation with a negative base.");
                return Double.NaN;
            }
            System.out.print("Enter the exponent: ");
            double exponent = scanner.nextDouble();

            if (Double.isNaN(Math.pow(base, exponent))) {
                System.out.println("Error: Invalid input for exponentiation.");
                return Double.NaN;
            }
            return Math.pow(base, exponent);

    }
    

    private static void runLogarithmicEquationSolver(Scanner scanner) {

        System.out.println("Logarithmic Equation Solver");

        System.out.println("Enter the base of the logarithm:");
        double base = scanner.nextDouble();

        System.out.println("Enter the result of the logarithm:");
        double result = scanner.nextDouble();

        // Solve the logarithmic equation
        double exponent = solveLogarithmicEquation(base, result);

        // Display the solution
        System.out.println("Exponent: " + exponent);

    }

    private static double solveLogarithmicEquation(double base, double result) {
        // Handling edge case where the result is very close to zero
        if (Math.abs(result) < 1e-10 || base <= 0 || base == 1) {
            System.out.println("Result is very close to zero. Returning NaN.");
            return Double.NaN;
        }

        // Using logarithmic property: log_base(result) = exponent
        double exponent = Math.log(result) / Math.log(base);
        return exponent;
    }

    private static void runLinearEquationSolver(Scanner scanner) {
        
            System.out.println("Linear Equation Solver");
    
            System.out.println("Enter the number of variables (2 or 3):");
            int numVariables = scanner.nextInt();
    
            if (numVariables < 2) {
                System.out.println("Invalid number of variables. Exiting.");
                return;
            }
    
            System.out.println("Enter the coefficients and constants of the equations:");
    
            // Coefficient matrix
            double[][] coefficients = new double[numVariables][numVariables];
    
            // Constants vector
            double[] constants = new double[numVariables];
    
            for (int i = 0; i < numVariables; i++) {
                System.out.println("Equation " + (i + 1) + ":");
                for (int j = 0; j < numVariables; j++) {
                    System.out.println("Enter coefficient for variable x" + (j + 1) + ":");
                    coefficients[i][j] = scanner.nextDouble();
                }
                System.out.println("Enter constant term:");
                constants[i] = scanner.nextDouble();
            }
    
            // Solve the system of linear equations
            double[] solution = solveLinearEquations(coefficients, constants);
    
            // Display the solution
            System.out.println("Solution:");
            for (int i = 0; i < numVariables; i++) {
                System.out.println("x" + (i + 1) + " = " + solution[i]);
            }
        // Scanner is automatically closed when the try block is exited
    }
    

    private static double[] solveLinearEquations(double[][] coefficients, double[] constants) {
        int numVariables = coefficients.length;
        double[] solution = new double[numVariables];

        double determinant = calculateDeterminant(coefficients);

        if (determinant == 0) {
            System.out.println("The system of equations is singular. No unique solution exists.");
            return solution;
        }

        for (int i = 0; i < numVariables; i++) {
            double[][] modifiedMatrix = modifyMatrix(coefficients, constants, i);
            solution[i] = calculateDeterminant(modifiedMatrix) / determinant;
        }

        return solution;
    }


    private static double[][] modifyMatrix(double[][] coefficients, double[] constants, int columnIndex) {
        int numVariables = coefficients.length;
        double[][] modifiedMatrix = new double[numVariables][numVariables];

        for (int i = 0; i < numVariables; i++) {
            for (int j = 0; j < numVariables; j++) {
                if (j == columnIndex) {
                    modifiedMatrix[i][j] = constants[i];
                } else {
                    modifiedMatrix[i][j] = coefficients[i][j];
                }
            }
        }
        return modifiedMatrix;
    }

    private static double calculateDeterminant(double[][] matrix) {
        int order = matrix.length;
        if (order == 1) {
            return matrix[0][0];
        } else if (order == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        } else {
            double determinant = 0;
            for (int i = 0; i < order; i++) {
                determinant += Math.pow(-1, i) * matrix[0][i] * calculateDeterminant(getSubmatrix1(matrix, 0, i));
            }
            return determinant;
        }
    }
    private static double[][] getSubmatrix1(double[][] matrix, int rowToRemove, int colToRemove) {
        int order = matrix.length - 1;
        double[][] submatrix = new double[order][order];
        int rowIndex = 0;
    
        for (int i = 0; i < matrix.length; i++) {
            if (i == rowToRemove) continue;
            int colIndex = 0;
    
            for (int j = 0; j < matrix[i].length; j++) {
                if (j == colToRemove) continue;
                submatrix[rowIndex][colIndex] = matrix[i][j];
                colIndex++;
            }
            rowIndex++;
        }
    
        return submatrix;
    }
    

    private static void runPolynomialEquationSolver(Scanner scanner) {
        
            System.out.println("Polynomial Equation Solver");
    
            System.out.println("Enter the degree of the polynomial equation (2 for quadratic, 3 for cubic, etc.):");
            int degree = scanner.nextInt();
    
            double[] coefficients = new double[degree + 1];
    
            if (degree != 2 && degree != 3) {
                System.out.println("Polynomial degree not supported. Please use a library or numerical methods for other degrees.");
                return;
            }
    
            System.out.println("Enter the coefficients of the polynomial equation from the highest degree to the constant term:");
    
            for (int i = degree; i >= 0; i--) {
                System.out.println("Enter the coefficient for x^" + i + ":");
                coefficients[i] = scanner.nextDouble();
            }
    
            System.out.println("Polynomial equation: " + generatePolynomialEquation(coefficients));
    
            // Solve the polynomial equation
            double[] roots = solvePolynomialEquation(coefficients);
    
            // Display the roots
            System.out.println("Roots:");
            for (double root : roots) {
                System.out.println(root);
            }
        // Scanner is automatically closed when the try block is exited
    }
    
    private static String generatePolynomialEquation(double[] coefficients) {
        StringBuilder equation = new StringBuilder();

        for (int i = coefficients.length - 1; i >= 0; i--) {
            if (coefficients[i] != 0) {
                if (equation.length() > 0) {
                    equation.append(" + ");
                }
                equation.append(coefficients[i]);
                if (i > 0) {
                    equation.append("x^").append(i);
                }
            }
        }

        return equation.toString();
    }

    private static double[] solvePolynomialEquation(double[] coefficients) {
        int degree = coefficients.length - 1;

        if (degree == 2) {
            // Quadratic equation, use quadratic formula
            double a = coefficients[2];
            double b = coefficients[1];
            double c = coefficients[0];

            double discriminant = b * b - 4 * a * c;

            if (discriminant > 0) {
                double root1 = (-b + Math.sqrt(discriminant)) / (2 * a);
                double root2 = (-b - Math.sqrt(discriminant)) / (2 * a);
                return new double[]{root1, root2};
            } else if (discriminant == 0) {
                double root = -b / (2 * a);
                return new double[]{root};
            } else {
                // Complex roots
                return new double[]{};
            }
        } else if (degree == 3) {
            // Cubic equation, use Cardano's method (an example, not exhaustive)
            double a = coefficients[3];
            double b = coefficients[2];
            double c = coefficients[1];
            double d = coefficients[0];

            // (Implementation of cubic equation solving can be complex)
            // Here, we use a simple example of Cardano's method
            // For a more accurate and comprehensive solution, consider using libraries

            double p = c / a - b * b / (3 * a * a);
            double q = 2 * b * b * b / (27 * a * a * a) - b * c / (3 * a * a) + d / a;

            double discriminant = q * q / 4 + p * p * p / 27;

            if (discriminant > 0) {
                // 1 real root and 2 complex roots
                double root1 = -2 * Math.cbrt(-q / 2 + Math.sqrt(discriminant)) - b / (3 * a);
                double root2 = Math.cbrt(-q / 2 - Math.sqrt(discriminant)) - b / (3 * a);
                return new double[]{root1, root2};
            } else if (discriminant == 0) {
                // 3 real roots (1 triple root)
                double root = -2 * Math.cbrt(-q / 2) - b / (3 * a);
                return new double[]{root};
            } else {
                // 3 real roots
                double r = Math.sqrt(-p * p * p / 27);
                double theta = Math.acos(-q / (2 * r));
                double root1 = 2 * Math.cbrt(r) * Math.cos(theta / 3) - b / (3 * a);
                double root2 = 2 * Math.cbrt(r) * Math.cos((theta + 2 * Math.PI) / 3) - b / (3 * a);
                double root3 = 2 * Math.cbrt(r) * Math.cos((theta + 4 * Math.PI) / 3) - b / (3 * a);
                return new double[]{root1, root2, root3};
            }
        } else {
            // For higher-degree polynomials, you may use numerical methods or external libraries.
            System.out.println("Polynomial degree not supported. Please use a library or numerical methods for higher degrees.");
            return new double[]{};
        }
    }
    private static void runMatrixOperations(Scanner scanner) {
        System.out.print("Enter the number of rows for the matrices: ");
        int rows = scanner.nextInt();

        System.out.print("Enter the number of columns for the matrices: ");
        int columns = scanner.nextInt();

        int[][] matrixA = new int[rows][columns];
        int[][] matrixB = new int[rows][columns];

        System.out.println("Enter elements for Matrix A:");
        inputMatrix(matrixA, scanner);

        System.out.println("Enter elements for Matrix B:");
        inputMatrix(matrixB, scanner);

        System.out.println("Matrix A:");
        printMatrix(matrixA);

        System.out.println("Matrix B:");
        printMatrix(matrixB);
        while(true){
            System.out.println("Select a matrix operation:");
            System.out.println("1. Matrix Addition");
            System.out.println("2. Matrix Subtraction");
            System.out.println("3. Matrix Multiplication");
            System.out.println("4. Matrix Transpose");
            System.out.println("5. Matrix Determinant");
            
            int operationChoice;
            try {
                operationChoice = scanner.nextInt();
            } catch (InputMismatchException e) {
                System.out.println("InputMismatchException");
                scanner.nextLine(); // Consume the invalid input
                continue; // Restart the loop
            }
            if(operationChoice==6){
                System.out.println("returning to main menu.");
                break;
            }

            
            switch (operationChoice) {
                case 1:
                    int[][] sumMatrix = addMatrices(matrixA, matrixB);
                    System.out.println("Sum Matrix:");
                    printMatrix(sumMatrix);
                    break;
                case 2:
                    int[][] differenceMatrix = subtractMatrices(matrixA, matrixB);
                    System.out.println("Difference Matrix:");
                    printMatrix(differenceMatrix);
                    break;
                case 3:
                    int[][] productMatrix = multiplyMatrices(matrixA, matrixB);
                    System.out.println("Product Matrix:");
                    printMatrix(productMatrix);
                    break;
                case 4:
                    int[][] transposedMatrixA = transposeMatrix(matrixA);
                    int[][] transposedMatrixB = transposeMatrix(matrixB);
                    System.out.println("Transposed Matrix A:");
                    printMatrix(transposedMatrixA);
                    System.out.println("Transposed Matrix B:");
                    printMatrix(transposedMatrixB);
                    break;
                case 5:
                    int determinantA = calculateDeterminant(matrixA);
                    int determinantB = calculateDeterminant(matrixB);
                    System.out.println("Determinant of Matrix A: " + determinantA);
                    System.out.println("Determinant of Matrix B: " + determinantB);
                    break;
                default:
                    System.out.println("Invalid choice. Exiting.");
            }
        }
    }
    private static void inputMatrix(int[][] matrix, Scanner scanner) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print("Enter element at position (" + (i + 1) + ", " + (j + 1) + "): ");
                matrix[i][j] = scanner.nextInt();
            }
        }
    }

    private static void printMatrix(int[][] matrix) {
        if (matrix == null) {
            System.out.println("Matrix is null. Cannot print.");
            return;
        }
    
        for (int[] row : matrix) {
            for (int value : row) {
                System.out.print(value + "\t");
            }
            System.out.println();
        }
    }
    

    private static int[][] addMatrices(int[][] matrixA, int[][] matrixB) {
        int rows = matrixA.length;
        int columns = matrixA[0].length;
        int[][] result = new int[rows][columns];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[i][j] = matrixA[i][j] + matrixB[i][j];
            }
        }
        return result;
    }

    private static int[][] subtractMatrices(int[][] matrixA, int[][] matrixB) {
        int rowsA = matrixA.length;
        int columnsA = matrixA[0].length;
        int columnsB = matrixB[0].length;
        int[][] result = new int[rowsA][columnsB];
    
        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < columnsB; j++) {
                result[i][j] = matrixA[i][j] - matrixB[i][j];  // Fix: Subtract matrixB elements
            }
        }
        return result;
    }
    
    private static int[][] multiplyMatrices(int[][] matrixA, int[][] matrixB) {
        try{
            int rowsA = matrixA.length;
            int columnsA = matrixA[0].length;
            int rowsB = matrixB.length;
            int columnsB = matrixB[0].length;
        
            // Print matrix dimensions for debugging
            System.out.println("Matrix A dimensions: " + rowsA + "x" + columnsA);
            System.out.println("Matrix B dimensions: " + rowsB + "x" + columnsB);
        
            // Check if matrix multiplication is defined
            if (columnsA != rowsB) {
                System.out.println("Error: Matrix multiplication is not defined. Number of columns in Matrix A must be equal to the number of rows in Matrix B.");
                
            }
        
            int[][] result = new int[rowsA][columnsB];
        
            for (int i = 0; i < rowsA; i++) {
                for (int j = 0; j < columnsB; j++) {
                    for (int k = 0; k < columnsA; k++) {
                        result[i][j] += matrixA[i][k] * matrixB[k][j];
                    }
                }
            }
             return result;
        }
       
        catch(Exception e){
            System.out.println("IllegalArgumentException");
            return null;
        }
    }
    
    

    private static int[][] transposeMatrix(int[][] matrix) {
        int rows = matrix.length;
        int columns = matrix[0].length;
        int[][] result = new int[columns][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    private static int calculateDeterminant(int[][] matrix) {
        int order = matrix.length;
        if (order == 1) {
            return matrix[0][0];
        } else if (order == 2) {
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        } else {
            int determinant = 0;
            for (int i = 0; i < order; i++) {
                determinant += Math.pow(-1, i) * matrix[0][i] * calculateDeterminant(getSubmatrix(matrix, 0, i));
            }
            return determinant;
        }
    }

    private static double[][] getSubmatrix(int[][] matrix, int rowToRemove, int colToRemove) {
        int order = matrix.length - 1;
        double[][] submatrix = new double[order][order];
        int rowIndex = 0;
    
        for (int i = 0; i < matrix.length; i++) {
            if (i == rowToRemove) continue;
            int colIndex = 0;
    
            for (int j = 0; j < matrix[i].length; j++) {
                if (j == colToRemove) continue;
                submatrix[rowIndex][colIndex] = matrix[i][j];
                colIndex++;
            }
            rowIndex++;
        }
    
        return submatrix;
    }
    

    private static void runTrigonometryCalculator(Scanner scanner) {

        System.out.println("Trigonometry Calculator");
    
        while (true) {
            System.out.println("Enter the angle in radians (or type 'exit' to quit):");
            String input = scanner.nextLine();
    
            if (input.equalsIgnoreCase("exit")) {
                System.out.println("Exiting the calculator. Goodbye!");
                break;
            }
    
            try {
                double angle = Double.parseDouble(input);
                double radians = Math.toRadians(angle);
    
                // Check if the angle is exactly 90 degrees
                if (angle % 180 == 90) {
                    System.out.println("Angle: " + angle + " degrees");
                    System.out.println("Sine: 1.0 (undefined)");
                    System.out.println("Cosine: 0.0");
                    System.out.println("Tangent: undefined (division by zero)");
                } else {
                    // Calculate trigonometric functions
                    double sine = Math.sin(radians);
                    double cosine = Math.cos(radians);
                    double tangent = Math.tan(radians);
    
                    // Display results
                    System.out.println("Angle: " + angle + " degrees");
                    System.out.println("Sine: " + sine);
                    System.out.println("Cosine: " + cosine);
                    System.out.println("Tangent: " + tangent);
                }
            } catch (NumberFormatException e) {
                System.out.println("Invalid input. Please enter a valid angle.");
            }
        }
    }
    

}
