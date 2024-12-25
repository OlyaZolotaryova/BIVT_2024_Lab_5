using System;
using System.Collections;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data.Common;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static Program;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }
    #region Level 1
    static int Factorial(int n)
    {
        int f = 1;
        for (int i = 2; i <= n; i++)
            f = f * i;
        return f;
    }
    static int Combinations(int n, int k)
    {
        int cnk = Factorial(n) / (Factorial(k) * Factorial(n - k));
        return cnk;
    }
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here

        if (k < 0 || n < 0) return 0;
        answer=Combinations (n, k);
        // create and use Combinations(n, k);
        // create and use Factorial(n);

        // end

        return answer;
    }
    static double GeronArea(double a, double b, double c)
    {
        if (a <= 0 || b <= 0 || c <= 0) return default;
        else
        {
            double p = (a + b + c) / 2;
            double S = Math.Sqrt(p * (p - a) * (p - b) * (p - c));
            return S;
        }
    }
    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        double a1=first[0];
        double b1 = first[1];
        double c1 = first[2];
        double a2 = second[0];
        double b2 = second[1];
        double c2 = second[2];
        if (a1<= 0 || b1<=0 || c1<=0 || a2 <= 0 || b2 <= 0 || c2 <= 0 || a1+b1<=c1 || a2 + b2 <= c2 || a1 + c1 <= b1 || a2 + c2 <= b2 || b1 + c1 <= a1 || b2 + c2 <= a2) return -1;
        if (GeronArea(a1, b1, c1) > GeronArea(a2, b2, c2)) return 1;
        if (GeronArea(a1, b1, c1) < GeronArea(a2, b2, c2)) return 2;
        else return 0;
        // create and use GeronArea(a, b, c);

        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }
    static double GetDistance(double v, double a, int t)
    {
        double s = v * t + a * t * t / 2;
        return s;
    }
    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here

        // create and use GetDistance(v, a, t); t - hours
        if (GetDistance(v1, a1, time) > GetDistance(v2, a2, time)) return 1;
        if (GetDistance(v1, a1, time) < GetDistance(v2, a2, time)) return 2;
        else return 0;
        // end

        // first = 1, second = 2, equal = 0
        return answer;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here
        for (int t=1; ;t++)
            if(GetDistance(v2, a2, t) >= GetDistance(v1, a1, t)) return t;
        // use GetDistance(v, a, t); t - hours

        // end

        return answer;
    }
    #endregion

    #region Level 2
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here
        int iA= 0, jA= 0, iB=0, jB=0;
        FindMaxIndex(A, out iA, out jA);
        int Amax= A[iA, jA];
        FindMaxIndex(B, out iB, out jB);
        int Bmax = B[iB, jB];
        A[iA, jA] = Bmax;
        B[iB, jB] = Amax;
        // create and use FindMaxIndex(matrix, out row, out column);

        // end
    }
    static void FindMaxIndex(int[,] matrix, out int row, out int column)
    {
        int max = matrix[0, 0];
        row = 0;
        column = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    row = i;
                    column = j;
                }
            }
    }
    public void Task_2_2(double[] A, double[] B)
    {
        // code here

        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!

        // end
    }
    static int FindDiagonalMaxIndex(int [,] matrix)
    {
        int ii = 0, max = matrix[0,0];
        for (int i=0; i<matrix.GetLength(0); i++)
            if (matrix[i,i] > max)
            {
                max = matrix[i,i];
                ii = i;
            }
        return ii;
    }
    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here
        int iB = FindDiagonalMaxIndex(B);
        int iC = FindDiagonalMaxIndex(C);
        int[,] b=new int[B.GetLength(0)-1, B.GetLength(0)];
        int[,] c = new int[C.GetLength(0) - 1, C.GetLength(0)];
        for (int i = 0; i < B.GetLength(0)-1; i++)
            for (int j = 0; j < B.GetLength(0); j++)
            {
                if (i < iB) b[i, j] = B[i, j];
                else b[i, j] = B[i + 1, j];
            }
        for (int i = 0; i < C.GetLength(0) - 1; i++)
            for (int j = 0; j < C.GetLength(0); j++)
            {
                if (i < iC) c[i, j] = C[i, j];
                else c[i, j] = C[i + 1, j];
            }
        B = b;
        C = c;
        //  create and use method FindDiagonalMaxIndex(matrix);

        // end
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3

        // end
    }
    static void FindMaxInColumn(int [,] matrix, int columnIndex, out int rowIndex)
    {
        int max = matrix[0, columnIndex];
        rowIndex = 0;
        for (int i = 0; i<matrix.GetLength(0); i++)
            if (matrix[i, columnIndex] >max)
            {
                max = matrix[i, columnIndex];
                rowIndex = i;
            }
    }
    public void Task_2_5(int[,] A, int[,] B)
    {
    // code here
    int column = 0;
        int ia = 0;
        int ib = 0;
        FindMaxInColumn(A, column, out ia);
        FindMaxInColumn(B, column, out ib);
        for (int j = 0; j < B.GetLength(0); j++)
        {
            int t = A[ia, j];
            A[ia, j] = B[ib, j];
            B[ib, j] = t;
        }
        // create and use FindMaxInColumn(matrix, columnIndex, out rowIndex);

        // end
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here

        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);

        // end
    }
    static int CountRowPositive(int[,] matrix, int rowIndex)
    {
        int  k= 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
            if (matrix[rowIndex, j] > 0) k++;
        return k;
    }
    static int CountColumnPositive(int[,] matrix, int colIndex)
    {
        int k= 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            if (matrix[i, colIndex] > 0) k++;
        return k;
    }
    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here
        int maxB = 0, maxrow=0;
        for (int i = 0; i < B.GetLength(0); i++)
        {
            int c = CountRowPositive(B, i); 
            if (c > maxB)
            {
                maxB = c;
                maxrow = i;
            }
        }
        int maxC = 0, maxcol=0;
        for (int j = 0; j < C.GetLength(1); j++)
        {
            int c = CountColumnPositive(C, j);
            if (c > maxC)
            {
                maxC = c;
                maxcol = j;
            }
        }
        int[,] b = new int[B.GetLength(0) + 1, B.GetLength(1)];

        for (int i = 0; i < B.GetLength(0); i++)
        {
            if (i <= maxrow)
                for (int j = 0; j < B.GetLength(1); j++)
                    b[i, j] = B[i, j];
            else
                for (int j = 0; j < B.GetLength(1); j++)
                    b[i + 1, j] = B[i, j];
        }
        for (int j = 0; j < B.GetLength(1); j++)
            b[maxrow+1, j] = C[j, maxcol];
        B = b;
        // create and use CountRowPositive(matrix, rowIndex);
        // create and use CountColumnPositive(matrix, colIndex);

        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here

        // create and use SortArrayPart(array, startIndex);

        // end
    }
    static int [] SumPositiveElementsInColumns(int [,] matrix)
    {
        int[] a = new int[matrix.GetLength(1)];
        int k = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            int s = 0;
            for (int i = 0; i < matrix.GetLength(0); i++)
                if (matrix[i, j] > 0)
                    s += matrix[i, j];
            a[k] = s;
            k++;
        }
        return a;
    }
    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here
        int[] a = new int[A.GetLength(1)];
        int[] c = new int[C.GetLength(1)];
        int[] b = new int[a.Length+c.Length];
        a = SumPositiveElementsInColumns(A);
        c = SumPositiveElementsInColumns(C);
        int k = 0;
        for (int i=0; i<a.Length+c.Length; i++)
        {
            if (i<a.Length) b[i] = a[i];
            else
            {
                b[i] = c[k];
                k++;
            }
        }
        answer = b;
        // create and use SumPositiveElementsInColumns(matrix);

        // end

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here

        // create and use RemoveColumn(matrix, columnIndex);

        // end
    }

    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here
        int iA = 0, jA = 0, iB = 0, jB = 0;
        FindMaxIndex(A, out iA, out jA);
        int Amax = A[iA, jA];
        FindMaxIndex(B, out iB, out jB);
        int Bmax = B[iB, jB];
        A[iA, jA] = Bmax;
        B[iB, jB] = Amax;
        // use FindMaxIndex(matrix, out row, out column); from Task_2_1

        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxColumnIndex(matrix);

        // end
    }
    static int [,] RemoveRow(int [,] matrix, int rowIndex)
    {
        int[,] a = new int[matrix.GetLength(0)-1, matrix.GetLength(1)];
        for (int i = 0; i < matrix.GetLength(0)-1; i++)
        {
            if (i < rowIndex)
                for (int j = 0; j < matrix.GetLength(1); j++)
                    a[i, j] = matrix[i, j];
            else
                for (int j = 0; j < matrix.GetLength(1); j++)
                    a[i, j] = matrix[i+1, j];
        }
        matrix = a;
        return matrix;
    }
    public void Task_2_13(ref int[,] matrix)
    {
        // code here
        int max = -1000000, min= 1000000, imax=0, imin=0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0;j < matrix.GetLength(1); j++)
            {
                if (matrix[i,j] > max)
                {
                    max = matrix[i,j];
                    imax = i;
                }
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    imin = i;
                }
            }
        if (imin == imax)
            matrix=RemoveRow(matrix, imin);
        else
        {
            if (imin > imax)
            {
                matrix=RemoveRow(matrix, imin);
                matrix=RemoveRow(matrix, imax);
            }
            else
            {
                matrix=RemoveRow(matrix, imax);
                matrix=RemoveRow(matrix, imin);
            }
        }
        // create and use RemoveRow(matrix, rowIndex);

        // end
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here

        // create and use SortRow(matrix, rowIndex);

        // end
    }
    static int GetAverageWithoutMinMax(int [,] matrix)
    {
        int max = matrix[0, 0];
        int min = matrix[0, 0];
        int s = 0;
        int sr = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                    max = matrix[i, j];
                if (matrix[i, j] < min)
                    min = matrix[i, j];
                s=s+matrix[i, j];
            }
        if (max == min) sr = (s-max) / (matrix.GetLength(0) * matrix.GetLength(1) - 1);
        else sr = (s-max-min) / (matrix.GetLength(0) * matrix.GetLength(1) - 2);
        return sr;
    }
    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here
        int[] a = new int[3];
        a[0] = GetAverageWithoutMinMax(A);
        a[1] = GetAverageWithoutMinMax(B);
        a[2] = GetAverageWithoutMinMax(C);
        if (a[0] > a[1] && a[1] > a[2])  return -1;
        if (a[0] < a[1] && a[1] < a[2]) return 1;
        else return 0;
            // create and use GetAverageWithoutMinMax(matrix);

            // end

            // 1 - increasing   0 - no sequence   -1 - decreasing
            return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here

        // create and use SortNegative(array);

        // end
    }
    static int [,] SortRowsByMaxElement(int [,] matrix)
    {
        for (int i = 1; i < matrix.GetLength(0);)
        {
            int max1 = -10000000, max2 = -100000;
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (i != 0)
                {
                    if (matrix[i, j] >= max2) max2 = matrix[i, j];
                    if (matrix[i - 1, j] >= max1) max1 = matrix[i - 1, j];
                }
            }
            if (i == 0 || max1 >= max2)
                i++;
            else
            {
                for (int j = 0; j < matrix.GetLength(1); j++)
                {
                    int t = matrix[i, j];
                    matrix[i, j] = matrix[i - 1, j];
                    matrix[i - 1, j] = t;
                }
                i--;
            }
        }
        return matrix;
    }
    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here
        A = SortRowsByMaxElement(A);
        B = SortRowsByMaxElement(B);
        // create and use SortRowsByMaxElement(matrix);

        // end
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here

        // create and use SortDiagonal(matrix);

        // end
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here
        for (int i= matrix.GetLength(0)-1; i>=0; i--)
        {
            int k = 0;
            for (int j = 0; j < matrix.GetLength(1); j++)
                if (matrix[i, j] == 0) k++;
            if (k!=0) matrix= RemoveRow(matrix, i);
        }
        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here

        // use RemoveColumn(matrix, columnIndex); from 2_10

        // end
    }
    static int [] CreateArrayFromMins(int [,] matrix)
    {
        int[] a = new int[matrix.GetLength(1)];
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int min= matrix[i,i];
            for (int j = i; j<matrix.GetLength(1); j++)
                if (matrix[i,j] < min) min = matrix[i,j];
            a[i] = min;
        }
        return a;
    }
    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here
        answerA= new int[A.GetLength(0)];
        answerB = new int[B.GetLength(0)];
        answerA = CreateArrayFromMins(A);
        answerB = CreateArrayFromMins(B);
        // create and use CreateArrayFromMins(matrix);

        // end
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here

        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);

        // end
    }
    static double [,] MatrixValuesChange(double [,] matrix)
    {
        double [] a = new double[matrix.GetLength(0)*matrix.GetLength(1)];
        int k = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                a[k] = matrix[i,j];
                k++;
            }
        for (int i = 1; i < a.Length;)
        {
            if (i == 0 || a[i] >= a[i - 1])
                i++;
            else
            {
                double t = a[i];
                a[i] = a[i - 1];
                a[i - 1] = t;
                i--;
            }
        }
        if (a.Length>5)
        {
            double[] b = new double[5];
            k = 0;
            for (int i = a.Length-5;i < a.Length; i++)
            {
                b[k] = a[i];
                k++;
            }
            a = b;
        }
        double[,] c = new double[matrix.GetLength(0), matrix.GetLength(1)];
        c = matrix;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (c[i, j] == a[0] || c[i, j] == a[1] || c[i, j] == a[2] || c[i, j] == a[3] || c[i, j] == a[4])
                    {
                        if (c[i, j] >= 0)
                            c[i, j] = 2 * c[i, j];
                        else
                            c[i, j] = c[i, j] / 2;
                    }
                else
                {
                    if (c[i, j] < 0)
                        c[i, j] = 2 * matrix[i, j];
                    else
                        c[i, j] = matrix[i, j] / 2;
                }
            }
        matrix = c;
        return matrix;
    }
    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here
        A= MatrixValuesChange(A);
        B= MatrixValuesChange(B);
        // create and use MatrixValuesChange(matrix);
        // end
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);

        // end
    }
    static int FindRowWithMaxNegativeCount(int [,] matrix)
    {
        int max = 0;
        int imax = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int k=CountNegativeInRow(matrix, i);
            if ( k> max)
            {
                max = k;
                imax = i;
            }
        }
        return imax;
    }
    static int CountNegativeInRow(int [,] matrix, int rowIndex)
    {
        int k = 0;
        for (int j = 0; j < matrix.GetLength(1); j++)
            if (matrix[rowIndex, j] < 0) k++;
        return k;
    }
    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here
        maxA = FindRowWithMaxNegativeCount(A);
        maxB = FindRowWithMaxNegativeCount(B);
        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22

        // end
    }
    static int FindRowMaxIndex(int [,] matrix, int rowIndex, out int columnIndex)
    {
            columnIndex = 0;
            int max = matrix[rowIndex, columnIndex];
            for(int j = 0;j < matrix.GetLength(1); j++)
                if (matrix[rowIndex, j] > max)
                {
                    max = matrix[rowIndex, j];
                    columnIndex = j;
                }
        return columnIndex;
    }
    static int ReplaceMaxElementOdd(int [,] matrix, int row, int column)
    {
        matrix[row, column] = matrix[row, column] * (column+1);
        return matrix[row, column];
    }
    static int ReplaceMaxElementEven (int [,] matrix, int row, int column)
    {
        matrix[row, column] = 0;
        return matrix[row, column];
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here
        for (int i = 0; i < A.GetLength(0); i++)
        {
            int columnIndex = 0;
            int Aj = FindRowMaxIndex(A, i, out columnIndex);
            if (i % 2 != 0)
                A[i, Aj] = ReplaceMaxElementEven(A, i, Aj);
            else
                A[i, Aj] = ReplaceMaxElementOdd(A, i, Aj);
        }
        for (int i = 0; i < B.GetLength(0); i++)
        {
            int columnIndex = 0;
            int Bj = FindRowMaxIndex(B, i, out columnIndex);
            if (i % 2 != 0)
                B[i, Bj] = ReplaceMaxElementEven(B, i, Bj);
            else
                B[i, Bj] = ReplaceMaxElementOdd(B, i, Bj);
        }
        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search

        // end
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here

        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search

        // end
    }
    #endregion
    public delegate double SumFunction(int i, double x, ref int c);
    public delegate double YFunction(double x);
    public double S1(int i, double x, ref int fac)
    {
        if (i > 0)
            fac = fac*i;
        double s= Math.Cos(i * x) / fac;
        return s;
    }
    public double S2(int i, double x, ref int k)
    {
        k =k* (-1);
        double s= k * Math.Cos(i * x) / (i * i);
        return s;
    }
    public double Y1(double x)
    {
        double y= Math.Exp(Math.Cos(x)) * Math.Cos(Math.Sin(x));
        return y;
    }
    public double Y2(double x)
    {
        double y= ((x * x) - Math.PI * Math.PI / 3) / 4;
        return y;
    }
    public double CalculateSum(SumFunction sumFunction, double x, int i)
    {
        double e = 0.0001, s = 0;
        int ch = 1;
        double c = sumFunction(i, x, ref ch);

        while (Math.Abs(c) > e)
        {
            s =s+ c;
            c = sumFunction(++i, x, ref ch);
        }
        return s;
    }
    public void GetSumAndY(SumFunction sFunction, YFunction yFunction, double a, double b, double h, double[,] sy, int st = 0)
    {
        for (int i = 0; i < (b - a) / h + 1; i=i+1)
        {
            double x = a + i * h;
            double s = CalculateSum(sFunction, x, st), y = yFunction(x);
            sy[i, 0] = s;
            sy[i, 1] = y;
        }
    }
    #region Level 3
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        //code here
        double a1 = 0.1, b1 = 1, h1 = 0.1;
        firstSumAndY = new double[(int)((b1 - a1) / h1) + 1, 2];
        GetSumAndY(S1, Y1, a1, b1, h1, firstSumAndY);
        double a2 = Math.PI / 5, b2 = Math.PI, h2 = Math.PI / 25;
        secondSumAndY = new double[(int)((b2 - a2) / h2) + 1, 2];
        GetSumAndY(S2, Y2, a2, b2, h2, secondSumAndY, 1);
        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }

    public void Task_3_2(int[,] matrix)
    {
        // SortRowStyle sortStyle = default(SortRowStyle); - uncomment me

        // code here

        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting

        // end
    }
    public delegate void SwapDirection(double[] array);
    public void SwapRight(double [] array)
    {
        for (int i = 0; i < array.Length - 1; i = i+2)
            (array[i], array[i + 1]) = (array[i + 1], array[i]);
    }
    public void SwapLeft(double[] array)
    {
        for (int i = array.Length-1; i >0; i= i-2)
            (array[i], array[i - 1]) = (array[i - 1], array[i]);
    }
    public double GetSum(double [] array, int start, int h)
    {
        double s = 0;
        for (int i = start; i < array.Length; i= i+h)
            s=s+ array[i];
        return s;
    }
    public double Task_3_3(double[] array)
    {
        double answer = 0;
        SwapDirection swapper = default(SwapDirection);
        double s = 0, sr = 0;
        for (int i = 0; i<array.Length; i = i+1)
            s= s+ array[i];
        sr=s/array.Length;
        if (array[0] > sr) swapper = SwapRight;
        else swapper = SwapLeft;
        swapper(array);
        answer = GetSum(array, 1, 2);


        // code here

        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here

        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end

        return answer;
    }
    static int F1(double a, double b, double h)
    {
        int k = 1, t=0;
        double y = a * a - Math.Sin(a);
        if (y < 0) t = -1;
        else t = 1;

        for (double x = a + h; x <= b; x = x+ h)
        {
            y = x * x + Math.Sin(x);
            if (((y < 0) && (t > 0)) || ((y > 0) && (t < 0))) k=k+1;
            if (y < 0) t = -1;
            else t = 1;
        }
        return k;
    }

    static int F2(double a, double b, double h)
    {
        int k = 1, t=0;
        double y = Math.Exp(a) - 1;
        if (y < 0) t = -1;
        else t = 1;
        for (double x = a + h; x <= b; x = x+h)
        {
            y = Math.Exp(x) - 1;
            if (((y< 0) && (t > 0)) || ((y > 0) && (t < 0))) k=k+1;
            if (y < 0) t= -1;
            else t = 1;
        }
        return k;
    }
    public delegate int Function(double a, double b, double h);
    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;
        //code here
        
        Function function = F1;
        func1 = function(0, 2, 0.1);
        Function function1 = F2;
        func2 = function1(-1, 1, 0.2);

        //use public delegate YFunction(x, a, b, h) from Task_3_1
        //create and use method CountSignFlips(YFunction, a, b, h);
        //create and use 2 methods for both functions
        //end
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here

        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);

        // end
    }
    static int [,] InsertColumn(int [,] B, CountPositive Row, int [,] C, CountPositive Column)
    {
        int maxB = 0, maxrow = 0;
        for (int i = 0; i < B.GetLength(0); i++)
        {
            int c = CountRowPositive(B, i);
            if (c > maxB)
            {
                maxB = c;
                maxrow = i;
            }
        }
        int maxC = 0, maxcol = 0;
        for (int j = 0; j < C.GetLength(1); j++)
        {
            int c = CountColumnPositive(C, j);
            if (c > maxC)
            {
                maxC = c;
                maxcol = j;
            }
        }
        int[,] b = new int[B.GetLength(0) + 1, B.GetLength(1)];

        for (int i = 0; i < B.GetLength(0); i++)
        {
            if (i <= maxrow)
                for (int j = 0; j < B.GetLength(1); j++)
                    b[i, j] = B[i, j];
            else
                for (int j = 0; j < B.GetLength(1); j++)
                    b[i + 1, j] = B[i, j];
        }
        for (int j = 0; j < B.GetLength(1); j++)
            b[maxrow + 1, j] = C[j, maxcol];
        B = b;
        return B;
    }
    public delegate int CountPositive(int [,] matrix, int index);
    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here
        CountPositive row = CountRowPositive;
        CountPositive column = CountColumnPositive;

        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);
        B = InsertColumn(B, row, C, column);
        // end
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // FindIndex searchArea = default(FindIndex); - uncomment me

        // code here

        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)

        // end
    }
    static int [,] RemoveRows(int [,] matrix, FindElementDelegate Max, FindElementDelegate Min)
    {
        int findMax = Max(matrix);
        int findMin = Min(matrix);
        if (findMin == findMax)
            matrix = RemoveRow(matrix, findMin);
        else
        {
            if (findMin > findMax)
            {
                matrix = RemoveRow(matrix, findMin);
                matrix = RemoveRow(matrix, findMax);
            }
            else
            {
                matrix = RemoveRow(matrix, findMax);
                matrix = RemoveRow(matrix, findMin);
            }
        }
        return matrix;
    }
    public delegate int FindElementDelegate(int[,] matrix);
    public void Task_3_13(ref int[,] matrix)
    {
        // code here
        FindElementDelegate findMax = FindMaxIndex;
        FindElementDelegate findMin = FindMinIndex;
        matrix = RemoveRows(matrix, findMax, findMin);
        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }
    static int FindMinIndex(int[,] matrix)
    {
        int min = matrix[0, 0];
        int row = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    row = i;
                }
            }
        return row;
    }
    static int FindMaxIndex(int[,] matrix)
    {
        int max = matrix[0, 0];
        int row = 0;
        for (int i = 0; i < matrix.GetLength(0); i++)
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    row = i;
                }
            }
        return row;
    }
    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {

        rows = null;
        cols = null;

        // code here

        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);

        // end
    }
    static int[,] EvenOddRowsTransform(int[,] matrix, ReplaceMaxElement Odd, ReplaceMaxElement Even)
    {
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            int columnIndex = 0;
            int Aj = FindRowMaxIndex(matrix, i, out columnIndex);
            if (i % 2 != 0)
                matrix[i, Aj] = Even(matrix, i, Aj);
            else
                matrix[i, Aj] = Odd(matrix, i, Aj);
        }
        return matrix;
    }   
    public delegate int ReplaceMaxElement(int[,] matrix, int rowIndex, int max);
    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here
        ReplaceMaxElement odd = ReplaceMaxElementOdd;
        ReplaceMaxElement even = ReplaceMaxElementEven;
        A = EvenOddRowsTransform(A, odd, even);
        B = EvenOddRowsTransform(B, odd, even);
        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method (EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);

        // end
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here

        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);

        // end
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
