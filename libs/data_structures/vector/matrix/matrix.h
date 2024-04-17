#ifndef INC_MATRIX_H
#define INC_MATRIX_H

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stddef.h>
#include <stdbool.h>
#include <string.h>

typedef struct matrix {
    int **values; // элементы матрицы
    int nRows; // количество рядов
    int nCols; // количество столбцов
} matrix;

typedef struct position {
    int rowIndex;
    int colIndex;
} position;

//размещает в динамической памяти матрицу размером nRows на nCols. Возвращает матрицу.
matrix getMemMatrix(int nRows, int nCols);

//размещает в динамической памяти массив из nMatrices матриц размером nRows на nCols.
// Возвращает указатель на нулевую матрицу.
matrix *getMemArrayOfMatrices(int nMatrices, int nRows, int nCols);

//освобождает память, выделенную под хранение матрицы m
void freeMemMatrix(matrix *m);

//освобождает память, выделенную под хранение массива ms из nMatrices матриц.
void freeMemMatrices(matrix *ms, int nMatrices);

// ввод матрицы m.
void inputMatrix(matrix *m);

// ввод массива из nMatrices матриц, хранящейся по адресу ms.
void inputMatrices(matrix *ms, int nMatrices);

//вывод матрицы m.
void outputMatrix(matrix m);

// считает количество цифр в числе.
int countDigits(int n);

//находит число в матрице m с самым большим количеством разрядов.
int findLongNumberInMatrix(matrix m);

//вывод массива из nMatrices матриц, хранящейся по адресу ms.
void outputMatrices(matrix *ms, int nMatrices);

//обмен строк с порядковыми номерами i1 и i2 в матрице m
void swapRows(matrix m, int i1, int i2);

//обмен колонок с порядковыми номерами j1 и j2 в матрице m.
void swapColumns(matrix m, int j1, int j2);

// выполняет сортировку вставками строк матрицы m по неубыванию значения функции criteria применяемой для строк.
void insertionSortRowsMatrixByRowCriteria(matrix m, int (*criteria)(int *, int));

// выполняет сортировку выбором столбцов матрицы m по неубыванию значения функции criteria применяемой для столбцов.
void selectionSortColsMatrixByColCriteria(matrix m, int (*criteria)(int *, int));

//возвращает значение ’истина’, если матрица m является квадратной, ложь – в противном случае.
bool isSquareMatrix(matrix *m);

// возвращает значение ’истина’, если матрицы m1 и m2 равны, ложь – в противном случае.
bool areTwoMatricesEqual(matrix *m1, matrix *m2);

// возвращает значение ’истина’, если матрица m является единичной, ложь – в противном случае.
bool isEMatrix(matrix *m);

//возвращает значение ’истина’, если матрица m является симметричной, ложь – в противном случае.
bool isSymmetricMatrix(matrix *m);

// транспонирует квадратную матрицу m.
void transposeSquareMatrix(matrix *m);

//  транспонирует матрицу m.
void transposeMatrix(matrix *m);

//  возвращает позицию минимального элемента матрицы m.
position getMinValuePos(matrix m);

// возвращает позицию максимального элемента матрицы m.
position getMaxValuePos(matrix m);

//возвращает матрицу размера nRows на nCols, построенную из элементов массива a.
matrix createMatrixFromArray(const int *a, size_t nRows, size_t nCols);

// возвращает указатель на нулевую матрицу массива из nMatrices матриц, размещенных
// в динамической памяти, построенных из элементов массива a.
matrix *createArrayOfMatrixFromArray(const int *values, size_t nMatrices, size_t nRows, size_t nCols);

//возвращает наибольший элемент в массиве.
int getMax(int *a, int n);

//возвращает наименьший элементв в массиве.
int getMin(int *a, int n);

// возводит матрицу в квадрат.
void squareMatrix(matrix *m);

// возвращает сумму массива.
int sumArray(int arr[], int size);

// проверяет массив на уникальность.
bool isUnique(const int *array, size_t size);

// умножает две матрицы.
matrix squareMatrixMultiplication(matrix *m1, matrix *m2);

// выполняет сортировку вставками строк матрицы m по неубыванию значения функции criteria применяемой для строк, для
// типа float.
void insertionSortRowsMatrixByRowCriteriaF(matrix m, float (*criteria)(int *, int));

// возвращает long long сумму.
long long sum(const int *array, int size);

// проверяет, является ли массив a неубывающей последовательностью.
bool isNonDescendingSorted(const int *a, int n);

//проверяет, содержит ли каждая строка матрицы m неубывающую последовательность.
bool hasAllNonDescendingRows(matrix m);

//считает количество вхождений заданного значения value в массив a длиной n.
int countValues(const int *a, int n, int value);

// подсчитывает количество строк в матрице m, где все элементы строки равны 0.
int countZeroRows(matrix m);

//вычисляет норму матрицы, которая представляет собой наибольшее абсолютное значение среди всех элементов матрицы.
int getMatrixNorm(matrix m);

#endif