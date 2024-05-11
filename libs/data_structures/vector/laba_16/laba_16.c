// № 1
#include "libs/data_structures/matrix/matrix.h"
#include <stdio.h>

int main(){
    int nRows, nCols;
    scanf("%d %d", &nRows, &nCols);

    matrix m = getMemMatrix(nRows, nCols);
    inputMatrix(m);

    position p1 = getMinValuePos(m);
    position p2 = getMaxValuePos(m);
    swapRows(m, p1.rowIndex, p2.rowIndex);

    outputMatrix(m);

    return 0;
}

//№2
getMax(int *a, int n) {
    int max = a[0];
    for (int i = 1; i < n; i++)
        if (a[i] > max)
            max = a[i];
    return max;
}

void sortRowsByMinElement(matrix m) {
    insertionSortRowsMatrixByRowCriteria(m, getMax);
}

int main() {
    int nRows, nCols;
    scanf("%d %d", &nRows, &nCols);

    matrix m = getMemMatrix(nRows, nCols);
    inputMatrix(m);

    sortRowsByMinElement(m);

    outputMatrix(m);

    return 0;
}

// №3
#include "libs/data_structures/matrix/matrix.h"
int getMin(int *a, int n) {
    int min = a[0];
    for (int i = 1; i < n; i++)
        if (a[i] < min)
            min = a[i];
    return min;
}
void sortRowsByMinElement(matrix m) {
    insertionSortColsMatrixByColCriteria(m, getMin);
}
int main() {
    int nRows, nCols;
    scanf("%d %d", &nRows, &nCols);

    matrix m = getMemMatrix(nRows, nCols);
    inputMatrix(m);

    sortRowsByMinElement(m);

    outputMatrix(m);

    return 0;
}


//№4
#include "libs/data_structures/matrix/matrix.h"
matrix mulMatrices(matrix m1, matrix m2) {
    assert(m1.nCols == m2.nRows);
    matrix m3 = getMemMatrix(m1.nRows, m2.nCols);
    for (int i = 0; i < m3.nRows; i++)
        for (int j = 0; j < m3.nCols; j++) {
            m3.values[i][j] = 0;
            for (size_t k = 0; k < m1.nCols; k++)
                m3.values[i][j] += m1.values[i][k] * m2.values[k][j];
        }
    return m3;
}

void getSquareOfMatrixIfSymmetric(matrix *m) {
    if (isSymmetricMatrix(*m))
        *m = mulMatrices(*m, *m);
}
int main() {
    int nRows, nCols;
    scanf("%d %d", &nRows, &nCols);

    matrix m = getMemMatrix(nRows, nCols);
    inputMatrix(m);

    getSquareOfMatrixIfSymmetric(&m);

    outputMatrix(m);

    return 0;
}

// №5
#include "libs/data_structures/matrix/matrix.h"

bool isUnique(int *a, int n) {
    for (size_t i = 0; i < n - 1; i++) {
        for (size_t j = i + 1; j < n; j++) {
            if (a[i] == a[j])
                return false;
        }
    }
    return true;
}

void transposeIfMatrixHasNotEqualSumOfRows(matrix m) {
    long long sumArray[m.nRows];
    for (int i = 0; i < m.nRows; i++)
        sumArray[i] = getSum(m.values[i], m.nCols);
    if (isUnique(sumArray, m.nRows))
        transposeSquareMatrix(m);
}

int main(){
    int nRows, nCols;
    scanf("%d %d", &nRows, &nCols);

    matrix m = getMemMatrix(nRows, nCols);
    inputMatrix(m);

    transposeIfMatrixHasNotEqualSumOfRows(m);

    outputMatrix(m);

    return 0;
}

//№6
#include "libs/data_structures/matrix/matrix.h"

matrix mulMatrices(matrix m1, matrix m2) {
    assert(m1.nCols == m2.nRows);
    matrix m3 = getMemMatrix(m1.nRows, m2.nCols);
    for (int i = 0; i < m3.nRows; i++)
        for (int j = 0; j < m3.nCols; j++) {
            m3.values[i][j] = 0;
            for (size_t k = 0; k < m1.nCols; k++)
                m3.values[i][j] += m1.values[i][k] * m2.values[k][j];
        }
    return m3;
}

bool isMutuallyInverseMatrices(matrix m1, matrix m2) {
    matrix m3 = mulMatrices(m1, m2);
    return isEMatrix(m3) ? 1 : 0;
}

int main() {
    int nRows1, nCols1;
    scanf("%d %d", &nRows1, &nCols1);
    matrix m1 = getMemMatrix(nRows1, nCols1);
    inputMatrix(m1);
    int nRows2, nCols2;
    scanf("%d %d", &nRows2, &nCols2);
    matrix m2 = getMemMatrix(nRows2, nCols2);
    inputMatrix(m2);
    printf("%d", isMutuallyInverseMatrices(m1, m2));
    return 0;
}

// №7
#include "libs/data_structures/matrix/matrix.h"

int getSum(int *a, int n) {
    int sum = 0;
    for (int i = 0; i < n; i++)
        sum += a[i];
    return sum;
}

int max(int a, int b) {
    return a > b ? a : b;
}

long long findSumOfMaxesOfPseudoDiagonal(matrix m) {
    int n = m.nRows + m.nCols - 2;
    int maxesOfPseudoDiagonal[n];
    for (size_t i = 0; i < n; i++)
        maxesOfPseudoDiagonal[i] = 0;
    int indexPseudoDiagonalElement;
    for (int i = 0; i < m.nRows; i++)
        for (int j = 0; j < m.nCols; j++)
            if (j != i) {
                if (i > j)
                    indexPseudoDiagonalElement = j - i + m.nRows - 1;
                else
                    indexPseudoDiagonalElement = j - i + m.nRows - 2;
                maxesOfPseudoDiagonal[indexPseudoDiagonalElement] =
                        max(maxesOfPseudoDiagonal[indexPseudoDiagonalElement],
                            m.values[i][j]);
            }
    return getSum(maxesOfPseudoDiagonal, n);
}

int main() {
    int nRows1, nCols1;
    scanf("%d %d", &nRows1, &nCols1);
    matrix m = getMemMatrix(nRows1, nCols1);
    inputMatrix(m);
    printf("%lld", findSumOfMaxesOfPseudoDiagonal(m));
    return 0;
}

// №8
#include "libs/data_structures/matrix/matrix.h"

int getMin(int *a, int left, int right) {
    \
 int min = a[left];
    for (int i = left + 1; i <= right; i++) {
        if (a[i] < min)
            min = a[i];
    }
    return min;
}

int getMinInArea(matrix m) {
    int minInArea;
    position maxElement = getMaxValuePos(m);
    int left = maxElement.colIndex;
    int right = maxElement.colIndex;
    for (int i = maxElement.rowIndex; i >= 0; i--) {
        int minInRow = getMin(m.values[i], left, right);
        if (i == maxElement.rowIndex || minInRow < minInArea)
            minInArea = minInRow;
        if (right < m.nCols - 1)
            right++;
        if (left > 0)
            left--;
    }
    return minInArea;
}

int main() {
    int nRows1, nCols1;
    scanf("%d %d", &nRows1, &nCols1);
    matrix m1 = getMemMatrix(nRows1, nCols1);
    inputMatrix(m1);
    printf("%d", getMinInArea(m1));
    return 0;
}

// №9
#include "libs/data_structures/matrix/matrix.h"
#include <math.h>

float getDistance(int *a, int n) {
    float distance = 0;
    for (int i = 0; i < n; i++)
        distance += a[i] * a[i];
    return sqrt(distance);
}

void insertionSortRowsMatrixByRowCriteriaF(matrix m, float (*criteria)(int *,
                                                                       int)) {
    float criteriaArray[m.nRows];
    for (int i = 0; i < m.nRows; i++)
        criteriaArray[i] = criteria(m.values[i], m.nCols);
    for (int i = 1; i < m.nRows; i++) {
        float t = criteriaArray[i];
        int j = i;
        while (j > 0 && criteriaArray[j - 1] > t) {
            criteriaArray[j] = criteriaArray[j - 1];
            swapRows(m, j, j - 1);
            j--;
        }
        criteriaArray[j] = t;
    }
}

void sortByDistances(matrix m) {
    insertionSortRowsMatrixByRowCriteriaF(m, getDistance);
}

int main() {
    int nRows1, nCols1;
    scanf("%d %d", &nRows1, &nCols1);
    matrix m1 = getMemMatrix(nRows1, nCols1);
    inputMatrix(m1);
    sortByDistances(m1);
    outputMatrix(m1);
    return 0;
}

// №10
#include "libs/data_structures/matrix/matrix.h"
#include <math.h>

bool isNonDescendingSorted(int *a, int n) {
    for (int i = 1; i < n; i++)
        if (a[i] < a[i - 1])
            return false;
    return true;
}

bool hasAllNonDescendingRows(matrix m) {
    for (int i = 0; i < m.nRows; i++)
        if (!isNonDescendingSorted(m.values[i], m.nCols))
            return false;
    return true;
}

int countNonDescendingRowsMatrices(matrix *ms, int nMatrix) {
    int countMatrix = 0;
    for (int i = 0; i < nMatrix; i++)
        if (hasAllNonDescendingRows(ms[i]))
            countMatrix += 1;
    return countMatrix;
}

int main() {
    int nRows1, nCols1;
    scanf("%d %d", &nRows1, &nCols1);
    int nMatrices;
    scanf("%d", &nMatrices);
    matrix *ms = getMemArrayOfMatrices(nMatrices, nRows1, nCols1);
    inputMatrices(ms, nMatrices);
    printf("%d", countNonDescendingRowsMatrices(ms, nMatrices));
    return 0;
}

// №11
#include "libs/data_structures/matrix/matrix.h"
#include <math.h>

position getLeftMin(matrix m) {
    int min = m.values[0][0];
    int minRow = 0;
    int minCol = 0;
    for (int i = 0; i < m.nCols; i++)
        for (int j = 0; j < m.nRows; j++)
            if (m.values[j][i] < min) {
                min = m.values[j][i];
                minRow = j;
                minCol = i;
            }
    return (position) {minRow, minCol};
}

void swapPenultimateRow(matrix m) {
    assert(m.nRows > 1);
    position posMin = getLeftMin(m);
    int minColEl[m.nRows];
    for (int i = 0; i < m.nRows; ++i)
        minColEl[i] = m.values[i][posMin.colIndex];
    for (size_t i = 0; i < m.nRows; i++)
        m.values[m.nRows - 2][i] = minColEl[i];
}

int main() {
    int nRows1, nCols1;
    scanf("%d %d", &nRows1, &nCols1);
    matrix m1 = getMemMatrix(nRows1, nCols1);
    inputMatrix(m1);
    swapPenultimateRow(m1);
    outputMatrix(m1);
    return 0;
}

