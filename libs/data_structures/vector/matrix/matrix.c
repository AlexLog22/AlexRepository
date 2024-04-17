#include "matrix.h"
#include "../../algorithms/array/array.h"

matrix getMemMatrix(int nRows, int nCols) {
    int **values = (int **) malloc(sizeof(int *) * nRows);
    for (int i = 0; i < nRows; i++)
        values[i] = (int *) malloc(sizeof(int) * nCols);
    return (matrix) {values, nRows, nCols};
}

matrix *getMemArrayOfMatrices(int nMatrices, int nRows, int nCols) {
    matrix *ms = (matrix *) malloc(sizeof(matrix) * nMatrices);
    for (int i = 0; i < nMatrices; i++)
        ms[i] = getMemMatrix(nRows, nCols);
    return ms;
}

void freeMemMatrix(matrix *m) {
    for (int i = 0; i < m->nRows; i++) {
        free(m->values[i]);
    }
    free(m->values);
    m->values = NULL;
    m->nRows = 0;
    m->nCols = 0;
}

void freeMemMatrices(matrix *ms, int nMatrices) {
    for (int i = 0; i < nMatrices; i++) {
        freeMemMatrix(&ms[i]);
    }
    free(ms);
    ms = NULL;
}

void inputMatrix(matrix *m) {
    for (int i = 0; i < m->nRows; i++) {
        for (int j = 0; j < m->nCols; j++) {
            scanf("%d", &(m->values[i][j]));
        }
    }
}

void inputMatrices(matrix *ms, int nMatrices) {
    for (int i = 0; i < nMatrices; i++)
        inputMatrix(&ms[i]);
}

void outputMatrix(matrix m) {
    int max_digits = findLongNumberInMatrix(m);

    for (int i = 0; i < m.nRows; i++) {
        printf("|");
        for (int j = 0; j < m.nCols; j++) {
            int num = m.values[i][j];

            int digits = countDigits(num);

            int spacesBefore = (max_digits - digits) / 2;

            int spacesAfter = max_digits - digits - spacesBefore;

            for (int k = 0; k < spacesBefore; k++) {
                printf(" ");
            }
            printf("%d", m.values[i][j]);

            for (int k = 0; k < spacesAfter; k++) {
                printf(" ");
            }
            printf(" ");
        }
        printf("|\n");
    }
    printf("\n");
}

int countDigits(int n) {
    int digits = 0;
    while (n != 0) {
        n /= 10;
        digits++;
    }
    return digits;
}

int findLongNumberInMatrix(matrix m) {
    int max_digits = 0;
    for (int i = 0; i < m.nRows; i++) {
        for (int j = 0; j < m.nCols; j++) {
            int num = m.values[i][j];
            int digits = 0;
            while (num != 0) {
                num /= 10;
                digits++;
            }
            if (digits > max_digits) {
                max_digits = digits;
            }
        }
    }
    return max_digits;
}

void outputMatrices(matrix *ms, int nMatrices) {
    for (int i = 0; i < nMatrices; i++)
        outputMatrix(ms[i]);
}

void swapRows(matrix m, int i1, int i2) {
    assert(i1 >= 0 && i1 < m.nRows && i2 >= 0 && i2 < m.nRows);
    int *temp = m.values[i1];
    m.values[i1] = m.values[i2];
    m.values[i2] = temp;
}

void swapColumns(matrix m, int j1, int j2) {
    assert(j1 >= 0 && j1 < m.nCols && j2 >= 0 && j2 < m.nCols);
    for (int i = 0; i < m.nRows; i++) {
        int temp = m.values[i][j1];
        m.values[i][j1] = m.values[i][j2];
        m.values[i][j2] = temp;
    }
}

void insertionSortRowsMatrixByRowCriteria(matrix m, int (*criteria)(int *, int)) {
    int *criteria_values = (int *) malloc(sizeof(int) * m.nRows);
    assert(criteria_values != NULL);
    for (int i = 0; i < m.nRows; i++) {
        criteria_values[i] = criteria(m.values[i], m.nCols);
    }
    for (int i = 1; i < m.nRows; i++) {
        int key_value = criteria_values[i];
        int j = i;
        while (j > 0 && criteria_values[j - 1] > key_value) {
            swapRows(m, j - 1, j);
            int temp = criteria_values[j - 1];
            criteria_values[j - 1] = criteria_values[j];
            criteria_values[j] = temp;
            j--;
        }
    }
    free(criteria_values);
}

void selectionSortColsMatrixByColCriteria(matrix m, int (*criteria)(int *, int)) {
    int *criteria_values = (int *) malloc(sizeof(int) * m.nCols);
    int *intermediate_array = (int *) malloc(sizeof(int) * m.nRows);
    assert(criteria_values != NULL && intermediate_array != NULL);

    for (int j = 0; j < m.nCols; j++) {
        for (int i = 0; i < m.nRows; i++) {
            intermediate_array[i] = m.values[i][j];
        }
        criteria_values[j] = criteria(intermediate_array, m.nRows);
    }
    for (int i = 0; i < m.nCols - 1; i++) {
        int min_index = i;
        for (int j = i + 1; j < m.nCols; j++) {
            if (criteria_values[j] < criteria_values[min_index]) {
                min_index = j;
            }
        }
        if (min_index != i) {
            int temp = criteria_values[i];
            criteria_values[i] = criteria_values[min_index];
            criteria_values[min_index] = temp;
            swapColumns(m, i, min_index);
        }
    }
    free(criteria_values);
    free(intermediate_array);
}

bool isSquareMatrix(matrix *m) {
    return (m->nRows == m->nCols);
}

bool areTwoMatricesEqual(matrix *m1, matrix *m2) {
    if (m1->nRows != m2->nRows || m1->nCols != m2->nCols) {
        return false;
    }
    int total_elements_m1 = m1->nRows * m1->nCols;
    int total_elements_m2 = m2->nRows * m2->nCols;
    int *values_m1 = (int *) malloc(total_elements_m1 * sizeof(int));
    int *values_m2 = (int *) malloc(total_elements_m2 * sizeof(int));
    if (values_m1 == NULL || values_m2 == NULL) {
        return false;
    }
    for (int i = 0; i < m1->nRows; i++) {
        memcpy(values_m1 + i * m1->nCols, m1->values[i], m1->nCols * sizeof(int));
        memcpy(values_m2 + i * m2->nCols, m2->values[i], m2->nCols * sizeof(int));
    }
    bool equal = memcmp(values_m1, values_m2, total_elements_m1 * sizeof(int)) == 0;
    free(values_m1);
    free(values_m2);

    return equal;
}

bool isEMatrix(matrix *m) {
    if (!isSquareMatrix(m))
        return false;
    for (int i = 0; i < m->nRows; i++) {
        for (int j = 0; j < m->nCols; j++) {
            if ((i == j && m->values[i][j] != 1) || (i != j && m->values[i][j] != 0)) {
                return false;
            }
        }
    }
    return true;
}

bool isSymmetricMatrix(matrix *m) {
    if (!isSquareMatrix(m))
        return false;
    for (int i = 0; i < m->nRows; i++) {
        for (int j = i + 1; j < m->nCols; j++) {
            if (m->values[i][j] != m->values[j][i])
                return false;
        }
    }
    return true;
}

void transposeSquareMatrix(matrix *m) {
    if (isSquareMatrix(m)) {
        for (int i = 0; i < m->nRows; i++) {
            for (int j = i + 1; j < m->nCols; j++) {
                int temp = m->values[i][j];
                m->values[i][j] = m->values[j][i];
                m->values[j][i] = temp;
            }
        }
    }
}

void transposeMatrix(matrix *m) {
    matrix transposed = getMemMatrix(m->nCols, m->nRows);
    for (int i = 0; i < m->nRows; ++i) {
        for (int j = 0; j < m->nCols; ++j) {
            transposed.values[j][i] = m->values[i][j];
        }
    }
    for (int i = 0; i < m->nRows; ++i) {
        free(m->values[i]);
    }
    free(m->values);
    m->nRows = transposed.nRows;
    m->nCols = transposed.nCols;
    m->values = transposed.values;
}

position getMinValuePos(matrix m) {
    int min_value = m.values[0][0];
    position min_pos = {0, 0};
    for (int i = 0; i < m.nRows; i++) {
        for (int j = 0; j < m.nCols; j++) {
            if (m.values[i][j] < min_value) {
                min_value = m.values[i][j];
                min_pos.rowIndex = i;
                min_pos.colIndex = j;
            }
        }
    }
    return min_pos;
}

position getMaxValuePos(matrix m) {
    int max_value = m.values[0][0];
    position max_pos = {0, 0};
    for (int i = 0; i < m.nRows; i++) {
        for (int j = 0; j < m.nCols; j++) {
            if (m.values[i][j] > max_value) {
                max_value = m.values[i][j];
                max_pos.rowIndex = i;
                max_pos.colIndex = j;
            }
        }
    }
    return max_pos;
}

matrix createMatrixFromArray(const int *a, size_t nRows, size_t nCols) {
    matrix m = getMemMatrix(nRows, nCols);
    int k = 0;
    for (int i = 0; i < nRows; i++)
        for (int j = 0; j < nCols; j++)
            m.values[i][j] = a[k++];
    return m;
}

matrix *createArrayOfMatrixFromArray(const int *values, size_t nMatrices, size_t nRows, size_t nCols) {
    matrix *ms = getMemArrayOfMatrices(nMatrices, nRows, nCols);
    int l = 0;
    for (size_t k = 0; k < nMatrices; k++)
        for (size_t i = 0; i < nRows; i++)
            for (size_t j = 0; j < nCols; j++)
                ms[k].values[i][j] = values[l++];
    return ms;
}

int getMax(int *a, int n) {
    int max = a[0];
    for (int i = 0; i < n; i++) {
        if (max < a[i])
            max = a[i];
    }
    return max;
}

int getMin(int *a, int n) {
    int min = a[0];
    for (int i = 0; i < n; i++) {
        if (min > a[i])
            min = a[i];
    }
    return min;
}

void squareMatrix(matrix *m) {
    matrix result = getMemMatrix(m->nRows, m->nCols);
    for (int i = 0; i < m->nRows; i++) {
        for (int j = 0; j < m->nCols; j++) {
            result.values[i][j] = 0;
            for (int k = 0; k < m->nCols; k++) {
                result.values[i][j] += m->values[i][k] * m->values[k][j];
            }
        }
    }
    for (int i = 0; i < m->nRows; i++) {
        for (int j = 0; j < m->nCols; j++) {
            m->values[i][j] = result.values[i][j];
        }
    }
    freeMemMatrix(&result);
}

int sumArray(int arr[], int size) {
    int sum = 0;
    for (int i = 0; i < size; i++) {
        sum += arr[i];
    }
    return sum;
}

bool isUnique(const int *array, size_t size) {
    for (size_t i = 0; i < size; ++i) {
        for (size_t j = i + 1; j < size; ++j) {
            if (array[i] == array[j]) {
                return false;
            }
        }
    }
    return true;
}

matrix squareMatrixMultiplication(matrix *m1, matrix *m2) {
    matrix result = getMemMatrix(m1->nRows, m1->nCols);
    for (int i = 0; i < m1->nRows; i++) {
        for (int j = 0; j < m1->nCols; j++) {
            result.values[i][j] = 0;
            for (int k = 0; k < m1->nCols; k++) {
                result.values[i][j] += m1->values[i][k] * m2->values[k][j];
            }
        }
    }
    return result;
}

void insertionSortRowsMatrixByRowCriteriaF(matrix m, float (*criteria)(int *, int)) {
    float *criteria_values = (float *) malloc(sizeof(float) * m.nRows);
    assert(criteria_values != NULL);
    for (int i = 0; i < m.nRows; i++) {
        criteria_values[i] = criteria(m.values[i], m.nCols);
    }
    for (int i = 1; i < m.nRows; i++) {
        float key_value = criteria_values[i];
        int j = i;
        while (j > 0 && criteria_values[j - 1] > key_value) {
            swapRows(m, j - 1, j);
            float temp = criteria_values[j - 1];
            criteria_values[j - 1] = criteria_values[j];
            criteria_values[j] = temp;
            j--;
        }
    }
    free(criteria_values);
}

long long sum(const int *array, int size) {
    long long sum = 0;

    for (int i = 0; i < size; ++i) {
        sum += array[i];
    }
    return sum;
}

bool isNonDescendingSorted(const int *a, int n) {
    for (int i = 0; i < n - 1; i++) {
        if (a[i + 1] < a[i]) {
            return false;
        }
    }
    return true;
}

bool hasAllNonDescendingRows(matrix m) {
    for (int i = 0; i < m.nRows; i++) {
        if (!isNonDescendingSorted(m.values[i], m.nCols)) {
            return false;
        }
    }
    return true;
}

int countValues(const int *a, int n, int value) {
    int amount = 0;
    for (int i = 0; i < n; i++) {
        if (a[i] == value) {
            amount++;
        }
    }
    return amount;
}

int countZeroRows(matrix m) {
    int amount = 0;
    for (int i = 0; i < m.nRows; i++) {
        if (countValues(m.values[i], m.nCols, 0) == m.nCols) {
            amount++;
        }
    }
    return amount;
}

int getMatrixNorm(matrix m) {
    int max_value = 0;
    for (int i = 0; i < m.nRows; i++) {
        int *row = m.values[i];
        for (int j = 0; j < m.nCols; j++) {
            max_value = max(max_value, abs(row[j]));
        }
    }
    return max_value;
}










