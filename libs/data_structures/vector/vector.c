#include "vector.h"

vector createVector(size_t n) {
    int *data = malloc(sizeof(int) * n);
    if (data == NULL) {
        fprintf(stderr, "bad alloc");
        exit(1);
    }
    return (vector) {data, 0, n};
}

void reserve(vector *v, size_t newCapacity) {
    if (newCapacity == 0) {
        free(v->data);
    } else {
        int *newData = realloc(v->data, newCapacity * sizeof(int));
        if (newData == NULL) {
            fprintf(stderr, "bad alloc");
            exit(1);
        } else {
            if (newCapacity <= v->size)
                v->size = newCapacity;
            v->data = newData;
            v->capacity = newCapacity;
        }
    }
}


void clear(vector *v) {
    v->size = 0;
}

void shrinkToFit(vector *v) {
    reserve(v, v->size);
}

void deleteVector(vector *v) {
    free(v->data);
}

bool isEmpty(vector *v) {
    return v->size == 0;
}

bool isFull(vector *v) {
    return v->size == v->capacity;
}

int getVectorValue(vector *v, size_t i) {
    return i < v->size ? v->data[i] : -1;
}

void pushBack(vector *v, int x) {
    if (isFull(v)) {
        reserve(v, v->capacity == 0 ? 1 : v->capacity * 2);
    }
    v->data[v->size++] = x;
}

void popBack(vector *v) {
    if (isEmpty(v)) {
        fprintf(stderr, "vector is empty");
        exit(1);
    } else {
        v->size--;
        shrinkToFit(v);
    }
}

int *atVector(vector *v, size_t index) {
    if (index >= v->size) {
        fprintf(stderr, "IndexError: a[%d] is not exists", index);
        exit(1);
    }
    return &v->data[index];
}

int *back(vector *v) {
    return &v->data[v->size - 1];
}

int *front(vector *v) {
    return &v->data[0];
}