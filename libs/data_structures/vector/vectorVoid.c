#include "vectorVoid.h"

vectorVoid createVectorV(size_t n, size_t baseTypeSize) {
    void *data = malloc(baseTypeSize * n);
    if (data == NULL) {
        fprintf(stderr, "bad alloc");
        exit(1);
    }
    return (vectorVoid) {data, 0, n, baseTypeSize};
}

void reserveV(vectorVoid *v, size_t newCapacity) {
    if (newCapacity == 0) {
        free(v->data);
    } else {
        void *newData = realloc(v->data, newCapacity * v->baseTypeSize);
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

void shrinkToFitV(vectorVoid *v) {
    reserveV(v, v->size);
}

void clearV(vectorVoid *v) {
    v->size = 0;
}

void deleteVectorV(vectorVoid *v) {
    free(v->data);
}

bool isEmptyV(vectorVoid *v) {
    return v->size == 0;
}

bool isFullV(vectorVoid *v) {
    return v->size == v->capacity;
}

void getVectorValueV(vectorVoid *v, size_t index, void *destination) {
    // Проверка на нулевые указатели
    if (v == NULL || destination == NULL) {
        fprintf(stderr, "Error: NULL pointer argument detected\n");
        exit(1);
    }
    // Проверка на выход за пределы размера вектора
    if (index >= v->size) {
        fprintf(stderr, "Error: Index %zu is out of bounds (vector size: %zu)\n", index, v->size);
        exit(1);
    }
    // Получение адреса начала данных для index-ого элемента вектора
    char *source = (char *) v->data + index * v->baseTypeSize;
    // Копирование данных из вектора в destination
    memcpy(destination, source, v->baseTypeSize);
}


void setVectorValueV(vectorVoid *v, size_t index, void *source) {
    // Проверка на нулевые указатели
    if (v == NULL || source == NULL) {
        fprintf(stderr, "Error: NULL pointer argument detected\n");
        exit(1);
    }
    // Проверка на выход за пределы размера вектора
    if (index >= v->size) {
        fprintf(stderr, "Error: Index %zu is out of bounds (vector size: %zu)\n", index, v->size);
        exit(1);
    }
    // Вычисление адреса для index-ого элемента вектора и копирование данных
    char *destination = (char *) (v->data + index * v->baseTypeSize);
    memcpy(destination, source, v->baseTypeSize);
}


void popBackV(vectorVoid *v) {
    if (isEmptyV(v)) {
        fprintf(stderr, "vector is empty");
        exit(1);
    } else {
        v->size--;
        shrinkToFitV(v);
    }
}

void pushBackV(vectorVoid *v, void *source) {
    if (isFullV(v)) {
        reserveV(v, v->capacity == 0 ? 1 : v->capacity * 2);
    }
    char *destination = (char *) (v->data + v->size++ * v->baseTypeSize);
    memcpy(destination, source, v->baseTypeSize);
}
