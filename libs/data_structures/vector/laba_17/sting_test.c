#include "libs/data_structures/laba_17/string_.h"

void test_findNonSpace() {
    char c[] = " a b c";
    char *res = findNonSpace(c);
    assert(*res == 'a');
    assert(res == &c[1]);
}

void test_findSpace() {
    char c[] = "1234 abcd";
    char *res = findSpace(c);
    assert(*res == ' ');
    assert(res == &c[4]);
}

void test_findNonSpaceReverse() {
    char c[] = "cool ";
    char *res = findNonSpaceReverse(&c[4], &c[-1]);
    assert(*res == 'l');
    assert(res == &c[3]);
}

void test_findSpaceReverse() {
    char c[] = "i love programming";
    char *res = findSpaceReverse(&c[17], &c[-1]);
    assert(*res == ' ');
    assert(res == &c[6]);
}

void test_strcmp_equal() {
    char c1[] = "abc";
    char c2[] = "abc";
    assert(strcmp(c1, c2) == 0);
}

void test_strcmp_lesser() {
    char c1[] = "abc";
    char c2[] = "ccd";
    assert(strcmp(c1, c2) == -2);
}

void test_strcmp_bigger() {
    char c1[] = "bcd";
    char c2[] = "abc";
    assert(strcmp(c1, c2) == 1);
}


void test_strcmp() {
    test_strcmp_equal();
    test_strcmp_lesser();
    test_strcmp_bigger();
}

void test_copy() {
    char c1[] = "mongol";
    char c2[7];
    char *cc = copy(c1, &c1[6], c2);
    assert(c2[0] == 'm');
    assert(c2[1] == 'o');
    assert(c2[2] == 'n');
    assert(c2[3] == 'g');
    assert(c2[4] == 'o');
    assert(c2[5] == 'l');
    assert(cc == &c2[6] && c2[6] != '\0');
}

void test_copyIf() {
    char c1[] = "big brain";
    char c2[10] = "1234";
    char *c3 = copyIf(c1, &c1[9], c2 + 4, isspace);
    assert(c2[4] == ' ');
    assert(c3 == &c2[5]);
}

bool isNumber(char letter) {
    return letter >= '0' && letter <= '9';
}

void test_copyIfReverse() {
    char c1[] = "1a2b3c";
    char c2[10];
    char *c3 = copyIfReverse(&c1[6], c1, c2, isNumber);
    assert(c2[0] == '3');
    assert(c3 = &c2[3]);
}

void test_string() {
    test_findNonSpace();
    test_findSpace();
    test_findNonSpaceReverse();
    test_findSpaceReverse();
    test_strcmp();
    test_copy();
    test_copyIf();
    test_copyIfReverse();
}

int main() {
    test_string();

    return 0;
}

