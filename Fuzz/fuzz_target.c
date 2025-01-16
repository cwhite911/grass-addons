#include <stdint.h>
#include <stddef.h>
#include "pavl.h"
#include "sort.h"
#include "local_proto.h"
#include "test_simulation.cpp"
#include "test_treatments.cpp"

extern "C" int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size)
{
    if (size < sizeof(int))
        return 0;

    int input = *(int *)data;

    // fuzz test for pavl_create function from pavl.h
    struct pavl_table *table = pavl_create(NULL, NULL, NULL);
    if (table) {
        pavl_destroy(table, NULL);
    }

    // fuzz test for Mergesort_increasing_smallest_azimuth function from sort.h
    struct node head;
    struct node tail;
    Mergesort_increasing_smallest_azimuth(&head, &tail);

    // Fuzz test for functions from local_proto.h
    float result = f_and(input, input, E);
    result += f_or(input, input, E);
    result += f_not(input, E);

    // Fuzz test for test_simulation.cpp
    ret += test_calling_all_functions();

    // Fuzz test for test_treatments.cpp
    ret += test_steering();

    return 0;
}
