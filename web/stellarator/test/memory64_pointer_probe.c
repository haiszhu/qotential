#include <stdint.h>

void *pointer_identity(void *pointer)
{
    return pointer;
}

void pointer_store(uint8_t *pointer, uint8_t value)
{
    *pointer = value;
}
