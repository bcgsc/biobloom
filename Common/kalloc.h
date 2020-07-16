#ifndef KALLOC_H
#define KALLOC_H

#define kmalloc(km, size) (malloc(size))
#define kcalloc(km, count, size) (calloc(count, size))
#define krealloc(km, ptr, size) (realloc(ptr, size))
#define kfree(km, ptr) (free(ptr))

#endif
