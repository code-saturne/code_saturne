#ifndef CALLBACKS_H
#define CALLBACKS_H

#include "fmi/fmi2Functions.h"
#include <functional>

void setFunctions(const fmi2CallbackFunctions *functions);

void* fmiAlloc(size_t size);

void* fmiAlloc(size_t size, size_t size2);

void fmiFree(void* pointer);

void log(const char* id, const char* category, const char* text);

void log(fmi2ComponentEnvironment, fmi2String, fmi2Status, fmi2String, fmi2String, void* pointer);

#endif // CALLBACKS_H
