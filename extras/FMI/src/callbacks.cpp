#include <src/callbacks.h>

const fmi2CallbackFunctions *functions;

void setFunctions(const fmi2CallbackFunctions *newFunctions){
    functions = newFunctions;
}

void* fmiAlloc(size_t size){
    return functions->allocateMemory(size, 1);
}

void* fmiAlloc(size_t size, size_t size2){
    return functions->allocateMemory(size, size2);
}

void fmiFree(void* pointer){
    functions->freeMemory(pointer);
}

void log(const char* id, const char* category, const char* text){
    functions->logger(nullptr, id, fmi2OK, category, text);
}

void log(fmi2ComponentEnvironment environment, fmi2String id, fmi2Status status, fmi2String category, fmi2String text, void* pointers){
    functions->logger(environment, id, status, category, text, pointers);
}
