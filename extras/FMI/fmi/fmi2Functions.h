#ifndef fmi2Functions_h
#define fmi2Functions_h

#ifdef __cplusplus
extern "C" {
#endif

#include "fmi2TypesPlatform.h"
#include "fmi2FunctionTypes.h"
#include <stdlib.h>

#if defined(_WIN32) || defined(__CYGWIN__)
	#define EXPORT __declspec(dllexport)
#else
	#define EXPORT __attribute__ ((visibility ("default")))
#endif

	EXPORT const char* fmi2GetTypesPlatform(void);
	EXPORT const char* fmi2GetVersion(void);

	EXPORT fmi2Component fmi2Instantiate(fmi2String instanceName,fmi2Type fmuType, fmi2String fmuGUID, fmi2String fmuResourceLocation, const fmi2CallbackFunctions* callbacks, fmi2Boolean visible, fmi2Boolean loggingOn);
	EXPORT void fmi2FreeInstance(fmi2Component component);
	EXPORT fmi2Status fmi2SetDebugLogging(fmi2Component component, fmi2Boolean loggingOn, size_t nCategories, const fmi2String categories[]);

	EXPORT fmi2Status fmi2SetupExperiment(fmi2Component component, fmi2Boolean toleranceDefined, fmi2Real tolerance, fmi2Real startTime, fmi2Boolean stopTimeDefined, fmi2Real stopTime);
	EXPORT fmi2Status fmi2EnterInitializationMode(fmi2Component component);
	EXPORT fmi2Status fmi2ExitInitializationMode(fmi2Component component);
	EXPORT fmi2Status fmi2Terminate(fmi2Component component);
	EXPORT fmi2Status fmi2Reset(fmi2Component component);

	EXPORT fmi2Status fmi2GetReal(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, fmi2Real values[]);
	EXPORT fmi2Status fmi2GetInteger(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, fmi2Integer values[]);
	EXPORT fmi2Status fmi2GetBoolean(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, fmi2Boolean values[]);
	EXPORT fmi2Status fmi2GetString(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, fmi2String values[]);
	EXPORT fmi2Status fmi2SetReal(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, const fmi2Real values[]);
	EXPORT fmi2Status fmi2SetInteger(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, const fmi2Integer values[]);
	EXPORT fmi2Status fmi2SetBoolean(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, const fmi2Boolean values[]);
	EXPORT fmi2Status fmi2SetString(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, const fmi2String values[]);

	EXPORT fmi2Status fmi2GetFMUstate(fmi2Component component, fmi2FMUstate* state);
	EXPORT fmi2Status fmi2SetFMUstate(fmi2Component component, fmi2FMUstate  state);
	EXPORT fmi2Status fmi2FreeFMUstate(fmi2Component component, fmi2FMUstate* state);

	EXPORT fmi2Status fmi2SerializedFMUstateSize(fmi2Component component, fmi2FMUstate state, size_t* size);
	EXPORT fmi2Status fmi2SerializeFMUstate(fmi2Component component, fmi2FMUstate state, fmi2Byte serializedState[], size_t size);
	EXPORT fmi2Status fmi2DeSerializeFMUstate(fmi2Component component, const fmi2Byte serializedState, size_t size, fmi2FMUstate* state);

	EXPORT fmi2Status fmi2GetDirectionalDerivative(fmi2Component component, const fmi2ValueReference unknownValueReferences[], size_t numberOfUnknowns, const fmi2ValueReference knownValueReferences[], fmi2Integer numberOfKnowns, fmi2Real knownDifferential[], fmi2Real unknownDifferential[]);

	EXPORT fmi2Status fmi2SetRealInputDerivatives(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, fmi2Integer orders[], const fmi2Real values[]);
	EXPORT fmi2Status fmi2GetRealOutputDerivatives(fmi2Component component, const fmi2ValueReference valueReferences[], size_t numberOfValueReferences, const fmi2Integer orders[], fmi2Real values[]);

	EXPORT fmi2Status fmi2DoStep(fmi2Component component, fmi2Real currentCommunicationPoint, fmi2Real communicationStepSize, fmi2Boolean noSetFMUStatePriorToCurrentPoint);
	EXPORT fmi2Status fmi2CancelStep(fmi2Component component);

	EXPORT fmi2Status fmi2GetStatus(fmi2Component component, const fmi2StatusKind statusKind, fmi2Status* status);
	EXPORT fmi2Status fmi2GetRealStatus(fmi2Component component, const fmi2StatusKind statusKind, fmi2Real* status);
	EXPORT fmi2Status fmi2GetIntegerStatus(fmi2Component component, const fmi2StatusKind statusKind, fmi2Integer* status);
	EXPORT fmi2Status fmi2GetBooleanStatus(fmi2Component component, const fmi2StatusKind statusKind, fmi2Boolean* status);
	EXPORT fmi2Status fmi2GetStringStatus(fmi2Component component, const fmi2StatusKind statusKind, fmi2String* status);

#ifdef __cplusplus
}
#endif

#endif /* fmi2Functions_h */
