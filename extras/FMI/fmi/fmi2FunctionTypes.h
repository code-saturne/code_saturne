#ifndef fmi2FunctionTypes_h
#define fmi2FunctionTypes_h

#include "fmi2TypesPlatform.h"

#ifdef __cplusplus
extern "C" {
#endif

/* make sure all compiler use the same alignment policies for structures */
#if defined _MSC_VER || defined __GNUC__
#pragma pack(push,8)
#endif

/* Include stddef.h, in order that size_t etc. is defined */
#include <stddef.h>


/* Type definitions */
typedef enum {
    fmi2OK,
    fmi2Warning,
    fmi2Discard,
    fmi2Error,
    fmi2Fatal,
    fmi2Pending
} fmi2Status;

typedef enum {
    fmi2ModelExchange,
    fmi2CoSimulation
} fmi2Type;

typedef enum {
    fmi2DoStepStatus,
    fmi2PendingStatus,
    fmi2LastSuccessfulTime,
    fmi2Terminated
} fmi2StatusKind;

typedef void (*fmi2CallbackLogger)(fmi2ComponentEnvironment, fmi2String, fmi2Status, fmi2String, fmi2String, ...);
typedef void* (*fmi2CallbackAllocateMemory)(size_t, size_t);
typedef void (*fmi2CallbackFreeMemory)(void*);
typedef void (*fmi2StepFinished)(fmi2ComponentEnvironment, fmi2Status);

typedef struct {
   const fmi2CallbackLogger         logger;
   const fmi2CallbackAllocateMemory allocateMemory;
   const fmi2CallbackFreeMemory     freeMemory;
   const fmi2StepFinished           stepFinished;
   const fmi2ComponentEnvironment   componentEnvironment;
} fmi2CallbackFunctions;

typedef struct {
	 fmi2Boolean newDiscreteStatesNeeded;
   fmi2Boolean terminateSimulation;
   fmi2Boolean nominalsOfContinuousStatesChanged;
   fmi2Boolean valuesOfContinuousStatesChanged;
   fmi2Boolean nextEventTimeDefined;
   fmi2Real    nextEventTime;
} fmi2EventInfo;


/* reset alignment policy to the one set before reading this file */
#if defined _MSC_VER || defined __GNUC__
#pragma pack(pop)
#endif


/* Define fmi2 function pointer types to simplify dynamic loading */

/***************************************************
Types for Common Functions
****************************************************/

/* Inquire version numbers of header files and setting logging status */
   typedef const char* fmi2GetTypesPlatformTYPE(void);
   typedef const char* fmi2GetVersionTYPE(void);
   typedef fmi2Status  fmi2SetDebugLoggingTYPE(fmi2Component, fmi2Boolean, size_t, const fmi2String[]);

/* Creation and destruction of FMU instances and setting debug status */
   typedef fmi2Component fmi2InstantiateTYPE (fmi2String, fmi2Type, fmi2String, fmi2String, const fmi2CallbackFunctions*, fmi2Boolean, fmi2Boolean);
   typedef void          fmi2FreeInstanceTYPE(fmi2Component);

/* Enter and exit initialization mode, terminate and reset */
   typedef fmi2Status fmi2SetupExperimentTYPE        (fmi2Component, fmi2Boolean, fmi2Real, fmi2Real, fmi2Boolean, fmi2Real);
   typedef fmi2Status fmi2EnterInitializationModeTYPE(fmi2Component);
   typedef fmi2Status fmi2ExitInitializationModeTYPE (fmi2Component);
   typedef fmi2Status fmi2TerminateTYPE              (fmi2Component);
   typedef fmi2Status fmi2ResetTYPE                  (fmi2Component);

/* Getting and setting variable values */
   typedef fmi2Status fmi2GetRealTYPE   (fmi2Component, const fmi2ValueReference[], size_t, fmi2Real   []);
   typedef fmi2Status fmi2GetIntegerTYPE(fmi2Component, const fmi2ValueReference[], size_t, fmi2Integer[]);
   typedef fmi2Status fmi2GetBooleanTYPE(fmi2Component, const fmi2ValueReference[], size_t, fmi2Boolean[]);
   typedef fmi2Status fmi2GetStringTYPE (fmi2Component, const fmi2ValueReference[], size_t, fmi2String []);

   typedef fmi2Status fmi2SetRealTYPE   (fmi2Component, const fmi2ValueReference[], size_t, const fmi2Real   []);
   typedef fmi2Status fmi2SetIntegerTYPE(fmi2Component, const fmi2ValueReference[], size_t, const fmi2Integer[]);
   typedef fmi2Status fmi2SetBooleanTYPE(fmi2Component, const fmi2ValueReference[], size_t, const fmi2Boolean[]);
   typedef fmi2Status fmi2SetStringTYPE (fmi2Component, const fmi2ValueReference[], size_t, const fmi2String []);

/* Getting and setting the internal FMU state */
   typedef fmi2Status fmi2GetFMUstateTYPE           (fmi2Component, fmi2FMUstate*);
   typedef fmi2Status fmi2SetFMUstateTYPE           (fmi2Component, fmi2FMUstate);
   typedef fmi2Status fmi2FreeFMUstateTYPE          (fmi2Component, fmi2FMUstate*);
   typedef fmi2Status fmi2SerializedFMUstateSizeTYPE(fmi2Component, fmi2FMUstate, size_t*);
   typedef fmi2Status fmi2SerializeFMUstateTYPE     (fmi2Component, fmi2FMUstate, fmi2Byte[], size_t);
   typedef fmi2Status fmi2DeSerializeFMUstateTYPE   (fmi2Component, const fmi2Byte[], size_t, fmi2FMUstate*);

/* Getting partial derivatives */
   typedef fmi2Status fmi2GetDirectionalDerivativeTYPE(fmi2Component, const fmi2ValueReference[], size_t,
                                                                   const fmi2ValueReference[], size_t,
                                                                   const fmi2Real[], fmi2Real[]);

/***************************************************
Types for Functions for FMI2 for Model Exchange
****************************************************/

/* Enter and exit the different modes */
   typedef fmi2Status fmi2EnterEventModeTYPE         (fmi2Component);
   typedef fmi2Status fmi2NewDiscreteStatesTYPE      (fmi2Component, fmi2EventInfo*);
   typedef fmi2Status fmi2EnterContinuousTimeModeTYPE(fmi2Component);
   typedef fmi2Status fmi2CompletedIntegratorStepTYPE(fmi2Component, fmi2Boolean, fmi2Boolean*, fmi2Boolean*);

/* Providing independent variables and re-initialization of caching */
   typedef fmi2Status fmi2SetTimeTYPE            (fmi2Component, fmi2Real);
   typedef fmi2Status fmi2SetContinuousStatesTYPE(fmi2Component, const fmi2Real[], size_t);

/* Evaluation of the model equations */
   typedef fmi2Status fmi2GetDerivativesTYPE               (fmi2Component, fmi2Real[], size_t);
   typedef fmi2Status fmi2GetEventIndicatorsTYPE           (fmi2Component, fmi2Real[], size_t);
   typedef fmi2Status fmi2GetContinuousStatesTYPE          (fmi2Component, fmi2Real[], size_t);
   typedef fmi2Status fmi2GetNominalsOfContinuousStatesTYPE(fmi2Component, fmi2Real[], size_t);


/***************************************************
Types for Functions for FMI2 for Co-Simulation
****************************************************/

/* Simulating the slave */
   typedef fmi2Status fmi2SetRealInputDerivativesTYPE (fmi2Component, const fmi2ValueReference [], size_t, const fmi2Integer [], const fmi2Real []);
   typedef fmi2Status fmi2GetRealOutputDerivativesTYPE(fmi2Component, const fmi2ValueReference [], size_t, const fmi2Integer [], fmi2Real []);

   typedef fmi2Status fmi2DoStepTYPE     (fmi2Component, fmi2Real, fmi2Real, fmi2Boolean);
   typedef fmi2Status fmi2CancelStepTYPE (fmi2Component);

/* Inquire slave status */
   typedef fmi2Status fmi2GetStatusTYPE       (fmi2Component, const fmi2StatusKind, fmi2Status* );
   typedef fmi2Status fmi2GetRealStatusTYPE   (fmi2Component, const fmi2StatusKind, fmi2Real*   );
   typedef fmi2Status fmi2GetIntegerStatusTYPE(fmi2Component, const fmi2StatusKind, fmi2Integer*);
   typedef fmi2Status fmi2GetBooleanStatusTYPE(fmi2Component, const fmi2StatusKind, fmi2Boolean*);
   typedef fmi2Status fmi2GetStringStatusTYPE (fmi2Component, const fmi2StatusKind, fmi2String* );


#ifdef __cplusplus
}  /* end of extern "C" { */
#endif

#endif /* fmi2FunctionTypes_h */
