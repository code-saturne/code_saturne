#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkElevationFilter.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkXMLUnstructuredGridReader.h"

#include <mpi.h>
#include <string>

class DataGenerator {
public:
  DataGenerator()
    {
    this->Sphere = vtkSmartPointer<vtkSphereSource>::New();
    this->Sphere->SetThetaResolution(30);
    this->Sphere->SetPhiResolution(30);
    int procId;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    this->Sphere->SetCenter(procId*4.0, 0, 0);
    this->Elevation = vtkSmartPointer<vtkElevationFilter>::New();
    this->Elevation->SetInputConnection(this->Sphere->GetOutputPort());
    this->Index = 0;
    }

  vtkSmartPointer<vtkPolyData> GetNext()
    {
    double radius = fabs(sin(0.1 * this->Index));
    this->Index++;
    this->Sphere->SetRadius(1.0 + radius);
    this->Elevation->Update();
    vtkSmartPointer<vtkPolyData> ret = vtkSmartPointer<vtkPolyData>::New();
    ret->DeepCopy(this->Elevation->GetOutput());
    return ret;
    }

protected:
  int Index;
  vtkSmartPointer<vtkSphereSource> Sphere;
  vtkSmartPointer<vtkElevationFilter> Elevation;
};

int main(int argc, char* argv[])
{
  if (argc < 3)
    {
    printf("Usage: %s <python coprocessing script> <number of time steps>\n", argv[0]);
    return 1;
    }
  // we assume that this is done in parallel
  MPI_Init(&argc, &argv);

  std::string cpPythonFile = argv[1];
  int nSteps = atoi(argv[2]);

  vtkCPProcessor* processor = vtkCPProcessor::New();
  processor->Initialize();
  vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();

  // read the coprocessing python file
  if(pipeline->Initialize(cpPythonFile.c_str()) == 0)
    {
    cout << "Problem reading the python script.\n";
    return 1;
    }

  processor->AddPipeline(pipeline);
  pipeline->Delete();

  if (nSteps == 0)
    {
    return 0;
    }

  // create a data source, typically this will come from the adaptor
  // but here we use generator to create it ourselves
  DataGenerator generator;

  // do coprocessing
  double tStart = 0.0;
  double tEnd = 1.0;
  double stepSize = (tEnd - tStart)/nSteps;

  vtkCPDataDescription* dataDesc = vtkCPDataDescription::New();
  dataDesc->AddInput("input");

  for (int i = 0; i < nSteps; ++i)
    {
    double currentTime = tStart + stepSize*i;
    // set the current time and time step
    dataDesc->SetTimeData(currentTime, i);

    // check if the script says we should do coprocessing now
    if(processor->RequestDataDescription(dataDesc) != 0)
      {
      // we are going to do coprocessing so use generator to
      // create our grid at this timestep and provide it to
      // the coprocessing library
      vtkSmartPointer<vtkDataObject> dataObject =
        generator.GetNext();

      dataDesc->GetInputDescriptionByName("input")->SetGrid(dataObject);
      processor->CoProcess(dataDesc);
      }
    }

  dataDesc->Delete();
  processor->Finalize();
  processor->Delete();

  MPI_Finalize();

  return 0;
}


