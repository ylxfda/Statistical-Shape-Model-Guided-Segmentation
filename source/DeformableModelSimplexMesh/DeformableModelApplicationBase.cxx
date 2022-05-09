


#include "DeformableModelApplicationBase.h"
#include "Surface.h"
extern Surface* Sur;

DeformableModelApplicationBase
::DeformableModelApplicationBase()
{
  m_VolumeReader      = VolumeReaderType::New();

  m_CastImage         = CastImageType::New();

  m_GradientAnisotropicImage  = GradientAnisotropicImageType::New();

  m_GradientMagnitude = GradientMagnitudeType::New();

  m_SigmoidImage = SigmoidImageType::New();

  m_RescaleIntensity  = RescaleIntensityFilterType::New();

  m_ITK2VTKAdaptor    = ITK2VTKAdaptorFilterType::New();
  m_ITK2VTKAdaptor->SetInput( m_RescaleIntensity->GetOutput() );

  m_RescaleIntensity->SetOutputMaximum( itk::NumericTraits< VisualizationPixelType >::max() );
  m_RescaleIntensity->SetOutputMinimum( itk::NumericTraits< VisualizationPixelType >::min() );

  m_SphereMeshSource        = SphereMeshSourceType::New();
  PointType center; 
  center.Fill(0);
 
  VectorType sphereRadius;
  sphereRadius.Fill( 10.0 );

  m_SphereMeshSource->SetCenter(center);
  m_SphereMeshSource->SetScale( sphereRadius );
  m_SphereMeshSource->SetResolution(2); 

  m_VTKImageExport = VTKImageExportType::New();

  m_ImageToVTKImage = ImageToVTKImageType::New();

  m_SimplexFilter  = SimplexFilterType::New();
  m_SimplexFilter->SetInput( m_SphereMeshSource->GetOutput() );

  //deformation stuff
  m_DeformFilter = DeformFilterType::New();

  m_TriangleMesh = TriangleMeshType::New();

  /////////////////////////added by me///////////////////////////////////////////////////////////
  m_TriangleMeshSource = TriangleMeshSourceType::New();
  m_SurfaceCutter = vtkCutter::New();
  m_cutStripter = vtkStripper::New();
  m_cutPlaneAxial = vtkPlane::New();
  m_cutPlaneSagittal = vtkPlane::New();
  m_cutPlaneCoronal = vtkPlane::New();
  m_vtkMesh = vtkPolyData::New();
  m_vtkContour = vtkPolyData::New();
  m_PreprocessedReader = VolumeReaderType::New();

  /* Level0StepLen=10;
   Level0StepNum=5;
   Level0Iter=3;
   Level1StepLen=10;
   Level1StepNum=10;
   Level1Iter=2;
   Level2StepLen=5;
   Level2StepNum=3;
   Level2Iter=1;
   Level3StepLen=3;
   Level3StepNum=2;
   Level3Iter=1;
   Level4StepLen=15;
   Level4StepNum=3;
   Level4Iter=0;*/

   Level0StepLen=30;
   //Level0StepNum=5;
   Level0StepNum=6;
   Level0Iter=3;
   Level1StepLen=30;
   Level1StepNum=6;
   Level1Iter=2;
   Level2StepLen=30;
   Level2StepNum=2;
   Level2Iter=2;
   Level3StepLen=25;
   Level3StepNum=2;
   Level3Iter=1;
   Level4StepLen=15;
   Level4StepNum=3;
   Level4Iter=0;
  ////////////////////////////////////////////////////////////////////////////////////////////////

  m_SimplexMeshFilter  = SimplexFilterType::New();

  m_GradientFilter = GradientFilterType::New();

  m_IterationObserver = IterationObserverType::New();
  m_IterationObserver->SetCallbackFunction( this, & DeformableModelApplicationBase::IterationCallback );

  m_DeformFilter->AddObserver( itk::IterationEvent(), m_IterationObserver );
  m_DeformFilter->AddObserver( itk::ProgressEvent(), m_IterationObserver );

  m_SimplexToTriangle = TriangleFilterType::New();

  m_TriangleToImage = TriangleMeshToBinaryImageFilterType::New();
  
  m_ImageWriter = ImageWriterType::New(); 

}


DeformableModelApplicationBase
::~DeformableModelApplicationBase()
{
	
}


void  
DeformableModelApplicationBase::SetSeedPoint( double x, double y, double z )
{
  m_SeedPoint[0] = x;
  m_SeedPoint[1] = y;
  m_SeedPoint[2] = z;
 
}


void 
DeformableModelApplicationBase::GetSeedPoint(double data[3])
{
  for(int i=0; i<3; i++)
    {
      data[i] = m_SeedPoint[i];
    }
}


void 
DeformableModelApplicationBase::IterationCallback()
{
  std::cout << "Iteration " << m_DeformFilter->GetProgress() << std::endl;
}

