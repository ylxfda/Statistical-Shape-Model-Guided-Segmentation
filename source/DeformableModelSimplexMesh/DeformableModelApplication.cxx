


#include "DeformableModelApplication.h"

#include "FL/Fl_File_Chooser.H"

#include "vtkImageShiftScale.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkImageMarchingCubes.h"
#include "vtkDecimatePro.h"

#include "vtkImageImport.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"
#include "vtkCellArray.h"

#include "ClickedPointEvent.h"
#include "itkImage.h"

#include "itkImageFileWriter.h"
#include "itkTimeProbe.h"

#ifndef vtkDoubleType
#define vtkDoubleType double
#endif

#ifndef vtkFloatingPointType
# define vtkFloatingPointType vtkFloatingPointType
typedef float vtkFloatingPointType;
#endif

/////  added by me////////////////////////////
#include "Surface.h"
extern Surface* Sur;
//#include"iostream.h"
#include <time.h>
//////////////////////////////////////////////

DeformableModelApplication
::DeformableModelApplication()
{
  m_AxialViewer.SetOrientation(    ImageSliceViewer::Axial    );
  m_CoronalViewer.SetOrientation(  ImageSliceViewer::Coronal  );
  m_SagittalViewer.SetOrientation( ImageSliceViewer::Sagittal );
  m_Dummy.SetOrientation(ImageSliceViewer::Axial);

  m_ShiftScaleImageFilter = vtkImageShiftScale::New();

  m_ShiftScaleImageFilter->SetInput(  m_ITK2VTKAdaptor->GetOutput() );
  m_ShiftScaleImageFilter->SetOutputScalarTypeToUnsignedChar();
  m_ShiftScaleImageFilter->ClampOverflowOn();

  m_AxialViewerCommand = itk::SimpleMemberCommand<DeformableModelApplication>::New();
  m_AxialViewerCommand->SetCallbackFunction(this, &DeformableModelApplication::ProcessAxialViewInteraction);
  m_AxialViewer.AddObserver(ClickedPointEvent(), m_AxialViewerCommand);

  m_CoronalViewerCommand = itk::SimpleMemberCommand<DeformableModelApplication>::New();
  m_CoronalViewerCommand->SetCallbackFunction(this, &DeformableModelApplication::ProcessCoronalViewInteraction);
  m_CoronalViewer.AddObserver(ClickedPointEvent(), m_CoronalViewerCommand);

  m_SagittalViewerCommand = itk::SimpleMemberCommand<DeformableModelApplication>::New();
  m_SagittalViewerCommand->SetCallbackFunction(this, &DeformableModelApplication::ProcessSagittalViewInteraction);
  m_SagittalViewer.AddObserver(ClickedPointEvent(), m_SagittalViewerCommand);

  
  //   m_SurfaceViewerCommand = itk::SimpleMemberCommand<DeformableModelApplication>::New();
  //  m_SurfaceViewerCommand->SetCallbackFunction(this, &DeformableModelApplication::ProcessSurfaceViewInteraction);
  //  m_SimplexMeshViewer.AddObserver(ClickedPointEvent(), m_SurfaceViewerCommand);


  const float alpha = 0.8; //internal forces
  const float beta  = 0.8; //external forces
  const float gamma = 0.35; 
  const int   range = 2; // how far to go through scan line algorithm
  const int   rigidity = 0; // regularization
  const int   iterations = 100;

//  sprintf(m_MessageString, "%4.1f", alpha);
//  m_InternalForceValueInput->value(m_MessageString);
//  sprintf(m_MessageString, "%4.1f", beta);
//  m_ExternalForceValueInput->value(m_MessageString);
//  sprintf(m_MessageString, "%5.2f", gamma);
//  m_GammaForceValueInput->value(m_MessageString);
//  sprintf(m_MessageString, "%d", range);
//  m_RangeForceValueInput->value(m_MessageString);
//  sprintf(m_MessageString, "%d", rigidity);
//  m_RigidityForceValueInput->value(m_MessageString);
//  sprintf(m_MessageString, "%d", iterations);
//  m_IterationsValueInput->value(m_MessageString);
 
  m_ImageLoaded = false;
  m_MeshLoaded = false;
  m_MeshCreated = false;
  m_PreprocessingFinished = false;
  m_PreprocessedLoaded = false;
  m_WhichImage = 0;

}

DeformableModelApplication
::~DeformableModelApplication()
{
  if( m_ShiftScaleImageFilter )
    {
      m_ShiftScaleImageFilter->Delete();
    }
}

void
DeformableModelApplication
::Show()
{
 
  mainWindow->show();
  axialView->show();
  coronalView->show();
  sagittalView->show();
  surfaceView->show();
  
  m_AxialViewer.SetInteractor( axialView );
  m_CoronalViewer.SetInteractor( coronalView );
  m_SagittalViewer.SetInteractor( sagittalView );
  m_SimplexMeshViewer.SetInteractor( surfaceView );
 
  axialView->Initialize();
  coronalView->Initialize();
  sagittalView->Initialize();
  surfaceView->Initialize();

}

void 
DeformableModelApplication
::Hide()
{
  mainWindow->hide();
}

void 
DeformableModelApplication
::Quit()
{
  this->Hide();
}

void 
DeformableModelApplication
::LoadMesh()
{
  if (m_ImageLoaded )
    {
      if (m_MeshCreated) 
        {
          return;
        }
    }
  else
    {
      return;
    }
  const char * filename = fl_file_chooser("Mesh filename","*.pri","");
 
  if( !filename  || strlen(filename) == 0 )
    {
    return;
    }
  
///////////////added by me//////////////////////////////////////
	Sur->ReadPrior(filename);
	//Sur->DWT(6,0);
	//Sur->MeshOut("test.wrl",4);
	//Sur->GetMean();
	/*clock_t start,finish;
	double totaltime;
	start=clock();*/
	/*finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"\n此程序的运行时间为"<<totaltime<<"秒！"<<endl;*/
	Sur->SetImage(m_LoadedPreprocessed);
	Sur->GetMean();
	Sur->SetBackup();
	Sur->Inv_DWT(0,5);
	//Sur->DWT(5,0);
	//Sur->SaveCof("test.cof");
	
	
	//Sur->MeshToITK(m_TriangleMeshSource,6);
////////////////////////////////////////////////////////////////

  //m_TriangleMesh = m_TriangleMeshSource->GetOutput();
	m_TriangleMesh = Sur->GenerateTriangleMesh(5);
 
  std::cout << "Number of Points =   " << m_TriangleMesh->GetNumberOfPoints() << std::endl;
  std::cout << "Number of Cells  =   " << m_TriangleMesh->GetNumberOfCells()  << std::endl;

  m_MeshLoaded = true;
  this->RefreshMeshVisualization();
   
  // force a redraw
  axialView->redraw();
  coronalView->redraw();
  sagittalView->redraw();
  surfaceView->redraw();

  Fl::check(); 

  cout<<"The initial fitness measurement is: "<<Sur->ComputeFitness()<<endl;

  //delete all VTK objects
//  vmarchingcubes->Delete();
//  vdecimate->Delete();
}

void 
DeformableModelApplication
::LoadMesh(char* filename)
{
	if (m_ImageLoaded )
	{
		if (m_MeshCreated) 
		{
			return;
		}
	}
	else
	{
		return;
	}
	

	///////////////added by me//////////////////////////////////////
	Sur->ReadPrior(filename);
	//Sur->DWT(6,0);
	//Sur->MeshOut("test.wrl",4);
	//Sur->GetMean();
	/*clock_t start,finish;
	double totaltime;
	start=clock();*/
	/*finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"\n此程序的运行时间为"<<totaltime<<"秒！"<<endl;*/
	Sur->SetImage(m_LoadedPreprocessed);
	Sur->GetMean();
	Sur->SetBackup();
	Sur->Inv_DWT(0,5);
	//Sur->DWT(5,0);
	//Sur->SaveCof("test.cof");


	//Sur->MeshToITK(m_TriangleMeshSource,6);
	////////////////////////////////////////////////////////////////

	//m_TriangleMesh = m_TriangleMeshSource->GetOutput();
	m_TriangleMesh = Sur->GenerateTriangleMesh(5);

	std::cout << "Number of Points =   " << m_TriangleMesh->GetNumberOfPoints() << std::endl;
	std::cout << "Number of Cells  =   " << m_TriangleMesh->GetNumberOfCells()  << std::endl;

	m_MeshLoaded = true;
	this->RefreshMeshVisualization();

	// force a redraw
	axialView->redraw();
	coronalView->redraw();
	sagittalView->redraw();
	surfaceView->redraw();

	Fl::check(); 

	cout<<"The initial fitness measurement is: "<<Sur->ComputeFitness()<<endl;

	//delete all VTK objects
	//  vmarchingcubes->Delete();
	//  vdecimate->Delete();
}

void 
DeformableModelApplication
::LoadCoefficients() 
{
	if (m_ImageLoaded )
	{
		if (m_MeshCreated) 
		{
			return;
		}
	}
	else
	{
		return;
	}
	const char * filename = fl_file_chooser("COEF filename","*.coef","");

	if( !filename  || strlen(filename) == 0 )
	{
		return;
	}

	///////////////added by me//////////////////////////////////////
	Sur->ReadCof(filename);
	Sur->Inv_DWT(0,5);
	m_TriangleMesh = Sur->GenerateTriangleMesh(5);

	std::cout << "Number of Points =   " << m_TriangleMesh->GetNumberOfPoints() << std::endl;
	std::cout << "Number of Cells  =   " << m_TriangleMesh->GetNumberOfCells()  << std::endl;

	m_MeshLoaded = true;
	this->RefreshMeshVisualization();

	// force a redraw
	axialView->redraw();
	coronalView->redraw();
	sagittalView->redraw();
	surfaceView->redraw();

	Fl::check(); 
}

void
DeformableModelApplication
::CreateMesh()
{
  
   if (m_ImageLoaded )
    {
      if (m_MeshLoaded) 
        {
          return;
        }
    }
  else
    {
      return;
    }
  m_SphereMeshSource->SetCenter(m_SeedPoint);
  m_SimplexFilter->Update();
  
  m_SimplexMesh = m_SimplexFilter->GetOutput();
  m_SimplexMesh->DisconnectPipeline();

  m_SimplexMeshToShow = m_SimplexMesh;
  m_MeshCreated = true;
  this->RefreshMeshVisualization();
   
  // force a redraw
  axialView->redraw();
  coronalView->redraw();
  sagittalView->redraw();
  surfaceView->redraw();

  Fl::check(); 
}

void
DeformableModelApplication
::RefreshMeshVisualization()
{
   int numPoints =  m_TriangleMesh->GetNumberOfPoints();
  
 if (numPoints == 0)
   {
     fl_alert( "no points in Grid ");
     return; 
   }
 
 // fill in here the conversion between itkMesh versus vtkUnstructuredGrid

 vtkPolyData* vgrid = vtkPolyData::New();

 // Create the vtkPoints object and set the number of points
 vtkPoints* vpoints = vtkPoints::New();
 vpoints->SetNumberOfPoints(numPoints);
  
 // iterate over all the points in the itk mesh filling in
 // the vtkPoints object as we go
 TriangleMeshType::PointsContainer::Pointer points = m_TriangleMesh->GetPoints();
 for(TriangleMeshType::PointsContainer::Iterator i = points->Begin();
     i != points->End(); ++i)
   {
     // Get the point index from the point container iterator
     int idx = i->Index();
     // Set the vtk point at the index with the the coord array from itk
     // itk returns a const pointer, but vtk is not const correct, so
     // we have to use a const cast to get rid of the const
     vtkFloatingPointType * pp = const_cast<vtkFloatingPointType*>(i->Value().GetDataPointer());
      
     vpoints->SetPoint(idx, pp);
     
   }

 // Set the points on the vtk grid
 vgrid->SetPoints(vpoints);
  
 // it is probably better not to add scalar to the vis mesh
 // as it makes it harder to visualize
 //  vgrid->GetPointData()->SetScalars(scalars);
 //  vgrid->GetPointData()->CopyAllOn();
 
 TriangleMeshType::CellType::MultiVisitor::Pointer mv =
   TriangleMeshType::CellType::MultiVisitor::New();

  LineVisitor::Pointer lv = LineVisitor::New();
  PolygonVisitor::Pointer pv = PolygonVisitor::New();
  TriangleVisitor::Pointer tv = TriangleVisitor::New();
 
  //set up the visitors
  int vtkCellCount = 0; // running counter for current cell inserted into vtk
  int numCells = m_TriangleMesh->GetNumberOfCells();
  int *types = new int[numCells]; //type array for vtk
  bool onlyTriangles = false;
  //create vtk cells and estimate the size
   vtkCellArray* cells = vtkCellArray::New();
     cells->EstimateSize(numCells, 4);

    lv->SetTypeArray(types);
    lv->SetCellCounter(&vtkCellCount);
    lv->SetCellArray(cells);

    pv->SetTypeArray(types);
    pv->SetCellCounter(&vtkCellCount);
    pv->SetCellArray(cells);

    tv->SetTypeArray(types);
    tv->SetCellCounter(&vtkCellCount);
    tv->SetCellArray(cells);
  

  // Set the TypeArray CellCount and CellArray for both visitors
  //lv->SetTypeArray(types);
  // tv->SetTypeArray(types);
  //lv->SetCellCounter(&vtkCellCount);
    //  tv->SetCellCounter(&vtkCellCount);
  //lv->SetCellArray(cells);
    //  tv->SetCellArray(cells);
 
  // mv->AddVisitor(lv);
  //mv->AddVisitor(tv);
  m_TriangleMesh->Accept(mv);


  //vgrid->SetLines(cells);
  //vgrid->SetStrips(cells);  
   
if (onlyTriangles) {
      mv->AddVisitor(tv);
      m_TriangleMesh->Accept(mv);
      vgrid->SetStrips(cells);  
    }
    else 
    {
      mv->AddVisitor(tv);
      mv->AddVisitor(lv);
      //      mv->AddVisitor(pv);
  //    mv->AddVisitor(qv);
      // Now ask the mesh to accept the multivisitor which
      // will Call Visit for each cell in the mesh that matches the
      // cell types of the visitors added to the MultiVisitor
  m_TriangleMesh->Accept(mv);
      // Now set the cells on the vtk grid with the type array and cell array
  vgrid->SetPolys(cells);  
      //vgrid->SetStrips(cells);  
  vgrid->SetLines(cells);
    // Clean up vtk objects (no vtkSmartPointer ... )
  }

  cells->Delete();
  vpoints->Delete();
   
  //m_AxialViewer.SetSimplexMesh(vgrid);
  //m_CoronalViewer.SetSimplexMesh(vgrid);
  //m_SagittalViewer.SetSimplexMesh(vgrid);
  m_SimplexMeshViewer.SetSimplexMesh(vgrid);

  this->m_vtkMesh = vgrid;

 return;
  
}


void 
DeformableModelApplication
::ComputeInternalForces()
{
  if ( !m_ImageLoaded || (!m_MeshLoaded && !m_MeshCreated))
    {
      return;
    }
  std::cout << " Performing Preprocessing ... " << std::endl;
  m_CastImage->SetInput( m_VolumeReader->GetOutput() );
  m_CastImage->Update();
  /* removed this as it was too slow
  m_GradientAnisotropicImage->SetInput( m_CastImage->GetOutput());
  m_GradientAnisotropicImage->SetNumberOfIterations(5);
  m_GradientAnisotropicImage->SetTimeStep(0.0625);
  m_GradientAnisotropicImage->SetConductanceParameter(3);
  */
  //m_GradientMagnitude->SetInput( m_GradientAnisotropicImage->GetOutput() );
  m_GradientMagnitude->SetInput(m_CastImage->GetOutput() );
  m_GradientMagnitude->SetSigma(0.5);
  
  m_SigmoidImage->SetInput( m_GradientMagnitude->GetOutput());
  m_SigmoidImage->SetOutputMinimum(0);
  m_SigmoidImage->SetOutputMaximum(1);
  m_SigmoidImage->SetAlpha(230);
  m_SigmoidImage->SetBeta(1300);
  
  m_GradientFilter->SetInput( m_SigmoidImage->GetOutput());
  m_GradientFilter->SetSigma( 0.5);

  m_GradientFilter->Update();
  std::cout << " Preprocessing DONE!... " << std::endl;
  m_DeformFilter->SetGradient( m_GradientFilter->GetOutput() );
  VolumeType::Pointer   Image = m_VolumeReader->GetOutput();
  Image->DisconnectPipeline();
  m_DeformFilter->SetImage(Image );
  m_PreprocessingFinished = true;

}


void 
DeformableModelApplication::IterationCallback()
{
  // The visualization cannot be refreshed at every iteration because the
  // deformable model filter do not fully recomputes its output at every
  // iteration. It only does it when the iterations are completed.
  // 
  // this->DeformableModelApplicationBase::IterationCallback();
  //
}

void
DeformableModelApplication
::SaveMesh()
{
 
  const char * filename = fl_file_chooser("Save Mesh As OFF mesh","*.off","");
 
  if( !filename  || strlen(filename) == 0 )
    {
    return;
    }

  Sur->GetBackup();
  Sur->Inv_DWT(0,4);
  Sur->CC_Subd(4,5);
  Sur->MeshOut(filename, 5);

  filename = fl_file_chooser("Save Mesh As Wavelet Coefficients","*.coef","");

  if( !filename  || strlen(filename) == 0 )
  {
	  return;
  }

  Sur->GetBackup();
  Sur->SaveCof(filename);
  /*vtkPolyDataWriter *vpolywriter = vtkPolyDataWriter::New();
  vpolywriter->SetInput(m_SimplexMeshViewer.GetSimplexMesh());
  vpolywriter->SetFileName(filename);
  vpolywriter->Write();

  vpolywriter->Delete();*/
      
}

void
DeformableModelApplication
::SaveMesh(char* meshFile, char* coefFile)
{

	Sur->GetBackup();
	Sur->Inv_DWT(0,4);
	Sur->CC_Subd(4,5);
	Sur->MeshOut(meshFile, 5);

	
	Sur->GetBackup();
	Sur->SaveCof(coefFile);
}

void
DeformableModelApplication
::SaveMask()
{
  if (!m_ImageLoaded || (!m_MeshLoaded && !m_MeshCreated) || !m_PreprocessingFinished)
    {
      return;
    }

  const char * filename = fl_file_chooser("Save Mask As","*.mhd","");
 
  if( !filename  || strlen(filename) == 0 )
    {
    return;
    }
  double orgn[3];
  orgn[0] = m_LoadedVolume->GetOrigin()[0];
  orgn[1] = m_LoadedVolume->GetOrigin()[1];
  orgn[2] = m_LoadedVolume->GetOrigin()[2];

  std::cout << "Converting Simplex Mesh to Triangle Mesh . . ." << std::endl;
  
  m_SimplexToTriangle->SetInput(m_DeformFilter->GetOutput());
  m_SimplexToTriangle->Update();
  
  TriangleMeshType::Pointer triangleMesh = m_SimplexToTriangle->GetOutput();
  triangleMesh->DisconnectPipeline();
  m_SimplexToTriangle->Delete();
  m_TriangleToImage->SetInput(triangleMesh);
  MeshPixelType::SizeType size;

  size[0] = m_LoadedVolume->GetBufferedRegion().GetSize()[0];
  size[1] = m_LoadedVolume->GetBufferedRegion().GetSize()[1];
  size[2] = m_LoadedVolume->GetBufferedRegion().GetSize()[2];
  m_TriangleToImage->SetSize(size);
  
 
  m_TriangleToImage->SetOrigin(orgn);
  // spacing remains (1,1,1) until we make a change to deformable model class
  float spacing[3];
  spacing[0] = 1;
  spacing[1] = 1;
  spacing[2] = 1;

  m_TriangleToImage->SetSpacing(spacing);

  m_TriangleToImage->Update();

  m_ImageWriter->SetInput(m_TriangleToImage->GetOutput() );
  m_ImageWriter->SetFileName(filename);
  m_ImageWriter->UseInputMetaDataDictionaryOn();
  m_ImageWriter->Update();
}

void 
DeformableModelApplication
::DeformMesh()
{
	//---------------for output the interim result------------------------//
	/*char IntrimResultName[50]; 
	sprintf(IntrimResultName,"Interim%d.off",IntrimResultCounter);
	IntrimResultCounter++;
	Sur->MeshOut(IntrimResultName,5);*/
	//---------------for output the interim result------------------------//

	/*clock_t start,finish;
	double totaltime;
	start=clock();*/
	int i;

	/*cout<<Level0Iter<<endl;
	cout<<Level0StepLen<<endl;
	cout<<Level0StepNum<<endl;*/

	
	for (i=1;i<=Level0Iter;i++)
	{
		Sur->DeformByVertax(0,int(Level0StepLen),int(Level0StepNum));
	}
	for (i=1;i<=Level1Iter;i++)
	{
		Sur->DeformByVertax(1,int(Level1StepLen),int(Level1StepNum));
	}
	for (i=1;i<=Level2Iter;i++)
	{
		Sur->DeformByVertax(2,int(Level2StepLen),int(Level2StepNum));
	}
	for (i=1;i<=Level3Iter;i++)
	{
		Sur->DeformByVertax(3,int(Level3StepLen),int(Level3StepNum));
	}
	for (i=1;i<=Level4Iter;i++)
	{
		Sur->DeformByVertax(4,int(Level4StepLen),int(Level4StepNum));
	}

	//compute the optimized surface
	Sur->GetBackup();
	//Sur->SaveCof("best.coef");
	Sur->Inv_DWT(0,4);
	Sur->CC_Subd(4,5);

	m_TriangleMesh = Sur->GenerateTriangleMesh(5);
	this->RefreshMeshVisualization();
	// force a redraw
	axialView->redraw();
	coronalView->redraw();
	sagittalView->redraw();
	surfaceView->redraw();

	Fl::check();

	/*finish=clock();
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<"\n此程序的运行时间为"<<totaltime<<"秒！"<<endl;*/
	
	
	/*if (!m_ImageLoaded || (!m_MeshLoaded && !m_MeshCreated) || !m_PreprocessingFinished)
    {
      return;
    }
  SimplexMeshType::Pointer simplexMesh2 = m_SimplexMesh;
  const unsigned int numberOfIterationsToGo = atoi(m_IterationsValueInput->value());
 
  for( unsigned int i=0; i<numberOfIterationsToGo; i++ )
    {
      std::cout << " Iteration   " << i << std::endl;
      m_SimplexMesh->DisconnectPipeline();
      
      m_DeformFilter->SetInput( m_SimplexMesh );

      m_DeformFilter->SetIterations(1); 
      m_DeformFilter->Update();

      simplexMesh2 =  m_DeformFilter->GetOutput();
     
      m_SimplexMeshToShow  = m_SimplexMesh;

      this->RefreshMeshVisualization();
     
      // force a redraw
      axialView->redraw();
      coronalView->redraw();
      sagittalView->redraw();
      surfaceView->redraw();

      Fl::check(); 
    }
*/
  
}
   

void 
DeformableModelApplication
::Load()
{

  const char * filename = fl_file_chooser("Volume filename","*.mha","");
 
  if( !filename  || strlen(filename) == 0 )
    {
      return;
    }
  
  m_VolumeReader->SetFileName( filename );
  m_VolumeReader->Update();

  m_LoadedOriginal = m_VolumeReader->GetOutput();
  m_LoadedOriginal->DisconnectPipeline();
  m_LoadedVolume = m_LoadedOriginal;

  m_WhichImage = 1;
  this->LoadPostProcessing();

  m_ImageLoaded =  true;

 }

void 
DeformableModelApplication
::Load(char* imaFile)
{

	m_VolumeReader->SetFileName( imaFile );
	m_VolumeReader->Update();

	m_LoadedOriginal = m_VolumeReader->GetOutput();
	m_LoadedOriginal->DisconnectPipeline();
	m_LoadedVolume = m_LoadedOriginal;

	m_WhichImage = 1;
	this->LoadPostProcessing();

	m_ImageLoaded =  true;

}

//added by me start----------------------------------------
void 
DeformableModelApplication
::LoadPreprocessed()
{
	const char * filename = fl_file_chooser("Preprocessed filename","*.hdr","");

	if( !filename  || strlen(filename) == 0 )
	{
		return;
	}

	m_PreprocessedReader->SetFileName( filename );
	m_PreprocessedReader->Update();

	m_LoadedPreprocessed = m_PreprocessedReader->GetOutput();
	m_LoadedPreprocessed->DisconnectPipeline();
	m_LoadedVolume = m_LoadedPreprocessed;

	Sur->SetImage(m_LoadedPreprocessed);

	m_WhichImage = 2;
	this->LoadPostProcessing();
	
	//------------------------------------------------------------------------------------
	m_PreprocessedLoaded =  true;
}

void 
DeformableModelApplication
::LoadPreprocessed(char* imgFile)
{
	
	m_PreprocessedReader->SetFileName( imgFile );
	m_PreprocessedReader->Update();

	m_LoadedPreprocessed = m_PreprocessedReader->GetOutput();
	m_LoadedPreprocessed->DisconnectPipeline();
	m_LoadedVolume = m_LoadedPreprocessed;

	Sur->SetImage(m_LoadedPreprocessed);

	m_WhichImage = 2;
	this->LoadPostProcessing();

	//------------------------------------------------------------------------------------
	m_PreprocessedLoaded =  true;
}

void 
DeformableModelApplication
::ToggleImage()
{
	if (m_WhichImage==2)
	{
		//m_LoadedVolume = m_VolumeReader->GetOutput();
		m_LoadedVolume = m_LoadedOriginal;
		m_WhichImage = 1;
		this->LoadPostProcessing();		
	}
	else
	{
		//m_LoadedVolume = m_PreprocessedReader->GetOutput();
		m_LoadedVolume = m_LoadedPreprocessed;
		m_WhichImage = 2;
		this->LoadPostProcessing();		
	}
}
//added by me end------------------------------------------

void 
DeformableModelApplication
::LoadPostProcessing()
{

  /*if (m_WhichImage==1)
  {
	  m_RescaleIntensity->SetInput( m_LoadedOriginal );
  } 
  else if (m_WhichImage==2)
  {
	  m_RescaleIntensity->SetInput( m_LoadedPreprocessed );
  }
  else
  {
	  return;
  }*/
  
  m_RescaleIntensity->SetInput( m_LoadedVolume );
  m_RescaleIntensity->Update();

  try
    {
	  m_ITK2VTKAdaptor->SetInput( m_RescaleIntensity->GetOutput() );
	  m_ITK2VTKAdaptor->Update(); 
	  m_ShiftScaleImageFilter->SetInput(  m_ITK2VTKAdaptor->GetOutput() );
	  m_ShiftScaleImageFilter->Update();
      m_ShiftScaleImageFilter->UpdateWholeExtent();
    }
  catch( itk::ExceptionObject & exp )
    {
      fl_message( exp.GetDescription() );
      return;
    }

  
  // of type ImageSliceViewer
  m_AxialViewer.SetInput(    m_ShiftScaleImageFilter->GetOutput() );
  m_CoronalViewer.SetInput(  m_ShiftScaleImageFilter->GetOutput() );
  m_SagittalViewer.SetInput( m_ShiftScaleImageFilter->GetOutput() );
  m_SimplexMeshViewer.SetInput( m_ShiftScaleImageFilter->GetOutput() );

  const unsigned int numberOfZslices = m_LoadedVolume->GetBufferedRegion().GetSize()[2];
  const unsigned int numberOfYslices = m_LoadedVolume->GetBufferedRegion().GetSize()[1];
  const unsigned int numberOfXslices = m_LoadedVolume->GetBufferedRegion().GetSize()[0];

  axialViewSlider->bounds( 0.0, numberOfZslices-1 );
  axialViewSlider->value( numberOfZslices / 2 );
  this->SelectAxialSlice( numberOfZslices / 2 );

  coronalViewSlider->bounds( 0.0, numberOfYslices-1 );
  coronalViewSlider->value( numberOfYslices / 2 );
  this->SelectCoronalSlice( numberOfYslices / 2 );

  sagittalViewSlider->bounds( 0.0, numberOfXslices-1 );
  sagittalViewSlider->value( numberOfXslices / 2 );
  this->SelectSagittalSlice( numberOfXslices / 2 );

  m_AxialViewer.Render();
  m_CoronalViewer.Render();
  m_SagittalViewer.Render();
  m_SimplexMeshViewer.Render();
  
}



void 
DeformableModelApplication
::SelectAxialSlice( int slice )
{
  double *center;
  double *normal;
	
  m_AxialViewer.SelectSlice( slice );
  //////////////////////////added by me////////////////////////////////////
  m_SimplexMeshViewer.m_PlaneWidgetZ->SetSliceIndex(slice);
  m_SimplexMeshViewer.Render();

  center = m_SimplexMeshViewer.m_PlaneWidgetZ->GetCenter();
  normal = m_SimplexMeshViewer.m_PlaneWidgetZ->GetNormal();
  this->m_cutPlaneAxial->SetOrigin(slice+0.5,slice+0.5,slice+0.5);
  this->m_cutPlaneAxial->SetNormal(normal[0],normal[1],normal[2]);
  this->m_SurfaceCutter->SetInput(m_vtkMesh);
  this->m_SurfaceCutter->SetCutFunction(m_cutPlaneAxial);
  //m_SurfaceCutter->Update();
  this->m_vtkContour = m_SurfaceCutter->GetOutput();
  m_AxialViewer.SetSimplexMesh(m_vtkContour);
  //m_AxialViewer.Render();
  /////////////////////////////////////////////////////////////////////////
  axialView->redraw();
  Fl::check();
}



void 
DeformableModelApplication
::SelectCoronalSlice( int slice )
{
  double *center;
  double *normal;
	
  m_CoronalViewer.SelectSlice( slice );
  //////////////////////////added by me////////////////////////////////////
  m_SimplexMeshViewer.m_PlaneWidgetY->SetSliceIndex(slice);
  m_SimplexMeshViewer.Render();

  center = m_SimplexMeshViewer.m_PlaneWidgetY->GetCenter();
  normal = m_SimplexMeshViewer.m_PlaneWidgetY->GetNormal();
  this->m_cutPlaneCoronal->SetOrigin(slice+0.5,slice+0.5,slice+0.5);
  this->m_cutPlaneCoronal->SetNormal(normal[0],normal[1],normal[2]);
  this->m_SurfaceCutter->SetInput(m_vtkMesh);
  this->m_SurfaceCutter->SetCutFunction(m_cutPlaneCoronal);
  //m_SurfaceCutter->Update();
  this->m_vtkContour = m_SurfaceCutter->GetOutput();
  m_CoronalViewer.SetSimplexMesh(m_vtkContour);
  /////////////////////////////////////////////////////////////////////////
  coronalView->redraw();
  Fl::check();
}



void 
DeformableModelApplication
::SelectSagittalSlice( int slice )
{
  double *center;
  double *normal;
	
  m_SagittalViewer.SelectSlice( slice );
  //////////////////////////added by me////////////////////////////////////
  m_SimplexMeshViewer.m_PlaneWidgetX->SetSliceIndex(slice);
  m_SimplexMeshViewer.Render();

  center = m_SimplexMeshViewer.m_PlaneWidgetX->GetCenter();
  normal = m_SimplexMeshViewer.m_PlaneWidgetX->GetNormal();
  this->m_cutPlaneSagittal->SetOrigin(slice+0.5,slice+0.5,slice+0.5);
  this->m_cutPlaneSagittal->SetNormal(normal[0],normal[1],normal[2]);
  this->m_SurfaceCutter->SetInput(m_vtkMesh);
  this->m_SurfaceCutter->SetCutFunction(m_cutPlaneSagittal);
  //m_SurfaceCutter->Update();
  this->m_vtkContour = m_SurfaceCutter->GetOutput();
  m_SagittalViewer.SetSimplexMesh(m_vtkContour);
  /////////////////////////////////////////////////////////////////////////
  sagittalView->redraw();
  Fl::check();
}



void 
DeformableModelApplication
::ProcessAxialViewInteraction( void )
{
  m_AxialViewer.GetSelectPoint( m_SeedPoint );
  this->SyncAllViews();
}


void  
DeformableModelApplication
::ProcessCoronalViewInteraction( void )
{
  
  m_CoronalViewer.GetSelectPoint( m_SeedPoint );
  this->SyncAllViews();
}


void  
DeformableModelApplication
::ProcessSagittalViewInteraction( void )
{
  
  m_SagittalViewer.GetSelectPoint( m_SeedPoint );
  this->SyncAllViews();
}

void
DeformableModelApplication
::ProcessSurfaceViewInteraction( void )
{
  
  //m_SagittalViewer.GetSelectPoint( m_SeedPoint );
  //this->SyncAllViews();
}

void 
DeformableModelApplication
::SyncAllViews(void)
{
  itk::Point< double, 3 > point;
  point[0] = m_SeedPoint[0];
  point[1] = m_SeedPoint[1];
  point[2] = m_SeedPoint[2];

  
  VisualizationVolumeType::IndexType index; 

  m_LoadedVolume->TransformPhysicalPointToIndex( point, index );

  // If the point is outside the volume, do not change slice selection
  if( m_LoadedVolume->GetLargestPossibleRegion().IsInside( index ) )
    {
      axialViewSlider->value(    index[2]  );
      coronalViewSlider->value(  index[1]  );
      sagittalViewSlider->value( index[0]  );

      this->SelectAxialSlice(    index[2] );
      this->SelectCoronalSlice(  index[1] );
      this->SelectSagittalSlice( index[0] );
    }

  // Sync the selected point in all the views even if it is outside the image
  m_AxialViewer.SelectPoint(    m_SeedPoint[0], m_SeedPoint[1], m_SeedPoint[2] );
  m_CoronalViewer.SelectPoint(  m_SeedPoint[0], m_SeedPoint[1], m_SeedPoint[2] );
  m_SagittalViewer.SelectPoint( m_SeedPoint[0], m_SeedPoint[1], m_SeedPoint[2] );

}
