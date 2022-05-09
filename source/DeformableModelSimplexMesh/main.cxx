

#include "DeformableModelApplication.h"
#include "Surface.h"

float Vertex[8][3] = 
{
	{-6, -6, -6},
	{6, -6, -6},
	{6, -6, 6},
	{-6, -6, 6},
	{6, 6, -6},
	{6, 6, 6},
	{-6, 6, 6},
	{-6, 6, -6}
};

short Faces[6][4] = 
{
	{0,1,2,3},
	{2,5,6,3},
	{4,7,6,5},
	{1,0,7,4},
	{1,4,5,2},
	{0,3,6,7}
};

int IntrimResultCounter=0; //a counter for output the interim results
//float Alpha=145.0; // the weighting factor
float Alpha=45.0; // the weighting factor

Surface* Sur = new Surface(6, 6, 8, Vertex, Faces);

int main( int argc,  char *argv[] )
{

  if (argc==1)
  {
	  DeformableModelApplication app;

	  try
	  {
		  app.Show();
		  Fl::run();
	  }
	  catch( std::exception & ex )
	  {
		  std::cerr << ex.what() << std::endl;
	  }

	  return 0;
  } 
  else
  {
	  int count;

	  printf("The command line has %d arguments:\n", argc - 1);
	  for(count = 1; count < argc; count++)
	  {
			printf("%d: %s\n", count, argv[count]);
	  }
	  
	  DeformableModelApplication app;

      app.Show();
	  
	  app.Load(argv[2]);
	  app.LoadPreprocessed(argv[3]);
	  app.LoadMesh(argv[4]);
	  app.DeformMesh();
	  app.SaveMesh(argv[5],argv[6]);
	  //Fl::run();
	  return 0;	
  }  
}


