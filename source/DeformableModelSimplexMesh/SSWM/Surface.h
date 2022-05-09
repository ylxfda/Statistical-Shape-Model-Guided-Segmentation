#ifndef Surface_h  
#define Surface_h       1

#define PI 3.1415926535897932
#include "Mytree.h"
#include "QTreeR.h"
#include "Global.h"
#include "QTreeR.h"
#include "Quadtree.h"
#include "SetList.h"
#include <string>
#include <sstream>
#include "itkVertexCell.h"
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "itkAutomaticTopologyMeshSource.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "DeformableModelApplicationBase.h"
#include "DeformableModelApplication.h"
using namespace std;

extern float Alpha;

/******************************************/
/*    NODE LOCATION INSIDE THE QUADTREE   */
/******************************************/

const static char _INSIDE_NODE       = 0;
const static char _TOP_BORDER_NODE   = 1;
const static char _BOT_BORDER_NODE   = 2;
const static char _LEFT_BORDER_NODE  = 3;
const static char _RIGHT_BORDER_NODE = 4;
const static char _TR_CORNER_NODE    = 5;
const static char _TL_CORNER_NODE    = 6;
const static char _BR_CORNER_NODE    = 7;
const static char _BL_CORNER_NODE    = 8;
const static char _ROOT_NODE        = 9;

typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
typedef itk::Mesh<double,3, TriangleMeshTraits> TriangleMeshType;
typedef itk::AutomaticTopologyMeshSource< TriangleMeshType > TriangleMeshSourceType; 

class Surface
{
public:
	int NumBaseMesh;
	int NumLevel;
	int NumVertex;
	int NumSample; //the num of samples in the prior model.
	int NumLevPrior; //the num of levels in the prior model.
	DeformableModelApplicationBase::VolumeType* Image; //the image to fit with.
	QTreeR* ForestMesh;
	
	Surface(int NumBaseMesh, int NumLevel, int NumVertex, float Vertex[][3], short Faces[][4]);
	~Surface();

	void CC_Subd(int StartLevel, int EndLevel);
	void CC_Subd(int StartLevel);

	void DWT(int StartLevel);
	void DWT(int StartLevel, int EndLevel);
	void Inv_DWT(int StartLevel);
	void Inv_DWT(int StartLevel, int EndLevel);

	void SaveObj(const string& filename); //save the object to a file *.sur
	void ReadObj(const char* filename); //read the object from a file *.sur

	void SaveCof(const char* filename); //save the coefficients to a file *.cof
	void ReadCof(const char* filename); //read the coefficients from a file *.cof

	VertexList *Vertices;  //list of points in every level.
	//PointSet VertexSets; //set of all the points in the finest level
	WaveCoeffiList *Coefficients; //list of Wavelet Coefficients in every level. For Level0 is the scale function coeffi

	////////////////function to output the mesh file///////////////////////
	////// FileName: the name of the output mesh file
	////// Level: the level of the mesh to output (start from '0')
	///////////////////////////////////////////////////////////////////////
	void MeshOut(const string& filename, int Level);

	////////////////function to read in the mesh file in wrl format///////////////////////
	////// FileName: the name of the output mesh file
	////// Level: the level of the mesh to read (start from '0')
	///////////////////////////////////////////////////////////////////////
	void MeshIn(const char* filename, int Level);

	////////////////function to read in the prior model ///////////////////////
	////// FileName: the name of the output mesh file
	///////////////////////////////////////////////////////////////////////
	void ReadPrior(const char* filename);

	////////////////function to reconstruct the mean shape from prior model///////////////
	////// NumLev: the num of level in the prior model
	//////////////////////////////////////////////////////////////////////////////////////
	void GetMean();

	////////////////function to reconstruct the shape the specified sample///////////////
	////// NumIndex: the ID of the sample to reconstruct (start from '0')
	//////////////////////////////////////////////////////////////////////////////////////
	void GetSample(int NumIndex);

	////////////////function to output surface to a itk mesh///////////////////////
	////// ITKMeshSource: the name of the itk::mesh
	////// NumofLevel: the number of level of the mesh to output (start from '1')
	///////////////////////////////////////////////////////////////////////
	void MeshToITK(TriangleMeshSourceType::Pointer &ITKMeshSource, int NumofLevel);

	////////////generate the triangle mesh, in which only triangle cell inserted//////////
	//////////////////////////////////////////////////////////////////////////////////////
	TriangleMeshType::Pointer GenerateTriangleMesh(int Level);
	
	////////////////function to compute the coefficient at a vertex///////////////////////
	////// Point3D* &Vert: the pointer to the vertax
	////// Scores: the scores in pc.
	//////////////////////////////////////////////////////////////////////////////////////
	void ComputeCoefVertax(Point3D* &Vert, float Scores[3]);

	////////////////function to compute the coefficient at a vertex and save to backup//////////
	////// Point3D* &Vert: the pointer to the vertax
	////// Scores: the scores in pc.
	//////////////////////////////////////////////////////////////////////////////////////
	void ComputeCoefVertaxBackup(Point3D* &Vert, float Scores[3]);

	////////////////function to get back the coefficients from backup///////////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	void GetBackup(void);

	////////////////function to set the back from the current coefficients///////////////////
	////////////////////////////////////////////////////////////////////////////////////////
	void SetBackup(void);

	////////////////function to set the image to fit with//////////////////////////////////////////
	//ImageToFit: the image to fit with
	////////////////////////////////////////////////////////////////////////////////////////////////
	void SetImage(DeformableModelApplicationBase::VolumeType* ImageToFit);

	////////////////function to compute the fitness between image and the surface///////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	float ComputeFitness(void);

	////////////////function to compute the prior term in objective function////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	float PriTerm(void);

	////////////////function to deform the surface to fit image one vertax (coefficient) by one vertax///////
	////////////////////////////////////////////////////////////////////////////////////////////////
	void DeformByVertax(int Lev, int StepLen, int StepNum);

	////////////////function to set the current backup as the new initial guess/////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	void SetNewMean(void);

	////////////////function to compute the area of a quadrangle/////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////
	float QuadArea(Mynode* TheNode);

protected:
private:

	EdgeList *Edges;
	FaceList *Facets;	

	Point3D V_f(int treeIndx, int NodeIndx);
	Point3D V_fChildLevel(int treeIndx, int NodeIndx);

	/////////treeIndex is the # of tree in the forest, NodeIndx is the index of the father node
	Point3D E_f(int treeIndx, int NodeIndx);


	Point3D V_e(int treeIndx, int NodeIndx, char EdgeID);
	Point3D V_e_ChildLevel(int treeIndx, int NodeIndx, char EdgeID);

	///////TargetPoint pointing to the point which the result should be saved to////
	Point3D F_e(int treeIndx, int NodeIndx, char EdgeID);

	void FindEdgeNeighbor(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID, QTreeR* &NeighTreeR, int &NeighNodeIndex, char &NeighEdgeID);

	void EdgePoint(QTreeR* TreeR, int NodeIndex, char EdgeID, Point3D** &Start, Point3D** &End);

	Point3D ChildEdgeMidPoint(QTreeR* TreeR, int NodeIndex, char EdgeID);//return the value
	void ChildEdgeMidPoint(QTreeR* TreeR, int NodeIndex, char EdgeID, Point3D* &Pt);//return the pointer

	Point3D V_v(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID);
	//Point3D E_v(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID);
	Point3D F_v(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID);

	//////////////find 1-ring of a vertex and return the valence of the point///////
	int Find1Ring(QTreeR* TreeR, int NodeIndex, char EdgeID, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList);

	Point3D V_v(int Valence, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList);
	Point3D E_v(int Valence, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList);
	Point3D F_v(int Valence, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList);
	double TransX_mean, TransX_std; //the mean X-translation and std.
	double TransY_mean, TransY_std;
	double TransZ_mean, TransZ_std;
	double ScaleX_mean, ScaleX_std; //the mean X-rescale coeffi and std.
	double ScaleY_mean, ScaleY_std;
	double ScaleZ_mean, ScaleZ_std; 
	double RotationPhi_mean, RotationPhi_std; //the mean rotation angle and std.see http://mathworld.wolfram.com/EulerAngles.html
	double RotationTheta_mean, RotationTheta_std;
	double RotationTau_mean, RotationTau_std;
};
 
#endif