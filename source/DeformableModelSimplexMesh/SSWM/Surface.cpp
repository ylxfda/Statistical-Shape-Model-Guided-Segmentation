
#include "Global.h"
#include "Surface.h"
#include "Piont3D.h"
#include "SetList.h"
#include "Quadtree.h"
#include <cmath>
#include <string>        // for strings
#include <iostream>      // for I/O
#include <fstream>       // for file I/O
#include <iomanip>       // for setw()
#include <cstdlib>       // for exit()
#include <stdio.h>
#include <stdlib.h>
#include <time.h>



Surface::Surface(int NumBaseMesh, int NumLevel, int NumVertex, float Vertex[][3], short Faces[][4])
{
	std::cout<<"Initializing..."<<std::endl;
	
	int i,j,k;
	Mynode* HeadNode;

	this->NumSample = 0;
	this->NumLevPrior = 0;
	this->NumBaseMesh = NumBaseMesh;
	this->NumLevel = NumLevel;
	this->NumVertex = NumVertex;

	Point3D* BaseVertex = new Point3D [NumVertex];
	for (i=0; i<NumVertex; i++)
	{
		BaseVertex[i].SetLocation(Vertex[i][0],Vertex[i][1],Vertex[i][2]);
	}
	
	ForestMesh = new QTreeR[NumBaseMesh];
	
	for (i=0; i<NumBaseMesh; i++) 
	{
		ForestMesh[i].tree = new Mytree(NumLevel);
		HeadNode = ForestMesh[i].tree->getHeadNodes();
		HeadNode->TopLeft = &BaseVertex[Faces[i][0]];
		HeadNode->BotLeft = &BaseVertex[Faces[i][1]];
		HeadNode->BotRigh = &BaseVertex[Faces[i][2]];
		HeadNode->TopRigh = &BaseVertex[Faces[i][3]];
	}
	//initialize the affine transfore parameter..........................................
	TransX_mean=0.0, TransX_std=0.0; //the mean X-translation and std.
	TransY_mean=0.0, TransY_std=0.0;
	TransZ_mean=0.0, TransZ_std=0.0;
	ScaleX_mean=1.0, ScaleX_std=0.0; //the mean X-rescale coeffi and std.
	ScaleY_mean=1.0, ScaleY_std=0.0;
	ScaleZ_mean=1.0, ScaleZ_std=0.0; 
	RotationPhi_mean=0.0, RotationPhi_std=0.0; //the mean rotation angle and std.see http://mathworld.wolfram.com/EulerAngles.html
	RotationTheta_mean=0.0, RotationTheta_std=0.0;
	RotationTau_mean=0.0, RotationTau_std=0.0;
///////////////find the west neibor/////////////////////////////////////////////
	for (i=0; i<NumBaseMesh; i++) //the base-face to find its neib
	{
		int found = 0; //flag show wether found
		int faceN = 0; //which face 
		int edgeN = 0; //which edge
		for ( j=0; j<NumBaseMesh; j++)//check in the jth face
		{
			int startP = 0;       //start point of the edge being checked
			int endP = 0;         //end point of the edge being checked
			for ( k=0; k<=3; k++) 
			{
				startP = k;
				endP = (startP+1)%4;
				if ((Faces[i][0]==Faces[j][endP])&&(Faces[i][1]==Faces[j][startP]))
				{
					faceN = j;
					edgeN = k;
					found = 1;
					break;
				}
			}
			if (found == 1)
			{
				break;
			}			
		}
		if (found==1)
		{
			ForestMesh[i].NeibWest = &ForestMesh[j];			
			ForestMesh[i].ID_NeibWest = k+1;
		}
	}

///////////////find the south neibor/////////////////////////////////////////////
	for (i=0; i<NumBaseMesh; i++) //the base-face to find its neib
	{
		int found = 0; //flag show wether found
		int faceN = 0; //which face 
		int edgeN = 0; //which edge
		for ( j=0; j<NumBaseMesh; j++)//check in the jth face
		{
			int startP = 0;       //start point of the edge being checked
			int endP = 0;         //end point of the edge being checked
			for ( k=0; k<=3; k++) 
			{
				startP = k;
				endP = (startP+1)%4;
				if ((Faces[i][1]==Faces[j][endP])&&(Faces[i][2]==Faces[j][startP]))
				{
					faceN = j;
					edgeN = k;
					found = 1;
					break;
				}
			}
			if (found == 1)
			{
				break;
			}			
		}
		if (found==1)
		{
			ForestMesh[i].NeibSouth = &ForestMesh[j];			
			ForestMesh[i].ID_NeibSouth = k+1;
		}
	}

///////////////find the east neibor/////////////////////////////////////////////
	for (i=0; i<NumBaseMesh; i++) //the base-face to find its neib
	{
		int found = 0; //flag show wether found
		int faceN = 0; //which face 
		int edgeN = 0; //which edge
		for ( j=0; j<NumBaseMesh; j++)//check in the jth face
		{
			int startP = 0;       //start point of the edge being checked
			int endP = 0;         //end point of the edge being checked
			for ( k=0; k<=3; k++) 
			{
				startP = k;
				endP = (startP+1)%4;
				if ((Faces[i][2]==Faces[j][endP])&&(Faces[i][3]==Faces[j][startP]))
				{
					faceN = j;
					edgeN = k;
					found = 1;
					break;
				}
			}
			if (found == 1)
			{
				break;
			}			
		}
		if (found==1)
		{
			ForestMesh[i].NeibEast = &ForestMesh[j];			
			ForestMesh[i].ID_NeibEast = k+1;
		}
	}

///////////////find the north neibor/////////////////////////////////////////////
	for (i=0; i<NumBaseMesh; i++) //the base-face to find its neib
	{
		int found = 0; //flag show wether found
		int faceN = 0; //which face 
		int edgeN = 0; //which edge
		for ( j=0; j<NumBaseMesh; j++)//check in the jth face
		{
			int startP = 0;       //start point of the edge being checked
			int endP = 0;         //end point of the edge being checked
			for ( k=0; k<=3; k++) 
			{
				startP = k;
				endP = (startP+1)%4;
				if ((Faces[i][3]==Faces[j][endP])&&(Faces[i][0]==Faces[j][startP]))
				{
					faceN = j;
					edgeN = k;
					found = 1;
					break;
				}
			}
			if (found == 1)
			{
				break;
			}			
		}
		if (found==1)
		{
			ForestMesh[i].NeibNorth = &ForestMesh[j];			
			ForestMesh[i].ID_NeibNorth = k+1;
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////
	//  initial the point3d pointers at every node/////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	Mynode* Quatrangle;
	int NumNodes;
	int valence;

	QTreeR* ThisTreeR = NULL;
	int ThisNodeIndex; 
	char ThisEdgeID;

	QTreeR* NeighTreeR = NULL; 
//	char NeighEdgeID;

	list<QTreeR*> TreeList;
	list<int> NodeIndexList; 
	list<char> EdgesList;
	list<QTreeR*>::iterator posTree;
	list<int>::iterator posNode;
	list<char>::iterator posEdge;

	Point3D** Start;
	Point3D** End;
	Point3D* Temp;

   //int test=0; 

	if (NumLevel>1)
	{
		NumNodes = (int)((pow(4.0,NumLevel)-1)/3);
		for (j=0;j<NumBaseMesh;j++)//the # of the trees
		{
			for (k=1;k<NumNodes;k++)//index of the node
			{
				Quatrangle = ForestMesh[j].tree->getNode(k);				
				ThisTreeR = &(ForestMesh[j]);
				ThisNodeIndex = k; 
				
				TreeList.clear();
				NodeIndexList.clear();
				EdgesList.clear();
				if (Quatrangle->TopLeft==NULL)
				{
					ThisEdgeID = NORTH;
					valence = Find1Ring(ThisTreeR, ThisNodeIndex, ThisEdgeID, TreeList, NodeIndexList, EdgesList);    
					Temp = new Point3D;

					//test++;

					posTree = TreeList.begin();
					posNode = NodeIndexList.begin();
					posEdge = EdgesList.begin();
					for (i=0;i<valence;i++)
					{
						EdgePoint(*posTree, *posNode, *posEdge, Start, End);
						(*End) = Temp;
						posTree++;
						posNode++;
						posEdge++;
					}
				}
				
				TreeList.clear();
				NodeIndexList.clear();
				EdgesList.clear();
				if (Quatrangle->TopRigh==NULL)
				{
					ThisEdgeID = EAST;
					valence = Find1Ring(ThisTreeR, ThisNodeIndex, ThisEdgeID, TreeList, NodeIndexList, EdgesList);    
					Temp = new Point3D;

					//test++;

					posTree = TreeList.begin();
					posNode = NodeIndexList.begin();
					posEdge = EdgesList.begin();
					for (i=0;i<valence;i++)
					{
						EdgePoint(*posTree, *posNode, *posEdge, Start, End);
						(*End) = Temp;
						posTree++;
						posNode++;
						posEdge++;
					}
				}

				TreeList.clear();
				NodeIndexList.clear();
				EdgesList.clear();
				if (Quatrangle->BotLeft==NULL)
				{
					ThisEdgeID = WEST;
					valence = Find1Ring(ThisTreeR, ThisNodeIndex, ThisEdgeID, TreeList, NodeIndexList, EdgesList);    
					Temp = new Point3D;

					//test++;

					posTree = TreeList.begin();
					posNode = NodeIndexList.begin();
					posEdge = EdgesList.begin();
					for (i=0;i<valence;i++)
					{
						EdgePoint(*posTree, *posNode, *posEdge, Start, End);
						(*End) = Temp;
						posTree++;
						posNode++;
						posEdge++;
					}
				}

				TreeList.clear();
				NodeIndexList.clear();
				EdgesList.clear();
				if (Quatrangle->BotRigh==NULL)
				{
					ThisEdgeID = SOUTH;
					valence = Find1Ring(ThisTreeR, ThisNodeIndex, ThisEdgeID, TreeList, NodeIndexList, EdgesList);    
					Temp = new Point3D;

					//test++;

					posTree = TreeList.begin();
					posNode = NodeIndexList.begin();
					posEdge = EdgesList.begin();
					for (i=0;i<valence;i++)
					{
						EdgePoint(*posTree, *posNode, *posEdge, Start, End);
						(*End) = Temp;
						posTree++;
						posNode++;
						posEdge++;
					}
				}
			}				
		}
	}

	//cout<<"test="<<test<<endl;

	/////////////////////////////////Create the Vertex Edge Face list for every level/////////////////////////
	Vertices = new VertexList [NumLevel];
	Edges = new EdgeList [NumLevel];
 	Facets = new FaceList [NumLevel];
	Coefficients = new WaveCoeffiList [NumLevel];
	
	PointSet SetVertex;// a set used for checking whether the element has already been added
	set<PointSet> SetEdge;
	PointSet ThisEdge;
	set<Mynode*> SetFaces;	
	int EndIndex;
	int StartIndex;

	Vertax Ver;
	Edge Edg;
	Facet Fac;
	WaveCoeffi Coe;

	Point3D* CoeffFound;

	for (i=0;i<NumLevel;i++)
	{
		SetVertex.clear();
		SetEdge.clear();
		SetFaces.clear();
		EndIndex = (int)((pow(4.0, i+1)-1)/3-1);
		StartIndex = (int)(EndIndex - pow(pow(2.0, i), 2) + 1);
		for (j=0;j<NumBaseMesh;j++)
		{
			for (k=StartIndex;k<=EndIndex;k++)
			{
				Quatrangle = ForestMesh[j].tree->getNode(k);
				Ver.TreeIndex = j;
				Ver.NodeIndex = k;
				Edg.TreeIndex = j;
				Edg.NodeIndex = k;
				Fac.TreeIndex = j;
				Fac.NodeIndex = k;
				Coe.TreeIndex = j;
				Coe.NodeIndex = k;
				////////////////insert vertex///////////////////////////
				if (SetVertex.insert(Quatrangle->TopLeft).second)
				{					
					Ver.Position = TL;
					Vertices[i].push_back(Ver);
					if (i==0)
					{
						Coe.TypeID = ScaleCoffi;
						Coe.Coeffi = Quatrangle->TopLeft;
						Coefficients[0].push_back(Coe);
					}
				}
				if (SetVertex.insert(Quatrangle->TopRigh).second)
				{					
					Ver.Position = TR;
					Vertices[i].push_back(Ver);
					if (i==0)
					{
						Coe.TypeID = ScaleCoffi;
						Coe.Coeffi = Quatrangle->TopRigh;
						Coefficients[0].push_back(Coe);
					}
				}
				if (SetVertex.insert(Quatrangle->BotLeft).second)
				{					
					Ver.Position = BL;
					Vertices[i].push_back(Ver);
					if (i==0)
					{
						Coe.TypeID = ScaleCoffi;
						Coe.Coeffi = Quatrangle->BotLeft;
						Coefficients[0].push_back(Coe);
					}
				}
				if (SetVertex.insert(Quatrangle->BotRigh).second)
				{					
					Ver.Position = BR;
					Vertices[i].push_back(Ver);
					if (i==0)
					{
						Coe.TypeID = ScaleCoffi;
						Coe.Coeffi = Quatrangle->BotRigh;
						Coefficients[0].push_back(Coe);
					}
				}

				///////// insert the edges ///////////////////////////
				ThisEdge.insert(Quatrangle->TopLeft);
				ThisEdge.insert(Quatrangle->TopRigh);
				if (SetEdge.insert(ThisEdge).second)
				{					
					Edg.EdgeID = NORTH;
					Edges[i].push_back(Edg);
					if (i<NumLevel-1)
					{
						Coe.TypeID = NorthEdgeCoeffi;
						ChildEdgeMidPoint(&ForestMesh[j], k, NORTH, CoeffFound);
						Coe.Coeffi = CoeffFound;
						Coefficients[i+1].push_back(Coe);
					}
				}
				ThisEdge.clear();

				ThisEdge.insert(Quatrangle->BotLeft);
				ThisEdge.insert(Quatrangle->BotRigh);
				if (SetEdge.insert(ThisEdge).second)
				{					
					Edg.EdgeID = SOUTH;
					Edges[i].push_back(Edg);
					if (i<NumLevel-1)
					{
						Coe.TypeID = SouthEdgeCoeffi;
						ChildEdgeMidPoint(&ForestMesh[j], k, SOUTH, CoeffFound);
						Coe.Coeffi = CoeffFound;
						Coefficients[i+1].push_back(Coe);
					}
				}
				ThisEdge.clear();

				ThisEdge.insert(Quatrangle->TopLeft);
				ThisEdge.insert(Quatrangle->BotLeft);
				if (SetEdge.insert(ThisEdge).second)
				{					
					Edg.EdgeID = WEST;
					Edges[i].push_back(Edg);
					if (i<NumLevel-1)
					{
						Coe.TypeID = WestEdgeCoeffi;
						ChildEdgeMidPoint(&ForestMesh[j], k, WEST, CoeffFound);
						Coe.Coeffi = CoeffFound;
						Coefficients[i+1].push_back(Coe);
					}
				}
				ThisEdge.clear();

				ThisEdge.insert(Quatrangle->BotRigh);
				ThisEdge.insert(Quatrangle->TopRigh);
				if (SetEdge.insert(ThisEdge).second)
				{					
					Edg.EdgeID = EAST;
					Edges[i].push_back(Edg);
					if (i<NumLevel-1)
					{
						Coe.TypeID = EastEdgeCoeffi;
						ChildEdgeMidPoint(&ForestMesh[j], k, EAST, CoeffFound);
						Coe.Coeffi = CoeffFound;
						Coefficients[i+1].push_back(Coe);
					}
				}
				ThisEdge.clear();
				////////////////////////insert the faces////////////////////////////
				if (SetFaces.insert(Quatrangle).second)
				{					
					Facets[i].push_back(Fac);
					if (i<NumLevel-1)
					{
						Coe.TypeID = FacetCoeffi;
						Quatrangle = ForestMesh[j].tree->getNode(ForestMesh[j].tree->getTLChildNode(k));
						Coe.Coeffi = Quatrangle->BotRigh;
						Coefficients[i+1].push_back(Coe);
					}
				}
			}
		}
	}

	/*cout<<Vertices[0].size()<<endl;
	cout<<Edges[0].size()<<endl;
	cout<<Facets[0].size()<<endl;*/
	//cout<<"No of coeffi in level 0: "<<Coefficients[0].size()<<endl;
	cout<<"Initialization successful!"<<endl;
	cout<<"There are "<<this->NumLevel<<" levels."<<endl;
	cout<<"There are "<<this->Vertices[this->NumLevel-1].size()<<" vertex in the finest level."<<endl;

	//creat the mesh points sets in the finest level...........................................
	/*
	EndIndex = (int)((pow(4.0, NumLevel-1+1)-1)/3-1);
	StartIndex = (int)(EndIndex - pow(pow(2.0, NumLevel-1), 2) + 1);

	for (i=0; i<NumBaseMesh; i++)
	{
		for (j=StartIndex; j<=EndIndex; j++) 
		{
			VertexSets.insert(this->ForestMesh[i].tree->getNode(j)->TopLeft); 
			VertexSets.insert(this->ForestMesh[i].tree->getNode(j)->TopRigh);
			VertexSets.insert(this->ForestMesh[i].tree->getNode(j)->BotLeft);
			VertexSets.insert(this->ForestMesh[i].tree->getNode(j)->BotRigh);
		}		
	}	
	cout<<VertexSets.size()<<endl;*/
}

Surface::~Surface()
{
	for (int i=0; i<NumBaseMesh; i++) 
	{
		delete ForestMesh[i].tree;
	}
	delete[] ForestMesh;
}

void Surface::CC_Subd(int StartLevel, int EndLevel)
{
	for (int i=StartLevel;i<EndLevel;i++)
	{
		CC_Subd(i);
	}
}

void Surface::CC_Subd(int StartLevel)
{
	Mynode *ChildNode;
	Mynode *ThisNode;
	int ThisNodeIndex;
	int ChildNodeIndex;
	int ThisTreeIndex;
	int Valence;
	Point3D* ChildVertex;
	char EdgeID;
	char PointPosition;

	list<QTreeR*> TreeLt;
	list<int> NodeIndexLt; 
	list<char> EdgesLt;

	VertexList::iterator posVertex;
	EdgeList::iterator posEdge;
	FaceList::iterator posFace;

///////////////////Compute the f point in next level///////////////////////
	for (posFace=Facets[StartLevel].begin();posFace!=Facets[StartLevel].end();posFace++)
	{
		ThisTreeIndex = (*posFace).TreeIndex;
		ThisNodeIndex = (*posFace).NodeIndex;
		ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
		ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
		ChildVertex = ChildNode->BotRigh;
		*ChildVertex = V_f(ThisTreeIndex, ThisNodeIndex);		
	}
////////////////////////////////////////////////////////////////////////////

/////////////////compute the e point in next level//////////////////////////
	for (posEdge=Edges[StartLevel].begin();posEdge!=Edges[StartLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = (V_e(ThisTreeIndex, ThisNodeIndex, EdgeID) + F_e(ThisTreeIndex, ThisNodeIndex, EdgeID))/2;		
	}
///////////////////////////////////////////////////////////////////////////

/////////////////compute the v point in next level/////////////////////////
	for (posVertex=Vertices[StartLevel].begin();posVertex!=Vertices[StartLevel].end();posVertex++)
	{
		ThisTreeIndex = (*posVertex).TreeIndex;
		ThisNodeIndex = (*posVertex).NodeIndex;
		PointPosition = (*posVertex).Position;
		ThisNode = ForestMesh[ThisTreeIndex].tree->getNode(ThisNodeIndex);
		TreeLt.clear();
		NodeIndexLt.clear();
		EdgesLt.clear();
		switch(PointPosition)
		{
		case TL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, NORTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft =(*(*ThisNode).TopLeft)*(Valence-2) + F_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft = *(*ChildNode).TopLeft + V_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft = (*(*ChildNode).TopLeft)/Valence;
			break;
		case TR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EAST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh =(*(*ThisNode).TopRigh)*(Valence-2) + F_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh = *(*ChildNode).TopRigh + V_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh = (*(*ChildNode).TopRigh)/Valence;
			break;
		case BL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, WEST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft =(*(*ThisNode).BotLeft)*(Valence-2) + F_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft = *(*ChildNode).BotLeft + V_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft = (*(*ChildNode).BotLeft)/Valence;
			break;
		case BR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, SOUTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh =(*(*ThisNode).BotRigh)*(Valence-2) + F_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh = *(*ChildNode).BotRigh + V_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh = (*(*ChildNode).BotRigh)/Valence;
			break;
		}		
	}
}

void Surface::DWT(int StartLevel)
{
	int TargetLevel=StartLevel-1;
	int ThisTreeIndex;
	int ThisNodeIndex;
	int Valence;
	Mynode *ThisNode = NULL;
	int ChildNodeIndex;
	Mynode *ChildNode = NULL;
	char PointPosition;
	char EdgeID;
	Point3D* ChildVertex;

	list<QTreeR*> TreeLt;
	list<int> NodeIndexLt;
	list<char> EdgesLt;

	VertexList::iterator posVertex;
	EdgeList::iterator posEdge;
	FaceList::iterator posFace;
	/*----------------------------- Lift 1------------------------------------------*/
	for (posVertex=Vertices[TargetLevel].begin();posVertex!=Vertices[TargetLevel].end();posVertex++)
	{
		ThisTreeIndex = (*posVertex).TreeIndex;
		ThisNodeIndex = (*posVertex).NodeIndex;
		PointPosition = (*posVertex).Position;
		TreeLt.clear();
		NodeIndexLt.clear();
		EdgesLt.clear();
		switch(PointPosition)
		{
		case TL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, NORTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft =*(*ChildNode).TopLeft +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).TopLeft = *(*ChildNode).TopLeft - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		case TR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EAST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh =*(*ChildNode).TopRigh +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).TopRigh = *(*ChildNode).TopRigh - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		case BL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, WEST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft =*(*ChildNode).BotLeft +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).BotLeft = *(*ChildNode).BotLeft - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		case BR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, SOUTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh =*(*ChildNode).BotRigh +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).BotRigh = *(*ChildNode).BotRigh - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		}		
	}

	/*------------------------------Lift 2------------------------------------------*/
	for (posEdge=Edges[TargetLevel].begin();posEdge!=Edges[TargetLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = *ChildVertex - (F_e(ThisTreeIndex, ThisNodeIndex, EdgeID)/2);
	}
	/*----------------------------- Lift 3------------------------------------------*/
	for (posFace=Facets[TargetLevel].begin();posFace!=Facets[TargetLevel].end();posFace++)
	{
		ThisTreeIndex = (*posFace).TreeIndex;
		ThisNodeIndex = (*posFace).NodeIndex;
		ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
		ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
		ChildVertex = ChildNode->BotRigh;
		*ChildVertex = *ChildVertex + V_fChildLevel(ThisTreeIndex, ThisNodeIndex)*4 - E_f(ThisTreeIndex, ThisNodeIndex)*4;		
	}
	/*------------------------------Lift 4------------------------------------------*/
	for (posEdge=Edges[TargetLevel].begin();posEdge!=Edges[TargetLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = *ChildVertex - (V_e_ChildLevel(ThisTreeIndex, ThisNodeIndex, EdgeID)*2);		
	}
	/*----------------------------- Lift 5------------------------------------------*/
	for (posVertex=Vertices[TargetLevel].begin();posVertex!=Vertices[TargetLevel].end();posVertex++)
	{
		ThisTreeIndex = (*posVertex).TreeIndex;
		ThisNodeIndex = (*posVertex).NodeIndex;
		PointPosition = (*posVertex).Position;
		TreeLt.clear();
		NodeIndexLt.clear();
		EdgesLt.clear();
		ThisNode = ForestMesh[ThisTreeIndex].tree->getNode(ThisNodeIndex);
		switch(PointPosition)
		{
		case TL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, NORTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft =(*(*ChildNode).TopLeft)*4 +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/16;
			*(*ChildNode).TopLeft = *(*ChildNode).TopLeft + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3;
			*(ThisNode->TopLeft) = *(*ChildNode).TopLeft;
			break;
		case TR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EAST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh =(*(*ChildNode).TopRigh)*4 +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/16;
			*(*ChildNode).TopRigh = *(*ChildNode).TopRigh + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3;
			*(ThisNode->TopRigh) = *(*ChildNode).TopRigh;
			break;
		case BL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, WEST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft =(*(*ChildNode).BotLeft)*4 +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/16;
			*(*ChildNode).BotLeft = *(*ChildNode).BotLeft + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3;
			*(ThisNode->BotLeft) = *(*ChildNode).BotLeft;
			break;
		case BR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, SOUTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh =(*(*ChildNode).BotRigh)*4 +	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/16;
			*(*ChildNode).BotRigh = *(*ChildNode).BotRigh + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3;
			*(ThisNode->BotRigh) = *(*ChildNode).BotRigh;
			break;
		}		
	}
	/*------------------------------Lift 6------------------------------------------*/
	for (posEdge=Edges[TargetLevel].begin();posEdge!=Edges[TargetLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = (*ChildVertex)*2 + (F_e(ThisTreeIndex, ThisNodeIndex, EdgeID)*3/4);		
	}
}

void Surface::DWT(int StartLevel, int EndLevel)
{
	for (int i=StartLevel;i>EndLevel;i--)
	{
		DWT(i);
	}
}

void Surface::Inv_DWT(int StartLevel)
{
	int TargetLevel=StartLevel+1;
	int ThisTreeIndex;
	int ThisNodeIndex;
	int Valence;
	Mynode *ThisNode = NULL;
	int ChildNodeIndex;
	Mynode *ChildNode = NULL;
	char PointPosition;
	char EdgeID;
	Point3D* ChildVertex;

	list<QTreeR*> TreeLt;
	list<int> NodeIndexLt;
	list<char> EdgesLt;

	VertexList::iterator posVertex;
	EdgeList::iterator posEdge;
	FaceList::iterator posFace;
	
	/*----------------------------- Lift 1------------------------------------------*/
	for (posEdge=Edges[StartLevel].begin();posEdge!=Edges[StartLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = (*ChildVertex)/2 - (F_e(ThisTreeIndex, ThisNodeIndex, EdgeID)*3/8);
	}
	/*------------------------------Lift 2------------------------------------------*/
	for (posVertex=Vertices[StartLevel].begin();posVertex!=Vertices[StartLevel].end();posVertex++)
	{
		ThisTreeIndex = (*posVertex).TreeIndex;
		ThisNodeIndex = (*posVertex).NodeIndex;
		PointPosition = (*posVertex).Position;
		TreeLt.clear();
		NodeIndexLt.clear();
		EdgesLt.clear();
		ThisNode = ForestMesh[ThisTreeIndex].tree->getNode(ThisNodeIndex);
		switch(PointPosition)
		{
		case TL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, NORTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft = *(*ThisNode).TopLeft/4 - F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/64;
			*(*ChildNode).TopLeft = *(*ChildNode).TopLeft - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3/4;
			break;
		case TR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EAST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh =*(*ThisNode).TopRigh/4 - F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/64;
			*(*ChildNode).TopRigh = *(*ChildNode).TopRigh - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3/4;
			break;
		case BL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, WEST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft =*(*ThisNode).BotLeft/4 - F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/64;
			*(*ChildNode).BotLeft = *(*ChildNode).BotLeft - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3/4;
			break;
		case BR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, SOUTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh =*(*ThisNode).BotRigh/4 - F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*9/64;
			*(*ChildNode).BotRigh = *(*ChildNode).BotRigh - E_v(Valence, TreeLt, NodeIndexLt, EdgesLt)*3/4;
			break;
		}		
	}
	/*----------------------------- Lift 3------------------------------------------*/
	for (posEdge=Edges[StartLevel].begin();posEdge!=Edges[StartLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = *ChildVertex + (V_e_ChildLevel(ThisTreeIndex, ThisNodeIndex, EdgeID)*2);		
	}
	/*------------------------------Lift 4------------------------------------------*/
	for (posFace=Facets[StartLevel].begin();posFace!=Facets[StartLevel].end();posFace++)
	{
		ThisTreeIndex = (*posFace).TreeIndex;
		ThisNodeIndex = (*posFace).NodeIndex;
		ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
		ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
		ChildVertex = ChildNode->BotRigh;
		*ChildVertex = *ChildVertex - V_fChildLevel(ThisTreeIndex, ThisNodeIndex)*4 + E_f(ThisTreeIndex, ThisNodeIndex)*4;		
	}
	/*----------------------------- Lift 5------------------------------------------*/
	for (posEdge=Edges[StartLevel].begin();posEdge!=Edges[StartLevel].end();posEdge++)
	{
		ThisTreeIndex = (*posEdge).TreeIndex;
		ThisNodeIndex = (*posEdge).NodeIndex;
		EdgeID = (*posEdge).EdgeID;
		ChildEdgeMidPoint(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EdgeID, ChildVertex);
		*ChildVertex = *ChildVertex + (F_e(ThisTreeIndex, ThisNodeIndex, EdgeID)/2);
	}
	/*------------------------------Lift 6------------------------------------------*/
	for (posVertex=Vertices[StartLevel].begin();posVertex!=Vertices[StartLevel].end();posVertex++)
	{
		ThisTreeIndex = (*posVertex).TreeIndex;
		ThisNodeIndex = (*posVertex).NodeIndex;
		PointPosition = (*posVertex).Position;
		TreeLt.clear();
		NodeIndexLt.clear();
		EdgesLt.clear();
		switch(PointPosition)
		{
		case TL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, NORTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopLeft =*(*ChildNode).TopLeft -	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).TopLeft = *(*ChildNode).TopLeft + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		case TR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getTRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, EAST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).TopRigh =*(*ChildNode).TopRigh -	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).TopRigh = *(*ChildNode).TopRigh + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		case BL:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBLChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, WEST, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotLeft =*(*ChildNode).BotLeft -	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).BotLeft = *(*ChildNode).BotLeft + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		case BR:
			ChildNodeIndex = ForestMesh[ThisTreeIndex].tree->getBRChildNode(ThisNodeIndex);
			ChildNode = ForestMesh[ThisTreeIndex].tree->getNode(ChildNodeIndex);
			Valence = Find1Ring(&ForestMesh[ThisTreeIndex], ThisNodeIndex, SOUTH, TreeLt, NodeIndexLt, EdgesLt);
			*(*ChildNode).BotRigh =*(*ChildNode).BotRigh -	F_v(Valence, TreeLt, NodeIndexLt, EdgesLt)/4;
			*(*ChildNode).BotRigh = *(*ChildNode).BotRigh + E_v(Valence, TreeLt, NodeIndexLt, EdgesLt);
			break;
		}		
	}
}

void Surface::Inv_DWT(int StartLevel, int EndLevel)
{
	for (int i=StartLevel;i<EndLevel;i++)
	{
		Inv_DWT(i);
	}
}

void Surface::MeshOut(const string& filename, int Level)
{
	int i,j;
	VertexList::iterator pos;
	PointSet TempSet;
	Mynode *QualNode;
	Point3D *Vert;
//	Vertax VertThis;

	// open output file
	ofstream file (filename.c_str());

	// file opened?
	if (! file) {
		// NO, abort program
		cerr << "can't open output file \"" << filename << "\""
			<< endl;
		exit (EXIT_FAILURE);
	}

	file <<"OFF"<<endl;
	file <<this->Vertices[Level].size()<<" "<<NumBaseMesh*pow(4.0,Level)*2<<" 0"<<endl;
	i=0;
	for (pos=Vertices[Level].begin(); pos!=Vertices[Level].end(); ++pos) 
	{
		QualNode=(ForestMesh[(*pos).TreeIndex].tree)->getNode((*pos).NodeIndex);
		switch((*pos).Position)
		{
		case TR:
			Vert=(*QualNode).TopRigh;
			break;
		case TL:
			Vert=(*QualNode).TopLeft;
			break;
		case BR:
			Vert=(*QualNode).BotRigh;
		    break;
		case BL:
			Vert=(*QualNode).BotLeft;
		    break;
		default:
		    break;
		}
		(*Vert).Index=i;
		TempSet.insert(Vert);
		i++;
					
		file <<showpoint<<fixed<<(*Vert).x<<" "<<(*Vert).y<<" "<<(*Vert).z<<endl;
	}

	int EndIndex = (pow(4.0, Level+1)-1)/3-1;
	int StartIndex = EndIndex - pow(pow(2.0, Level), 2) + 1;
	for (i=0; i<NumBaseMesh; i++)
	{
		for (j=StartIndex; j<=EndIndex; j++) 
		{
			file<<"3 ";
			file<<(*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->TopLeft))).Index<<" ";
			file<<(*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->BotLeft))).Index<<" ";
			file<<(*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->BotRigh))).Index<<endl; 

			file<<"3 ";
			file<<(*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->BotRigh))).Index<<" ";
			file<<(*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->TopRigh))).Index<<" ";
			file<<(*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->TopLeft))).Index<<endl;
		}		
	}	
}

Point3D Surface::V_f(int treeIndx, int NodeIndx)
{
	Mynode V;
	V = *(ForestMesh[treeIndx].tree->getNode(NodeIndx));
	return (*V.BotLeft+*V.BotRigh+*V.TopLeft+*V.TopRigh)/4;
}

Point3D Surface::V_fChildLevel(int treeIndx, int NodeIndx)
{
	Mynode* TLChild;
	Mynode* BRChild;
	Mynode* TRChild;
	Mynode* BLChild;


	//index of the children of this node
	int TLChildInd;
	int BRChildInd;
	int TRChildInd;
	int BLChildInd;

	TLChildInd = this->ForestMesh[treeIndx].tree->getTLChildNode(NodeIndx);
	BRChildInd = this->ForestMesh[treeIndx].tree->getBRChildNode(NodeIndx);
	TRChildInd = this->ForestMesh[treeIndx].tree->getTRChildNode(NodeIndx);
	BLChildInd = this->ForestMesh[treeIndx].tree->getBLChildNode(NodeIndx);

	TLChild = this->ForestMesh[treeIndx].tree->getNode(TLChildInd);	
	BRChild = this->ForestMesh[treeIndx].tree->getNode(BRChildInd);
	TRChild = this->ForestMesh[treeIndx].tree->getNode(TRChildInd);	
	BLChild = this->ForestMesh[treeIndx].tree->getNode(BLChildInd);

	return (*(TLChild->TopLeft) + *(BLChild->BotLeft) + *(BRChild->BotRigh) + *(TRChild->TopRigh))/4;
}

Point3D Surface::E_f(int treeIndx, int NodeIndx)
{
	Mynode* TLChild;
	Mynode* BRChild;

	//index of the children of this node
	int TLChildInd;
	int BRChildInd;
	
	TLChildInd = this->ForestMesh[treeIndx].tree->getTLChildNode(NodeIndx);
	BRChildInd = this->ForestMesh[treeIndx].tree->getBRChildNode(NodeIndx);
	
	TLChild = this->ForestMesh[treeIndx].tree->getNode(TLChildInd);	
	BRChild = this->ForestMesh[treeIndx].tree->getNode(BRChildInd);
	
	return (*(TLChild->TopRigh) + *(TLChild->BotLeft) + *(BRChild->BotLeft) + *(BRChild->TopRigh))/4;
}

Point3D Surface::V_e(int treeIndx, int NodeIndx, char EdgeID)
{
	Mynode V;
	Point3D P;
	V = *(ForestMesh[treeIndx].tree->getNode(NodeIndx));

	switch(EdgeID) 
	{
	case NORTH:
		P = (*V.TopLeft+*V.TopRigh)/2; 
		break;
	case SOUTH:
		P = (*V.BotLeft+*V.BotRigh)/2;
		break;
	case WEST:
		P = (*V.TopLeft+*V.BotLeft)/2;
		break;
	case EAST:
		P = (*V.TopRigh+*V.BotRigh)/2;
		break;
	}
	return P;
}

Point3D Surface::F_e(int treeIndx, int NodeIndx, char EdgeID)
{
	int NeighNodeIndex; //index of the neighbor node;
	int ThisNodeIndex = NodeIndx;  //index of this node;
	QTreeR* NeighTreeR = NULL; //pointer to the neighbor tree root;
	QTreeR* ThisTreeR = &(ForestMesh[treeIndx]);

	char ThisEdgeID = EdgeID;
	char NeighEdgeID;
	////pointer to the node contain point This_F_Point, one of the children of NodeIndx;
	Mynode* ThisChildNode = NULL;
	Mynode* NeighChildNode = NULL;

	Point3D* This_F_Point = NULL;
	Point3D* Neigh_F_Point = NULL;	

	ThisChildNode = ThisTreeR->tree->getNode(ThisTreeR->tree->getTLChildNode(ThisNodeIndex));
	This_F_Point = ThisChildNode->BotRigh;

	FindEdgeNeighbor(ThisTreeR, ThisNodeIndex, ThisEdgeID, NeighTreeR, NeighNodeIndex, NeighEdgeID);
	NeighChildNode = NeighTreeR->tree->getNode(NeighTreeR->tree->getTLChildNode(NeighNodeIndex));
	Neigh_F_Point = NeighChildNode->BotRigh;	
	return (*This_F_Point+*Neigh_F_Point)/2;
}

void Surface::FindEdgeNeighbor(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID, QTreeR* &NeighTreeR, int &NeighNodeIndex, char &NeighEdgeID)
{
	char NeighTreeID;  //its location relative to its neighbor

	int coordThis[2], coordNeigh[2];
	int levelNode;

	ThisTreeR->tree->getNodeCoordinates(ThisNodeIndex, coordThis, &levelNode);
	
	switch(ThisTreeR->tree->nodeLocation(ThisNodeIndex)) 
	{
	case _INSIDE_NODE:	

		NeighTreeR = ThisTreeR;

		switch(ThisEdgeID)
		{
		case NORTH:
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _NORTH);
			NeighEdgeID = SOUTH;
			break;

		case SOUTH:
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _SOUTH);
			NeighEdgeID = NORTH;
			
			break;
		case EAST:
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _EAST);
			NeighEdgeID = WEST;

			break;
		case WEST:
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _WEST);
			NeighEdgeID = EAST;

			break;
		}
		break;
	case _TOP_BORDER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR =  ThisTreeR->NeibNorth;

			NeighTreeID = ThisTreeR->ID_NeibNorth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}			

			break;
		case SOUTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _SOUTH);
			NeighEdgeID = NORTH;
			
			break;
		case WEST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _WEST);
			NeighEdgeID = EAST;

			break;
		case EAST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _EAST);
			NeighEdgeID = WEST;

			break;
		}
		break;
	case _BOT_BORDER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _NORTH);
			NeighEdgeID = SOUTH;

			break;
		case SOUTH:
			NeighTreeR =  ThisTreeR->NeibSouth;

			NeighTreeID = ThisTreeR->ID_NeibSouth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case EAST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _EAST);
			NeighEdgeID = WEST;

			break;
		case WEST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _WEST);
			NeighEdgeID = EAST;

			break;
		}
		break;
	case _LEFT_BORDER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _NORTH);
			NeighEdgeID = SOUTH;

			break;
		case SOUTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _SOUTH);
			NeighEdgeID = NORTH;

			break;
		case EAST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _EAST);
			NeighEdgeID = WEST;

			break;
		case WEST:
			NeighTreeR =  ThisTreeR->NeibWest;

			NeighTreeID = ThisTreeR->ID_NeibWest;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		}
		break;
	case _RIGHT_BORDER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _NORTH);
			NeighEdgeID = SOUTH;

			break;
		case SOUTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _SOUTH);
			NeighEdgeID = NORTH;

			break;
		case EAST:
			NeighTreeR =  ThisTreeR->NeibEast;

			NeighTreeID = ThisTreeR->ID_NeibEast;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case WEST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _WEST);
			NeighEdgeID = EAST;

			break;
		}
		break;
	case _TR_CORNER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR =  ThisTreeR->NeibNorth;

			NeighTreeID = ThisTreeR->ID_NeibNorth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case SOUTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _SOUTH);
			NeighEdgeID = NORTH;

			break;
		case EAST:
			NeighTreeR =  ThisTreeR->NeibEast;

			NeighTreeID = ThisTreeR->ID_NeibEast;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case WEST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _WEST);
			NeighEdgeID = EAST;

			break;
		}
		break;
	case _TL_CORNER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR =  ThisTreeR->NeibNorth;

			NeighTreeID = ThisTreeR->ID_NeibNorth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case SOUTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _SOUTH);
			NeighEdgeID = NORTH;

			break;
		case EAST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _EAST);
			NeighEdgeID = WEST;

			break;
		case WEST:
			NeighTreeR =  ThisTreeR->NeibWest;

			NeighTreeID = ThisTreeR->ID_NeibWest;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		}
		break;
	case _BR_CORNER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _NORTH);
			NeighEdgeID = SOUTH;

			break;
		case SOUTH:
			NeighTreeR =  ThisTreeR->NeibSouth;

			NeighTreeID = ThisTreeR->ID_NeibSouth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case EAST:
			NeighTreeR =  ThisTreeR->NeibEast;

			NeighTreeID = ThisTreeR->ID_NeibEast;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case WEST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _WEST);
			NeighEdgeID = EAST;

			break;
		}
		break;
	case _BL_CORNER_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _NORTH);
			NeighEdgeID = SOUTH;

			break;
		case SOUTH:
			NeighTreeR =  ThisTreeR->NeibSouth;

			NeighTreeID = ThisTreeR->ID_NeibSouth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case EAST:
			NeighTreeR = ThisTreeR;
			NeighNodeIndex = ThisTreeR->tree->getNeighbor(ThisNodeIndex, _EAST);
			NeighEdgeID = WEST;

			break;
		case WEST:
			NeighTreeR =  ThisTreeR->NeibWest;

			NeighTreeID = ThisTreeR->ID_NeibWest;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		}
		break;
	case _ROOT_NODE:
		switch(ThisEdgeID)
		{
		case NORTH:
			NeighTreeR =  ThisTreeR->NeibNorth;

			NeighTreeID = ThisTreeR->ID_NeibNorth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case SOUTH:
			NeighTreeR =  ThisTreeR->NeibSouth;

			NeighTreeID = ThisTreeR->ID_NeibSouth;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[0];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[0];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[0];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case EAST:
			NeighTreeR =  ThisTreeR->NeibEast;

			NeighTreeID = ThisTreeR->ID_NeibEast;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		case WEST:
			NeighTreeR =  ThisTreeR->NeibWest;

			NeighTreeID = ThisTreeR->ID_NeibWest;

			switch(NeighTreeID) 
			{
			case NORTH:
				coordNeigh[1] = 0;
				coordNeigh[0] = coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = NORTH;

				break;
			case SOUTH:
				coordNeigh[1] = pow(2.0, levelNode) - 1;
				coordNeigh[0] = pow(2.0, levelNode) - 1 - coordThis[1];
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = SOUTH;

				break;
			case WEST:
				coordNeigh[1] = pow(2.0, levelNode) - 1 - coordThis[1];
				coordNeigh[0] = 0;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = WEST;

				break;
			case EAST:
				coordNeigh[1] = coordThis[1];
				coordNeigh[0] = pow(2.0, levelNode) - 1;
				NeighNodeIndex = NeighTreeR->tree->getNodeIndex(coordNeigh, levelNode);
				NeighEdgeID = EAST;

				break;
			}

			break;
		}	   
		break;
	default:
		cout<<"error!"<<endl;
	} 
}

void Surface::EdgePoint(QTreeR* TreeR, int NodeIndex, char EdgeID, Point3D** &Start, Point3D** &End)
{
	Mynode* Quadangle = TreeR->tree->getNode(NodeIndex);

	switch(EdgeID) 
	{
	case NORTH:
		Start = &(Quadangle->TopRigh);
		End = &(Quadangle->TopLeft);
		break;
	case SOUTH:
		Start = &(Quadangle->BotLeft);
		End = &(Quadangle->BotRigh);
		break;
	case WEST:
		Start = &(Quadangle->TopLeft);
		End = &(Quadangle->BotLeft);
		break;
	case EAST:
		Start = &(Quadangle->BotRigh);
		End = &(Quadangle->TopRigh);
		break;
	}
}

Point3D Surface::ChildEdgeMidPoint(QTreeR* TreeR, int NodeIndex, char EdgeID)
{
	Mynode* Quadangle = NULL;
	Point3D Pt;
	int ChildNodeIndex;

	switch(EdgeID) 
	{
	case NORTH:
		ChildNodeIndex = TreeR->tree->getTLChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = *(Quadangle->TopRigh);
		break;
	case SOUTH:
		ChildNodeIndex = TreeR->tree->getBLChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = *(Quadangle->BotRigh);
		break;
	case WEST:
		ChildNodeIndex = TreeR->tree->getTLChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = *(Quadangle->BotLeft);
		break;
	case EAST:
		ChildNodeIndex = TreeR->tree->getBRChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = *(Quadangle->TopRigh);
		break;
	} 
	return Pt;
}

void Surface::ChildEdgeMidPoint(QTreeR* TreeR, int NodeIndex, char EdgeID, Point3D* &Pt)
{
	Mynode* Quadangle = NULL;
	int ChildNodeIndex;

	switch(EdgeID) 
	{
	case NORTH:
		ChildNodeIndex = TreeR->tree->getTLChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = (Quadangle->TopRigh);
		break;
	case SOUTH:
		ChildNodeIndex = TreeR->tree->getBLChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = (Quadangle->BotRigh);
		break;
	case WEST:
		ChildNodeIndex = TreeR->tree->getTLChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = (Quadangle->BotLeft);
		break;
	case EAST:
		ChildNodeIndex = TreeR->tree->getBRChildNode(NodeIndex);
		Quadangle = TreeR->tree->getNode(ChildNodeIndex);
		Pt = (Quadangle->TopRigh);
		break;
	} 
}

Point3D Surface::V_v(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID)
{
	int valence;
	Point3D Vertex;

	list<QTreeR*> TreeList;
	list<int> NodeIndexList; 
	list<char> EdgesList;
	list<QTreeR*>::iterator posTree;
	list<int>::iterator posNode;
	list<char>::iterator posEdge;

	Point3D** Start;
	Point3D** End;

	valence = Find1Ring(ThisTreeR, ThisNodeIndex, ThisEdgeID, TreeList, NodeIndexList, EdgesList);

	posTree = TreeList.begin();
	posNode = NodeIndexList.begin();
	posEdge = EdgesList.begin();

	for (int i=0;i<valence;i++)
	{
		EdgePoint(*posTree, *posNode, *posEdge, Start, End);
		Vertex = Vertex + (**Start);
		posTree++;
		posNode++;
		posEdge++;
	}
	return Vertex/valence;
}

//Point3D Surface::E_v(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID)
//{
//	
//}

Point3D Surface::F_v(QTreeR* ThisTreeR, int ThisNodeIndex, char ThisEdgeID)
{
	int valence;
	Point3D Vertex;
	Mynode* Quadrangle;
	int ChildNodeIndex;

	list<QTreeR*> TreeList;
	list<int> NodeIndexList; 
	list<char> EdgesList;
	list<QTreeR*>::iterator posTree;
	list<int>::iterator posNode;
	list<char>::iterator posEdge;

	valence = Find1Ring(ThisTreeR, ThisNodeIndex, ThisEdgeID, TreeList, NodeIndexList, EdgesList);

	posTree = TreeList.begin();
	posNode = NodeIndexList.begin();
	posEdge = EdgesList.begin();

	for (int i=0;i<valence;i++)
	{
		ChildNodeIndex = (*posTree)->tree->getTLChildNode(*posNode);
		Quadrangle = (*posTree)->tree->getNode(ChildNodeIndex);
		Vertex = Vertex + *(Quadrangle->BotRigh);
		posTree++;
		posNode++;
		posEdge++;
	}
	return Vertex/valence;
}

int Surface::Find1Ring(QTreeR* TreeR, int NodeIndex, char EdgeID, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList)
{
	int Valence = 0;

	QTreeR* CurrentTree = TreeR;
	int CurrentNodeIndex = NodeIndex;
	char CurrentEdgeID = EdgeID;

	QTreeR* NextTree = NULL;
	int NextNodeIndex = 0;

	char NeighEdgeID;

	TreeList.push_back(TreeR);
	NodeIndexList.push_back(NodeIndex);
	EdgesList.push_back(EdgeID);
	Valence++;

	FindEdgeNeighbor(TreeR, NodeIndex, EdgeID, NextTree, NextNodeIndex, NeighEdgeID);

	while (!((TreeR==NextTree)&&(NodeIndex==NextNodeIndex)))
	{
		TreeList.push_back(NextTree);
		NodeIndexList.push_back(NextNodeIndex);
		Valence++;
		
		switch(NeighEdgeID)
		{
		case NORTH:
			EdgesList.push_back(EAST);
			CurrentTree = NextTree;
			CurrentNodeIndex = NextNodeIndex;
			CurrentEdgeID = EAST;
			FindEdgeNeighbor(CurrentTree, CurrentNodeIndex, CurrentEdgeID, NextTree, NextNodeIndex, NeighEdgeID);
			break;
		case SOUTH:
			EdgesList.push_back(WEST);
			CurrentTree = NextTree;
			CurrentNodeIndex = NextNodeIndex;
			CurrentEdgeID = WEST;
			FindEdgeNeighbor(CurrentTree, CurrentNodeIndex, CurrentEdgeID, NextTree, NextNodeIndex, NeighEdgeID);
			break;
		case WEST:
			EdgesList.push_back(NORTH);
			CurrentTree = NextTree;
			CurrentNodeIndex = NextNodeIndex;
			CurrentEdgeID = NORTH;
			FindEdgeNeighbor(CurrentTree, CurrentNodeIndex, CurrentEdgeID, NextTree, NextNodeIndex, NeighEdgeID);
			break;
		case EAST:
			EdgesList.push_back(SOUTH);
			CurrentTree = NextTree;
			CurrentNodeIndex = NextNodeIndex;
			CurrentEdgeID = SOUTH;
			FindEdgeNeighbor(CurrentTree, CurrentNodeIndex, CurrentEdgeID, NextTree, NextNodeIndex, NeighEdgeID);
			break;
		}
	}
	return Valence;	
}

Point3D Surface::V_v(int Valence, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList)
	{
		Point3D Vertex;
		
		list<QTreeR*>::iterator posTree;
		list<int>::iterator posNode;
		list<char>::iterator posEdge;

		Point3D** Start;
		Point3D** End;

		posTree = TreeList.begin();
		posNode = NodeIndexList.begin();
		posEdge = EdgesList.begin();

		for (int i=0;i<Valence;i++)
		{
			EdgePoint(*posTree, *posNode, *posEdge, Start, End);
			Vertex = Vertex + (**Start);
			posTree++;
			posNode++;
			posEdge++;
		}
		return Vertex/Valence;
	}

	Point3D Surface::E_v(int Valence, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList)
	{
		Point3D Vertex;
		//Mynode* Quadrangle;
		//int ChildNodeIndex;

		list<QTreeR*>::iterator posTree;
		list<int>::iterator posNode;
		list<char>::iterator posEdge;

		posTree = TreeList.begin();
		posNode = NodeIndexList.begin();
		posEdge = EdgesList.begin();

		for (int i=0;i<Valence;i++)
		{
			Vertex = Vertex + ChildEdgeMidPoint(*posTree, *posNode, *posEdge);
			posTree++;
			posNode++;
			posEdge++;
		}
		return Vertex/Valence;
	}

	Point3D Surface::F_v(int Valence, list<QTreeR*> &TreeList, list<int> &NodeIndexList, list<char> &EdgesList)
	{
		Point3D Vertex;
		Mynode* Quadrangle;
		int ChildNodeIndex;

		list<QTreeR*>::iterator posTree;
		list<int>::iterator posNode;
		list<char>::iterator posEdge;

		posTree = TreeList.begin();
		posNode = NodeIndexList.begin();
		posEdge = EdgesList.begin();

		for (int i=0;i<Valence;i++)
		{
			ChildNodeIndex = (*posTree)->tree->getTLChildNode(*posNode);
			Quadrangle = (*posTree)->tree->getNode(ChildNodeIndex);
			Vertex = Vertex + *(Quadrangle->BotRigh);
			posTree++;
			posNode++;
			posEdge++;
		}
		return Vertex/Valence;
	}

	Point3D Surface::V_e_ChildLevel(int treeIndx, int NodeIndx, char EdgeID)
	{
		Mynode C1, C2;
		Point3D P;

		switch(EdgeID) 
		{
		case NORTH:
			C1 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getTLChildNode(NodeIndx));
			C2 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getTRChildNode(NodeIndx));
			P = (*C1.TopLeft+*C2.TopRigh)/2; 
			break;
		case SOUTH:
			C1 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getBLChildNode(NodeIndx));
			C2 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getBRChildNode(NodeIndx));
			P = (*C1.BotLeft+*C2.BotRigh)/2;
			break;
		case WEST:
			C1 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getTLChildNode(NodeIndx));
			C2 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getBLChildNode(NodeIndx));
			P = (*C1.TopLeft+*C2.BotLeft)/2;
			break;
		case EAST:
			C1 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getTRChildNode(NodeIndx));
			C2 = *ForestMesh[treeIndx].tree->getNode(ForestMesh[treeIndx].tree->getBRChildNode(NodeIndx));
			P = (*C1.TopRigh+*C2.BotRigh)/2;
			break;
		}
		return P;
	}

	void Surface::SaveObj(const string& filename)
	{
		int i, j;
		Point3D* Point;

		// open output file.................................................
		ofstream file (filename.c_str());

		// file opened?......................................................
		if (! file) {
			// NO, abort program
			cerr << "can't open output file \"" << filename << "\""
				<< endl;
			exit (EXIT_FAILURE);
		}

		//write the file header...............................................
		file <<"SUR"<<endl;
		file <<"NumBaseMesh "<<this->NumBaseMesh<<" "<<"NumLevel "<<this->NumLevel<<endl;
		file <<"TransX_mean "<<showpoint<<fixed<<this->TransX_mean<<" "<<"TransX_std "<<this->TransX_std<<endl; 
		file <<"TransY_mean "<<this->TransY_mean<<" "<<"TransY_std "<<this->TransY_std<<endl;
		file <<"TransZ_mean "<<this->TransZ_mean<<" "<<"TransZ_std "<<this->TransZ_std<<endl;
		file <<"ScaleX_mean "<<this->ScaleX_mean<<" "<<"ScaleX_std "<<this->ScaleX_std<<endl;
		file <<"ScaleY_mean "<<this->ScaleY_mean<<" "<<"ScaleY_std "<<this->ScaleY_std<<endl;
		file <<"ScaleZ_mean "<<this->ScaleZ_mean<<" "<<"ScaleZ_std "<<this->ScaleZ_std<<endl;
		file <<"RotationPhi_mean "<<this->RotationPhi_mean<<" "<<"RotationPhi_std "<<this->RotationPhi_std<<endl;
		file <<"RotationTheta_mean "<<this->RotationTheta_mean<<" "<<"RotationTheta_std "<<this->RotationTheta_std<<endl;
		file <<"RotationTau_mean "<<this->RotationTau_mean<<" "<<"RotationTau_std "<<this->RotationTau_std<<endl;

		//begin to write the coefficients...........................................
		//"x y z MeanX StdX MeanY StdY MeanZ StdZ" is the format....................
		int EndIndex = (int)((pow(4.0, this->NumLevel)-1)/3-1);
		int StartIndex = 0;

		for (i=0; i<this->NumBaseMesh; i++)
		{
			for (j=StartIndex; j<=EndIndex; j++) 
			{
				Point=this->ForestMesh[i].tree->getNode(j)->TopLeft;
				file <<Point->x<<" ";
				file <<Point->y<<" ";
				file <<Point->z<<" ";
				file <<Point->MeanX<<" ";
				file <<Point->StdX<<" ";
				file <<Point->MeanY<<" ";
				file <<Point->StdY<<" ";
				file <<Point->MeanZ<<" ";
				file <<Point->StdZ<<endl;

				Point=this->ForestMesh[i].tree->getNode(j)->TopRigh;
				file <<Point->x<<" ";
				file <<Point->y<<" ";
				file <<Point->z<<" ";
				file <<Point->MeanX<<" ";
				file <<Point->StdX<<" ";
				file <<Point->MeanY<<" ";
				file <<Point->StdY<<" ";
				file <<Point->MeanZ<<" ";
				file <<Point->StdZ<<endl;

				Point=this->ForestMesh[i].tree->getNode(j)->BotLeft;
				file <<Point->x<<" ";
				file <<Point->y<<" ";
				file <<Point->z<<" ";
				file <<Point->MeanX<<" ";
				file <<Point->StdX<<" ";
				file <<Point->MeanY<<" ";
				file <<Point->StdY<<" ";
				file <<Point->MeanZ<<" ";
				file <<Point->StdZ<<endl;

				Point=this->ForestMesh[i].tree->getNode(j)->BotRigh;
				file <<Point->x<<" ";
				file <<Point->y<<" ";
				file <<Point->z<<" ";
				file <<Point->MeanX<<" ";
				file <<Point->StdX<<" ";
				file <<Point->MeanY<<" ";
				file <<Point->StdY<<" ";
				file <<Point->MeanZ<<" ";
				file <<Point->StdZ<<endl;
			}		
		}		
	}

	void Surface::ReadObj(const char* filename)
	{
		int i,j;
		FILE *fp;
		char head[10];
		int NbaseMesh=0, NLevel=0;
		//float data;
		Point3D* Point;

		fp=fopen(filename,"r");
		if (fp==NULL)
		{
			printf("cannot open file");
			exit(0);
		}

		fscanf(fp,"%s\n",head); //read in the file header 'SUR' 
		fscanf(fp,"%s %d %s %d\n",head,&NbaseMesh,head,&NLevel);
		if (((this->NumBaseMesh)!=NbaseMesh)||((this->NumLevel)!=NLevel))
		{
			printf("The dimension of the surface object doesn't match.");
			exit(0);
		}
		//begin read in data...............................................
		fscanf(fp,"%s %f %s %f\n",head,&(this->TransX_mean),head,&(this->TransX_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->TransY_mean),head,&(this->TransY_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->TransZ_mean),head,&(this->TransZ_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->ScaleX_mean),head,&(this->ScaleX_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->ScaleY_mean),head,&(this->ScaleY_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->ScaleZ_mean),head,&(this->ScaleZ_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->RotationPhi_mean),head,&(this->RotationPhi_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->RotationTheta_mean),head,&(this->RotationTheta_std));
		fscanf(fp,"%s %f %s %f\n",head,&(this->RotationTau_mean),head,&(this->RotationTau_std));
		//begin to read the coefficients...........................................
		//"x y z MeanX StdX MeanY StdY MeanZ StdZ" is the format....................
		int EndIndex = (pow(4.0, this->NumLevel)-1)/3-1;
		int StartIndex = 0;

		for (i=0; i<this->NumBaseMesh; i++)
		{
			for (j=StartIndex; j<=EndIndex; j++) 
			{
				Point=this->ForestMesh[i].tree->getNode(j)->TopLeft;
				fscanf(fp,"%f %f %f %f %f %f %f %f %f \n",&(Point->x),&(Point->y),&(Point->z),&(Point->MeanX),&(Point->StdX),&(Point->MeanY),&(Point->StdY),&(Point->MeanZ),&(Point->StdZ));

				Point=this->ForestMesh[i].tree->getNode(j)->TopRigh;
				fscanf(fp,"%f %f %f %f %f %f %f %f %f \n",&(Point->x),&(Point->y),&(Point->z),&(Point->MeanX),&(Point->StdX),&(Point->MeanY),&(Point->StdY),&(Point->MeanZ),&(Point->StdZ));

				Point=this->ForestMesh[i].tree->getNode(j)->BotLeft;
				fscanf(fp,"%f %f %f %f %f %f %f %f %f \n",&(Point->x),&(Point->y),&(Point->z),&(Point->MeanX),&(Point->StdX),&(Point->MeanY),&(Point->StdY),&(Point->MeanZ),&(Point->StdZ));

				Point=this->ForestMesh[i].tree->getNode(j)->BotRigh;
				fscanf(fp,"%f %f %f %f %f %f %f %f %f \n",&(Point->x),&(Point->y),&(Point->z),&(Point->MeanX),&(Point->StdX),&(Point->MeanY),&(Point->StdY),&(Point->MeanZ),&(Point->StdZ));
			}		
		}	

		fclose(fp); /*close the file*/
	}

	void Surface::MeshIn(const char* filename, int Level)
	{
		VertexList::iterator pos;
		Mynode *QualNode;
		Point3D *Vert;
		int numVertices;
		int numFacets;
		int temp;
				
		//skip the file header
		char buffer[256];
		ifstream myfile (filename);
		while (! myfile.eof() )
		{
			myfile.getline (buffer,100);
			if(strstr(buffer,"OFF"))
			{
				break;
			}
		}

		myfile>>numVertices;
		myfile>>numFacets;
		myfile>>temp;

		//begin to read in the vertex
		for (pos=Vertices[Level].begin(); pos!=Vertices[Level].end(); ++pos) 
		{
			QualNode=(ForestMesh[(*pos).TreeIndex].tree)->getNode((*pos).NodeIndex);
			switch((*pos).Position)
			{
			case TR:
				Vert=(*QualNode).TopRigh;
				break;
			case TL:
				Vert=(*QualNode).TopLeft;
				break;
			case BR:
				Vert=(*QualNode).BotRigh;
				break;
			case BL:
				Vert=(*QualNode).BotLeft;
				break;
			default:
				break;
			}
			
			myfile>>(*Vert).x;		
			myfile>>(*Vert).y;
			myfile>>(*Vert).z;
		}
	}

	void Surface::SaveCof(const char* filename)
	{
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;

		//write the file header
		ofstream myfile (filename);
		myfile.precision(12);
		myfile<<"COF\n";

		//begin to read in the vertex
		for (Level=0;Level<this->NumLevel;Level++)
		{
			myfile<<"level "<<Level<<" ";
			myfile<<Coefficients[Level].size()<<"\n";
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				myfile<<(*Vert).x<<" ";		
				myfile<<(*Vert).y<<" ";
				myfile<<(*Vert).z<<"\n";
			}		
		}		
	}

	void Surface::ReadCof(const char* filename)
	{
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;
		string buffer;

		//read the file header
		ifstream myfile (filename);
		getline(myfile,buffer) ;

		//begin to read in the vertex
		if (myfile.is_open())
		{
			for (Level=0;Level<this->NumLevel;Level++)
			{
				getline(myfile,buffer);
				for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
				{
					Vert=(*pos).Coeffi;

					getline(myfile,buffer);
					stringstream(buffer)>>(*Vert).x>>(*Vert).y>>(*Vert).z;
					//cout<<(*Vert).x<<" "<<(*Vert).y<<" "<<(*Vert).z<<endl;
				}		
			}

			myfile.close();
		}
		else 
		{
			cout << "Unable to open file"; 
		}
	}

	void Surface::MeshToITK(TriangleMeshSourceType::Pointer &ITKMeshSource, int NumofLevel)
	{
		int i,j;
		int EndIndex = (pow(4.0, NumofLevel)-1)/3-1;
		int StartIndex = EndIndex - pow(pow(2.0, NumofLevel-1), 2) + 1;

		for (i=0; i<NumBaseMesh; i++)
		{
			for (j=StartIndex; j<=EndIndex; j++) 
			{
				ITKMeshSource->AddTriangle(
					ITKMeshSource->AddPoint(this->ForestMesh[i].tree->getNode(j)->TopLeft->x, this->ForestMesh[i].tree->getNode(j)->TopLeft->y, this->ForestMesh[i].tree->getNode(j)->TopLeft->z),
					ITKMeshSource->AddPoint(this->ForestMesh[i].tree->getNode(j)->BotLeft->x, this->ForestMesh[i].tree->getNode(j)->BotLeft->y, this->ForestMesh[i].tree->getNode(j)->BotLeft->z),
					ITKMeshSource->AddPoint(this->ForestMesh[i].tree->getNode(j)->BotRigh->x, this->ForestMesh[i].tree->getNode(j)->BotRigh->y, this->ForestMesh[i].tree->getNode(j)->BotRigh->z) 
										   );

				ITKMeshSource->AddTriangle(
					ITKMeshSource->AddPoint(this->ForestMesh[i].tree->getNode(j)->TopLeft->x, this->ForestMesh[i].tree->getNode(j)->TopLeft->y, this->ForestMesh[i].tree->getNode(j)->TopLeft->z),
					ITKMeshSource->AddPoint(this->ForestMesh[i].tree->getNode(j)->BotRigh->x, this->ForestMesh[i].tree->getNode(j)->BotRigh->y, this->ForestMesh[i].tree->getNode(j)->BotRigh->z),
					ITKMeshSource->AddPoint(this->ForestMesh[i].tree->getNode(j)->TopRigh->x, this->ForestMesh[i].tree->getNode(j)->TopRigh->y, this->ForestMesh[i].tree->getNode(j)->TopRigh->z) 
					                      );
			}		
		}	
	}

	void Surface::ReadPrior(const char* filename)
	{
		int Level,NumLev;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;
		string buffer;
		string temp;
		stringstream ss (stringstream::in | stringstream::out);

		//read the file header
		ifstream myfile (filename);
		getline(myfile,buffer);
		ss.clear();
		ss << buffer;
		ss>>temp>>this->NumLevPrior>>this->NumSample;

		NumLev = this->NumLevPrior;

		//begin to read in the vertex
		if (myfile.is_open())
		{
			for (Level=0;Level<NumLev;Level++)
			{
				getline(myfile,buffer);
				for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
				{
					Vert=(*pos).Coeffi;

					getline(myfile,buffer);
					ss.clear();
					ss<<buffer;

					//read in the pc
					ss>>temp>>(*Vert).PriorPC[0][0]>>(*Vert).PriorPC[0][1]>>(*Vert).PriorPC[0][2];
					ss>>temp>>(*Vert).PriorPC[1][0]>>(*Vert).PriorPC[1][1]>>(*Vert).PriorPC[1][2];
					ss>>temp>>(*Vert).PriorPC[2][0]>>(*Vert).PriorPC[2][1]>>(*Vert).PriorPC[2][2];

					//read in the mean
					ss>>temp>>(*Vert).PriorMEAN[0]>>(*Vert).PriorMEAN[1]>>(*Vert).PriorMEAN[2];
					
					//read in the min score
					ss>>temp>>(*Vert).PriorMINScore[0]>>(*Vert).PriorMINScore[1]>>(*Vert).PriorMINScore[2];
					
					//read in the max score
					ss>>temp>>(*Vert).PriorMAXScore[0]>>(*Vert).PriorMAXScore[1]>>(*Vert).PriorMAXScore[2];

					//read in the std
					ss>>temp>>(*Vert).PriorSTD[0]>>(*Vert).PriorSTD[1]>>(*Vert).PriorSTD[2];

					//read in the scores of samples
					/*for (int i=0; i<NUM_Sample; i++)
					{
						ss>>temp>>(*Vert).Scores[i][0]>>(*Vert).Scores[i][1]>>(*Vert).Scores[i][2];
					}*/
					//cout<<(*Vert).PriorPC[0][0]<<" "<<(*Vert).PriorPC[0][1]<<" "<<(*Vert).PriorPC[0][2]<<endl;//test
				}		
			}

			myfile.close();
		}
		else 
		{
			cout << "Unable to open file"; 
		}
	}

	////////////////function to reconstruct the mean shape from prior model///////////////
	////// Level: the level of the surface to reconstruct (start from '0')
	//////////////////////////////////////////////////////////////////////////////////////
	void Surface::GetMean()
	{
		int NumLev = this->NumLevPrior;
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;
		
		//begin to read in the mean from prior
		for (Level=0;Level<NumLev;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//read in the mean value
				(*Vert).x=(*Vert).PriorMEAN[0];
				(*Vert).y=(*Vert).PriorMEAN[1];
				(*Vert).z=(*Vert).PriorMEAN[2];					
			}		
		}
		
		//write the coefficients in the rest levels to 0
		for (Level=NumLev;Level<this->NumLevel;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//read in the mean value
				(*Vert).x=0.0;
				(*Vert).y=0.0;
				(*Vert).z=0.0;					
			}		
		}
	}

	////////////////function to reconstruct the shape the specified sample///////////////
	////// NumIndex: the ID of the sample to reconstruct (start from '0')
	//////////////////////////////////////////////////////////////////////////////////////
	void Surface::GetSample(int NumIndex)
	{
		int NumLev = this->NumLevPrior;
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;

		//begin to read in the mean from prior
		for (Level=0;Level<NumLev;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//compute the value
				this->ComputeCoefVertax(Vert, Vert->Scores[NumIndex-1]);
			}		
		}

		//write the coefficients in the rest levels to 0
		for (Level=NumLev;Level<this->NumLevel;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//read in the mean value
				(*Vert).x=0.0;
				(*Vert).y=0.0;
				(*Vert).z=0.0;					
			}		
		}
	}

	TriangleMeshType::Pointer Surface::GenerateTriangleMesh(int Level) 
	{
		typedef TriangleMeshType::CellType CellType;
		typedef itk::TriangleCell< CellType > TriangleType;
		typedef CellType::CellAutoPointer CellAutoPointer;

		TriangleMeshType::PointType Point;
		CellAutoPointer Triangle;

		int i,j,k;
		VertexList::iterator pos;
		PointSet TempSet;
		Mynode *QualNode;
		Point3D *Vert;

		TriangleMeshType::Pointer TriangleMesh = TriangleMeshType::New(); 
		//traverse the vertex in mesh at Level
		i=0;
		for (pos=Vertices[Level].begin(); pos!=Vertices[Level].end(); ++pos) 
		{
			QualNode=(ForestMesh[(*pos).TreeIndex].tree)->getNode((*pos).NodeIndex);
			switch((*pos).Position)
			{
			case TR:
				Vert=(*QualNode).TopRigh;
				break;
			case TL:
				Vert=(*QualNode).TopLeft;
				break;
			case BR:
				Vert=(*QualNode).BotRigh;
				break;
			case BL:
				Vert=(*QualNode).BotLeft;
				break;
			default:
				break;
			}
			(*Vert).Index=i;  //set an ID to the point
			TempSet.insert(Vert);			

			//insert the point into the mesh
			Point[0]=(*Vert).x; Point[1]=(*Vert).y; Point[2]=(*Vert).z; 
			TriangleMesh->SetPoint(i,Point);

			i++;
		}

		//j = m_TriangleMesh->GetNumberOfPoints();
		k=0;
		int EndIndex = (pow(4.0, Level+1)-1)/3-1;
		int StartIndex = EndIndex - pow(pow(2.0, Level), 2) + 1;
		for (i=0; i<NumBaseMesh; i++)
		{
			for (j=StartIndex; j<=EndIndex; j++) 
			{
				//add triangle 1
				Triangle.TakeOwnership( new TriangleType ); 
				Triangle->SetPointId( 0, (*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->TopLeft))).Index );
				Triangle->SetPointId( 1, (*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->BotLeft))).Index );
				Triangle->SetPointId( 2, (*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->BotRigh))).Index );
				TriangleMesh->SetCell( k, Triangle );
				k++;

				//add triangle 2
				Triangle.TakeOwnership( new TriangleType ); 
				Triangle->SetPointId( 0, (*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->BotRigh))).Index );
				Triangle->SetPointId( 1, (*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->TopRigh))).Index );
				Triangle->SetPointId( 2, (*(*TempSet.find(this->ForestMesh[i].tree->getNode(j)->TopLeft))).Index );
				TriangleMesh->SetCell( k, Triangle );
				k++;
			}		
		}

		return TriangleMesh;
	}

	void Surface::ComputeCoefVertax(Point3D* &Vert, float Scores[3])
	{
		float temp[3];
		temp[0]=temp[1]=temp[2]=0.0;

		for (int i=0; i<3; i++)
		{
			temp[0]=temp[0]+Vert->PriorPC[i][0]*Scores[i];
			temp[1]=temp[1]+Vert->PriorPC[i][1]*Scores[i];
			temp[2]=temp[2]+Vert->PriorPC[i][2]*Scores[i];
		}

		Vert->x = temp[0]+Vert->PriorMEAN[0];
		Vert->y = temp[1]+Vert->PriorMEAN[1];
		Vert->z = temp[2]+Vert->PriorMEAN[2];
	}

	void Surface::ComputeCoefVertaxBackup(Point3D* &Vert, float Scores[3])
	{
		float temp[3];
		temp[0]=temp[1]=temp[2]=0.0;

		for (int i=0; i<3; i++)
		{
			temp[0]=temp[0]+Vert->PriorPC[i][0]*Scores[i];
			temp[1]=temp[1]+Vert->PriorPC[i][1]*Scores[i];
			temp[2]=temp[2]+Vert->PriorPC[i][2]*Scores[i];
		}

		Vert->X = temp[0]+Vert->PriorMEAN[0];
		Vert->Y = temp[1]+Vert->PriorMEAN[1];
		Vert->Z = temp[2]+Vert->PriorMEAN[2];
		//cout<<Vert->PriorMEAN[2]<<endl;//test
	}

	void Surface::GetBackup(void)
	{
		int NumLev = this->NumLevel;
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;

		//traversal the vertex
		for (Level=0;Level<NumLev;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//get back the value
				Vert->x=Vert->X;
				Vert->y=Vert->Y;
				Vert->z=Vert->Z;
			}		
		}
	}

	void Surface::SetBackup(void)
	{
		int NumLev = this->NumLevel;
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;

		//traversal the vertex
		for (Level=0;Level<NumLev;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//backup the value
				Vert->X=Vert->x;
				Vert->Y=Vert->y;
				Vert->Z=Vert->z;
			}		
		}
	}

	float Surface::ComputeFitness(void)
	{
		
		VertexList::iterator pos;
		Mynode *QualNode;
		Point3D *Vert;
		DeformableModelApplicationBase::VolumeType::IndexType pixelIndex;
		DeformableModelApplicationBase::VolumeType::PixelType pixelValue;
		int Level=this->NumLevel-1;

		float fitness=0;

		for (pos=Vertices[Level].begin(); pos!=Vertices[Level].end(); ++pos) 
		{
			QualNode=(ForestMesh[(*pos).TreeIndex].tree)->getNode((*pos).NodeIndex);
			switch((*pos).Position)
			{
			case TR:
				Vert=(*QualNode).TopRigh;
				break;
			case TL:
				Vert=(*QualNode).TopLeft;
				break;
			case BR:
				Vert=(*QualNode).BotRigh;
				break;
			case BL:
				Vert=(*QualNode).BotLeft;
				break;
			default:
				break;
			}
			
			
			pixelIndex[0] = floor((*Vert).x+0.5); // x position
			pixelIndex[1] = floor((*Vert).y+0.5); // y position
			pixelIndex[2] = floor((*Vert).z+0.5); // z position

			pixelValue = this->Image->GetPixel( pixelIndex );
			if (pixelValue>0.7)
			{
				pixelValue=15;
				//pixelValue=85;
			}	

			fitness=fitness+pixelValue;
			//cout<<pixelIndex[0]<<" "<<pixelIndex[1]<<" "<<pixelIndex[2]<<": "<<pixelValue<<endl;
		}		
		//cout<<"the image term = "<<fitness<<endl;
		return fitness+Alpha*this->PriTerm(); 
		

		//int i,j;
		//Point3D Center;//,Topleft,Topright,Botleft,Botright;
		//int Level=this->NumLevel-1;
		//Mynode *TheNode;
		//DeformableModelApplicationBase::VolumeType::IndexType pixelIndex;
		//DeformableModelApplicationBase::VolumeType::PixelType pixelValue;

		//float fitness=0;
		//float Area=0.0;
		////float im=0.0;
		//
		//int EndIndex = (pow(4.0, Level+1)-1)/3-1;
		//int StartIndex = EndIndex - pow(pow(2.0, Level), 2) + 1;
		//for (i=0; i<NumBaseMesh; i++)
		//{
		//	for (j=StartIndex; j<=EndIndex; j++) 
		//	{
		//		TheNode=this->ForestMesh[i].tree->getNode(j);
		//		Center = ((*(TheNode->TopLeft)) + (*(TheNode->TopRigh)) + (*(TheNode->BotLeft)) + (*(TheNode->BotRigh)))/4.0;
		//		pixelIndex[0] = floor(Center.x+0.5); // x position
		//		pixelIndex[1] = floor(Center.y+0.5); // y position
		//		pixelIndex[2] = floor(Center.z+0.5); // z position
		//		pixelValue = this->Image->GetPixel( pixelIndex );
		//		if (pixelValue>0.6)
		//		{
		//			pixelValue=50;
		//		}				

		//		Area=QuadArea(TheNode);

		//		fitness=fitness+pixelValue*Area;
		//	}		
		//}	
		//return fitness;
	}

	void Surface::SetImage(DeformableModelApplicationBase::VolumeType* ImageToFit)
	{
        this->Image=ImageToFit;
	}

	void Surface::DeformByVertax(int Lev, int StepLen, int StepNum)
	{
		int NumLev = this->NumLevPrior;
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;
		//just for output the interim results----------------start---------------
		char IntrimResultName[50]; 
		//just for output the interim results----------------start---------------
		
		
		float Fitness;               //the current fitness
		float FitnessMin;            //the best fitness achieved
		int   FitnessIndex;          //the step index with best fiting
		float ScoreBest[3];          //the scores corresponding to best fit

		float step[3];
		float scores[3];

		//write the file header
		//ofstream myfile ("test.txt");
		//myfile.precision(12);

		//begin to traverse the vertex
		for (Level=Lev; Level<Lev+1; Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;
				//compute the step length
				step[0] = (Vert->PriorMAXScore[0]-Vert->PriorMINScore[0])/StepLen;
				step[1] = (Vert->PriorMAXScore[1]-Vert->PriorMINScore[1])/StepLen;
				step[2] = (Vert->PriorMAXScore[2]-Vert->PriorMINScore[2])/StepLen;
				/*step[0] = Vert->PriorSTD[0]/StepLen;
				step[1] = Vert->PriorSTD[1]/StepLen;
				step[2] = Vert->PriorSTD[2]/StepLen;*/

				//search the line by step
				//----------  1st pc--------------------------------------------------------
				FitnessMin = 999999; 
				for (int i=-StepNum; i<=StepNum; i++)
				{
					scores[0]=i*step[0];
					scores[1]=0;
					scores[2]=0;
					this->ComputeCoefVertaxBackup(Vert, scores);
					this->GetBackup();
					this->Inv_DWT(0,4);
					this->CC_Subd(4,5);
					//cout<<Vert->X<<" "<<Vert->Y<<" "<<Vert->Z<<": "<<this->ComputeFitness()<<endl;
					Fitness = this->ComputeFitness();
					//myfile<<Fitness<<" ";
					if (Fitness<FitnessMin)
					{
						FitnessMin=Fitness;
						FitnessIndex=i;
						ScoreBest[0]=scores[0];
					}					
				}
				//Set to best fit
				scores[0]=ScoreBest[0];
				scores[1]=0;
				scores[2]=0;
				this->ComputeCoefVertaxBackup(Vert, scores);
				//myfile<<"best fit: "<<FitnessMin<<" best index: "<<FitnessIndex<<endl;
				
				//----------  2nd pc--------------------------------------------------------
				FitnessMin = 999999;
				for (int i=-StepNum; i<=StepNum; i++)
				{
					scores[0]=ScoreBest[0];
					scores[1]=i*step[1];
					scores[2]=0;
					this->ComputeCoefVertaxBackup(Vert, scores);
					this->GetBackup();
					this->Inv_DWT(0,4);
					this->CC_Subd(4,5);
					//cout<<Vert->X<<" "<<Vert->Y<<" "<<Vert->Z<<": "<<this->ComputeFitness()<<endl;
					Fitness = this->ComputeFitness();
					//myfile<<Fitness<<" ";
					if (Fitness<FitnessMin)
					{
						FitnessMin=Fitness;
						FitnessIndex=i;
						ScoreBest[1]=scores[1];
					}					
				}
				//Set to best fit
				scores[0]=ScoreBest[0];
				scores[1]=ScoreBest[1];
				scores[2]=0;
				this->ComputeCoefVertaxBackup(Vert, scores);
				//myfile<<"best fit: "<<FitnessMin<<" best index: "<<FitnessIndex<<endl;
				
				//----------  3rd pc--------------------------------------------------------
				FitnessMin = 999999;
				for (int i=-StepNum; i<=StepNum; i++)
				{
					scores[0]=ScoreBest[0];
					scores[1]=ScoreBest[1];
					scores[2]=i*step[2];
					this->ComputeCoefVertaxBackup(Vert, scores);
					this->GetBackup();
					this->Inv_DWT(0,4);
					this->CC_Subd(4,5);
					//cout<<Vert->X<<" "<<Vert->Y<<" "<<Vert->Z<<": "<<this->ComputeFitness()<<endl;
					Fitness = this->ComputeFitness();
					//myfile<<Fitness<<" ";
					if (Fitness<FitnessMin)
					{
						FitnessMin=Fitness;
						FitnessIndex=i;
						ScoreBest[2]=scores[2];
					}					
				}
				//Set to best fit
				scores[0]=ScoreBest[0];
				scores[1]=ScoreBest[1];
				scores[2]=ScoreBest[2];
				Vert->OpiScores[0]=ScoreBest[0];
				Vert->OpiScores[1]=ScoreBest[1];
				Vert->OpiScores[2]=ScoreBest[2];
				this->ComputeCoefVertaxBackup(Vert, scores);
				//myfile<<"best fit: "<<FitnessMin<<" best index: "<<FitnessIndex<<endl;

				//myfile<<"best scores: "<<ScoreBest[0]<<" "<<ScoreBest[1]<<" "<<ScoreBest[2]<<" "<<endl;

				//just for output the interim results----------------start---------------
				/*this->GetBackup(); //test
				this->Inv_DWT(0,4);//test
				this->CC_Subd(4,5);//test
				sprintf(IntrimResultName,"Interim%d.off",IntrimResultCounter);
				IntrimResultCounter++;
				this->MeshOut(IntrimResultName,5);*/
				//just for output the interim results----------------end------------------
			}
		}
		this->SetNewMean();		
	}

	void Surface::SetNewMean(void)
	{
		int NumLev = this->NumLevel;
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;

		//traversal the vertex
		for (Level=0;Level<NumLev;Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;

				//backup the value
				Vert->PriorMEAN[0]=Vert->X;
				Vert->PriorMEAN[1]=Vert->Y;
				Vert->PriorMEAN[2]=Vert->Z;
			}		
		}
	}

	float Surface::QuadArea(Mynode* TheNode)
	{
		float area=0;
		float a1,b1,c1,s1,a2,b2,c2,s2;

		a1=sqrt( pow(TheNode->TopLeft->x-TheNode->BotLeft->x,2)
			    +pow(TheNode->TopLeft->y-TheNode->BotLeft->y,2)
				+pow(TheNode->TopLeft->z-TheNode->BotLeft->z,2) );
		b1=sqrt( pow(TheNode->BotLeft->x-TheNode->BotRigh->x,2)
				+pow(TheNode->BotLeft->y-TheNode->BotRigh->y,2)
				+pow(TheNode->BotLeft->z-TheNode->BotRigh->z,2) );
		c1=sqrt( pow(TheNode->BotRigh->x-TheNode->TopLeft->x,2)
				+pow(TheNode->BotRigh->y-TheNode->TopLeft->y,2)
				+pow(TheNode->BotRigh->z-TheNode->TopLeft->z,2) );

		a2=sqrt( pow(TheNode->TopLeft->x-TheNode->TopRigh->x,2)
				+pow(TheNode->TopLeft->y-TheNode->TopRigh->y,2)
				+pow(TheNode->TopLeft->z-TheNode->TopRigh->z,2) );
		b2=sqrt( pow(TheNode->TopLeft->x-TheNode->BotRigh->x,2)
				+pow(TheNode->TopLeft->y-TheNode->BotRigh->y,2)
				+pow(TheNode->TopLeft->z-TheNode->BotRigh->z,2) );
		c2=sqrt( pow(TheNode->BotRigh->x-TheNode->TopRigh->x,2)
				+pow(TheNode->BotRigh->y-TheNode->TopRigh->y,2)
				+pow(TheNode->BotRigh->z-TheNode->TopRigh->z,2) );

		s1=(a1+b1+c1)/2;
		s2=(a2+b2+c2)/2;

		area=sqrt(s1*(s1-a1)*(s1-b1)*(s1-c1))+sqrt(s2*(s2-a2)*(s2-b2)*(s2-c2));

		return area;
	}

	float Surface::PriTerm(void)
	{
		int Level;
		WaveCoeffiList::iterator pos;
		Point3D *Vert;

		float PriorTerm=0.0;
		
		for (Level=0; Level<3; Level++)
		{
			for (pos=Coefficients[Level].begin(); pos!=Coefficients[Level].end(); ++pos) 
			{
				Vert=(*pos).Coeffi;
				/*PriorTerm=PriorTerm+pow((Vert->OpiScores[0]-Vert->PriorMEAN[0]),2)/Vert->PriorSTD[0];
				PriorTerm=PriorTerm+pow((Vert->OpiScores[1]-Vert->PriorMEAN[1]),2)/Vert->PriorSTD[1];
				PriorTerm=PriorTerm+pow((Vert->OpiScores[2]-Vert->PriorMEAN[2]),2)/Vert->PriorSTD[2];*/
				PriorTerm=PriorTerm+pow(Vert->OpiScores[0],2)/Vert->PriorSTD[0];
				PriorTerm=PriorTerm+pow(Vert->OpiScores[1],2)/Vert->PriorSTD[1];
				PriorTerm=PriorTerm+pow(Vert->OpiScores[2],2)/Vert->PriorSTD[2];
			}
		}
		//cout<<"PriorTerm = "<<PriorTerm*Alpha<<endl;
		return PriorTerm;		
	}