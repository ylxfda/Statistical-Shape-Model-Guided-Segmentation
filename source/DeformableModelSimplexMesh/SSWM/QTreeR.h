
#ifndef QTreeR_h  
#define QTreeR_h       1

#include "Global.h"
#include "Mytree.h"

class QTreeR
{
public:
	Mytree* tree;
	
	QTreeR* NeibEast;
	QTreeR* NeibWest;
	QTreeR* NeibNorth;
	QTreeR* NeibSouth;

	char ID_NeibEast;  //its location seen by its neibor
	char ID_NeibWest;
	char ID_NeibNorth;
	char ID_NeibSouth;

	QTreeR();

	~QTreeR();

protected:
private:
};


#endif