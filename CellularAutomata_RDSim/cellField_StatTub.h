#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include "vec3d.h"
#include <list>
#include <cstring>

#ifndef CELLFIELD
#define CELLFIELD

//#define _TIRF_IMAGE
//#define _LOAD_3D_PNGSTACK


using namespace std;

const int DIM=3;

const int numSpecies=9;

extern unsigned char *imageData;  /* image char input matrix */
extern double *initVals,*normVals;
extern mxArray *Params;

template <typename Comp_Type> class cellField {
	public:
	  cellField(int* _size);
      ~cellField();
	  void showField(ostream &out=cout);
	  Comp_Type* showField(Comp_Type *Out);
	  void diffuse();
      void diffuse_varGeom();
	  void interact();
      void deallocClass();
      
	private:
	  int size[DIM];
	  int *anIso_pos[DIM];
	  int *anIso_neg[DIM];
	  int *maxRecip;
	  int numStates, numComparts;

      long int PMArea,cytArea,ERArea,REArea,intArea;
	  long int CellsSize,SpeciesSize,SubCellsSize;
      int SubSize,**SubSizes,*SubCellsOffset,*SubCellsVol,*SubCellsConstits;
	  Comp_Type *Cells;
//      Comp_Type Stat,Stat_p,Stat_p_Tub,Tub,Stat_Tub;
	  int *State;
	  vec3d *Center;
	  int *NRecip;
	  int *NNeighb;
	  int *NNeighbC;
	  long int **Neighb;
	  Comp_Type maxVal;
	  Comp_Type *SubCells;

	  void initSubGrid();
	  void initNeighb(double *DiffRad);
	  void subSum();
	  void subDiv();
	  void cleanUpNeighb(int Depth, int tSize, long int *tNeighb, long int &Pos);
	  double anIsoCheck(vec3d direction, int state);
      void clearCenter();
      void findPM();
      void oneWayConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere);
      void exchangeConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere);
      void dissociateConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere, Comp_Type &andWhere);
      void associateConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andWhere, Comp_Type &toWhere);
      void asDissConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andWhere, Comp_Type &toWhere);
      void asymAsDissConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andFromWhere, Comp_Type &toWhere, Comp_Type &andToWhere);
};

#include "constructor_StatTub.h"
#include "helper_deallocClass.h"
#include "helper_clearCenter.h"
#include "helper_findPM.h"
#include "helper_initSubGrid.h"    
#include "helper_cleanUpNeighb.h"
#include "helper_anIsoCheck.h"
#include "helper_initNeighb_StatTub.h"
#include "helper_subSum.h"
#include "helper_subDiv.h"
#include "helper_diffuse.h"
#include "helper_diffuse_varGeom.h"
#include "helper_interact_StatTub.h"


template <typename Comp_Type> cellField<Comp_Type>::~cellField(){
    deallocClass();
}

	

template <typename Comp_Type> Comp_Type* cellField<Comp_Type>::showField(Comp_Type *Out) {
#ifdef _MATLAB_MEX
	std::memcpy (Out, Cells, SpeciesSize * sizeof(Comp_Type));

    return Out;
#else
// Strange output-scheme; disregard!
    for (int i=0; i<size[0]*size[1]; i++) {
		Out[i]=Cells[i];
		for (int sp=1; sp<numSpecies; sp++) Out[i]+=Cells[i+CellsSize];
	}
	return &Out[size[0]*size[1]];
#endif
}

template <typename Comp_Type> void cellField<Comp_Type>::showField(ostream &out) {
	out.write((char*)Cells, SpeciesSize*sizeof(Comp_Type));
/*
static int ImgCounter=0;
static Image ImgOut(size[0],size[1],1);

// configure maxVal!!!

Pos=0;
	for (int k=0; k<size[2]; k++) {
		for (int j=0; j<size[1]; j++) 
			for (int i=0; i<size[0]; i++) {
				int color(255.*Cells[Pos++]/maxVal);
				ImgOut.ImgData()[i][j][0]=static_cast<char>(color);
			}
stringstream Filename;
Filename << OutNameBase;
Filename << setw(6) << setfill('0') << ImgCounter++;
Filename << ".png";

ImgOut.writeToPNG(Filename.str().c_str(),9);
	}
*/
};
#endif
