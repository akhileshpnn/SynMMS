template <typename Comp_Type> void cellField<Comp_Type>::initSubGrid() {
	SubSizes=new int*[SubSize];
	for (int i=0;i<SubSize;++i) {
		SubSizes[i]=new int[DIM];
		for (int j=0; j<DIM; j++) SubSizes[i][j]=0;
	}
	SubCellsOffset=new int[SubSize+2];
	SubCellsOffset[0]=0;
	SubCellsVol=new int[SubSize+1];
	SubCellsConstits=new int[SubSize+1];
	SubCellsVol[0]=1;
	SubCellsSize=size[0]*size[1]*size[2]*numSpecies;
	SubCellsOffset[1]=SubCellsSize;
	for (int j=0; j<DIM; j++) SubSizes[0][j]=size[j];
	for (int i=1;i<SubSize;++i) {
		SubCellsVol[i]=1;
		for (int j=0; j<DIM; ++j) {
			SubSizes[i][j]=(SubSizes[i-1][j]/2);
			if (SubSizes[i][j]) SubCellsVol[i]*=2;
		}
		if (SubCellsVol[i]>1) {
			SubCellsOffset[i+1]=1;
			for (int j=0; j<DIM; ++j) {
				if (SubSizes[i][j]) SubCellsOffset[i+1]*=SubSizes[i][j];
			}
			SubCellsOffset[i+1]*=numSpecies;
			SubCellsSize+=SubCellsOffset[i+1];
			SubCellsOffset[i+1]+=SubCellsOffset[i];			
		}
		else {
			SubSize=i;
		}
	}
	SubCellsOffset[SubSize+1]=-1;

	SubCells = new Comp_Type[SubCellsSize];
	SubCellsConstits[0]=1;
	for (int i=1; i<SubSize; i++) {
		SubCellsConstits[i]=SubCellsConstits[i-1]*SubCellsVol[i];
	}
};
