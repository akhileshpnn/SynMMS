template <> void cellField<long int>::diffuse() {
	static int *Rest;
	static bool init=true;
	static int debug_cntr=0;
	static int debug_max=1;
	static ofstream debug;
	
	if (init) {Rest=new int[CellsSize];}
	if (debug_cntr==0) debug.open("debug.log");
	
	for (long int i=0; i<CellsSize; ++i) {
		if (NNeighb[i]) SubCells[i]=Cells[i]/NNeighb[i];
		else SubCells[i]=0;
		Rest[i]=Cells[i]-SubCells[i]*NNeighb[i];
		Cells[i]=0;
	}
	for (long int i=CellsSize; i<SubCellsSize; ++i) SubCells[i]=0;

	subSum();

	for (int i=0; i<CellsSize; ++i) {
		long int *Index = Neighb[i];
		for (int j=0; j<NNeighbC[i]; ++j) {
			Cells[i]+=SubCells[Index[j]];
		}
	}
	
	for (int i=0; i<SubCellsSize; ++i) SubCells[i]=0;
	for (int i=0; i<CellsSize; ++i) {
		long int *Index = Neighb[i];
		int Depth,k,l;
		for (int n=0; (n<Rest[i] && NNeighbC[i]); n++) {
			k=(int)((NNeighb[i]*rand_norm)*rand());
			if (k>=NNeighb[i]) k=NNeighb[i]-1;
			if (SubSize>1) {
				Depth=0;l=NNeighbC[i]-1;
				while (k>=SubCellsConstits[Depth] && l>0) {
					k-=SubCellsConstits[Depth];
					l--;
					if (Depth<SubSize-1) if (Index[l]>SubCellsOffset[Depth+1]) Depth++;
				}
			}
			else l=NNeighbC[i]-1-k;
			SubCells[Index[l]]++;
//			Cells[Index[k]]++;
		}
	}
if (init) {init=false;}

//	maxVal=0;
	for (int i=0; i<CellsSize; ++i) {
		Cells[i]+=SubCells[i];
//		maxVal=MAX(maxVal,Cells[i]);
	}

};
