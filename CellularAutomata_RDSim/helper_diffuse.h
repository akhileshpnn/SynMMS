template <typename Comp_Type> void cellField<Comp_Type>::diffuse() {
	int mN, state=0;
	for (int i=0; i<SpeciesSize; ++i) {
		if (NRecip[i]) {
            if (State[i%CellsSize]==0) mexPrintf("SHIT");
//			state=(State[i%CellsSize]-1)+i/CellsSize*numComparts;
//			mN=maxRecip[state];
			SubCells[i]=Cells[i]/NRecip[i];
			Cells[i]=0.;
		}
		else {
			SubCells[i]=Cells[i];
            Cells[i]=0.;
        }
    }
	for (int i=SpeciesSize; i<SubCellsSize; ++i) SubCells[i]=0;
	subSum();

	long int *Index;
	for (int i=0; i<SpeciesSize; ++i) {
		Index = Neighb[i];
		for (int j=0; j<NNeighbC[i]; ++j) {
			Cells[i]+=SubCells[Index[j]];
		}
	}
};
