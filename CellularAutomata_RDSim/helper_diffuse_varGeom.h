template <typename Comp_Type> void cellField<Comp_Type>::diffuse_varGeom() {
	int specPos=0, mN, state=0;

    for (int sp=0; sp<numSpecies; sp++) {
        for (int pos=0; pos<CellsSize; ++pos) {
            if (NRecip[specPos]) {
                mN=maxRecip[sp];
                SubCells[specPos]=Cells[specPos]/mN;
                Cells[specPos]=(mN-NRecip[specPos])*SubCells[specPos];
            }
            else {
                SubCells[specPos]=Cells[specPos];
                Cells[specPos]=0.;
            }
            specPos++;
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
