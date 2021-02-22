template <typename Comp_Type> void cellField<Comp_Type>::subSum() {
	int yBase, zBase, zOffSet, offSet;
	int Pos, Sub_Pos;
	int *SSs;
	int dim_ext[DIM];
	bool OddCells,OddRows;

	int spPosShift, spSubPosShift;
	spSubPosShift=(SubCellsOffset[1]-SubCellsOffset[0])/numSpecies;
	for (int i=1; i<SubSize; ++i) {
		SSs=SubSizes[i];
		for (int j=0; j<DIM; j++) dim_ext[j]=(SSs[j]>0)?(2*SSs[j]):1;
		OddCells=((SubSizes[i-1][0]&1) && (dim_ext[0]>1));
		OddRows =((SubSizes[i-1][1]&1) && (dim_ext[1]>1));
		yBase=SSs[0];
		zBase=SSs[1];
		if (yBase) {
			if (zBase) zBase*=yBase;
			else zBase=yBase;
		}
		else if (zBase) yBase=1;
		else zBase=1;
		offSet=SubCellsOffset[i]; zOffSet=SubCellsOffset[i];
		Pos=SubCellsOffset[i-1];
		spPosShift=spSubPosShift;
		spSubPosShift=(SubCellsOffset[i+1]-SubCellsOffset[i])/numSpecies;

		for (int z=0;z<dim_ext[2];z++) {
			for (int y=0;y<dim_ext[1];y++) {
				Sub_Pos=offSet;
				for (int x=0;x<dim_ext[0];x++) {
					SubCells[Sub_Pos]+=SubCells[Pos];
					for (int sp=1; sp<numSpecies; sp++) SubCells[Sub_Pos+sp*spSubPosShift]+=SubCells[Pos+sp*spPosShift];
					Pos++;
					if (x&1) Sub_Pos++;
				}
				if (OddCells) Pos++; 
				if (y&1) offSet+=yBase;
			}
			if (OddRows)  Pos+=SubSizes[i-1][0]; 
			if (z&1) zOffSet+=zBase;
			offSet=zOffSet;
		}
	}
};
