template <typename Comp_Type> void cellField<Comp_Type>::subDiv() {
	int yBase, zBase;
	long int Pos;
    int *SSs;
    int Vol;
	int dim_ext[DIM];
	bool OddCells,OddRows;

	for (int i=SubSize-1; i>0; --i) {
		SSs=SubSizes[i-1];
		for (int j=0; j<DIM; j++) {
			dim_ext[j]=(SSs[j]>1)?2:1;
		}
		OddCells=((SSs[0]&1) && (SSs[0]>1));
		OddRows=((SSs[1]&1) && (SSs[1]>1));
		yBase=SSs[0];
		zBase=SSs[1];
		if (yBase) zBase*=yBase;
		else if (zBase) yBase=1;
		else zBase=1;

		Pos=SubCellsOffset[i];
		Vol=SubCellsVol[i];
		long int *SubPos = new long int[Vol];

		Vol=0;
		for (int z=0;z<dim_ext[2];z++) {
			for (int y=0;y<dim_ext[1];y++) {
				for (int x=0;x<dim_ext[0];x++) {
					SubPos[Vol++]=SubCellsOffset[i-1]+x+yBase*y+zBase*z;
				}
			}
		}

		SSs=SubSizes[i];
		for (int j=0; j<DIM; j++) {
			dim_ext[j]=(SSs[j]>0)?(SSs[j]):1;
		}
		int add;
		for (int z=0;z<dim_ext[2];z++) {
			for (int y=0;y<dim_ext[1];y++) {
				for (int x=0;x<dim_ext[0];x++) {
					if (SubCells[Pos]) {
						int S=SubCells[Pos]/Vol;
						if (S) {
							for (int j=0; j<Vol; j++) SubCells[SubPos[j]]+=S;
							S=SubCells[Pos]-S*Vol;
						}
						else S=SubCells[Pos];
						if (S)
							for (int j=0; j<S; j++) {
								int k=(int)(Vol*(rand()*rand_norm));
								if (k==Vol) k--;
								SubCells[SubPos[k]]++;
							}
					}
					if (dim_ext[0]>1) for (int j=0; j<Vol; j++) SubPos[j]+=2;
					Pos++;
				}
				if (OddCells) add=1; else add=0; 
				if (dim_ext[1]>1) add+=yBase;
				for (int j=0; j<Vol; j++) SubPos[j]+=add;
			}
			if (OddRows) add=yBase; else add=0; 
			if (dim_ext[2]>1) add+=zBase;
			for (int j=0; j<Vol; j++) SubPos[j]+=add;
		}
        delete []SubPos;
	}
};
