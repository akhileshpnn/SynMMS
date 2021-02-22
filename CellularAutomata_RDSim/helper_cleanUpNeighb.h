template <typename Comp_Type> void cellField<Comp_Type>::cleanUpNeighb(int Depth, int tSize, long int *tNeighb, long int &Pos) {
	if (Depth<SubSize && tSize>=SubCellsVol[Depth]) {

		long int *subNeighb=new long int[tSize];
		int *SSs=SubSizes[Depth-1],*Ss=SubSizes[Depth],nSize=0;

	    int dim_ext[DIM+2];
		dim_ext[0]=(SSs[0]>0)?SSs[0]:1;
		dim_ext[1]=(SSs[1]>0)?SSs[1]:1;
		dim_ext[2]=dim_ext[0]*dim_ext[1];
		dim_ext[3]=(Ss[0]>0)?Ss[0]:1;
		dim_ext[4]=((Ss[1]>0)?Ss[1]:1)*dim_ext[3];

		int x,y,z;
		long int val,subSubOff,subOff;
		int species=Pos/CellsSize;

		subOff=(SubCellsOffset[Depth]-SubCellsOffset[Depth-1])*species/numSpecies+SubCellsOffset[Depth-1];
		subSubOff=(SubCellsOffset[Depth+1]-SubCellsOffset[Depth])*species/numSpecies+SubCellsOffset[Depth];
		
		for (int i=0; i<tSize; i++) {
			val=tNeighb[i]-subOff;
			x=(val%dim_ext[0])/2; y=((val/dim_ext[0])%dim_ext[1])/2; z=(val/dim_ext[2])/2;
			if ((x==Ss[0] && (Ss[0])) || (y==Ss[1] && (Ss[1])) || (z==Ss[2] && (Ss[2]))) subNeighb[i]=-1;
			else subNeighb[i]=(x+y*dim_ext[3]+z*dim_ext[4])+subSubOff;
		}
		
		for (int i=0; i<tSize; i++) {
			if (subNeighb[i]>=0) {
				int num=1;
				for (int j=i+1; (j<tSize && num<SubCellsVol[Depth]); j++) {
					if (subNeighb[i]==subNeighb[j]) {
						num++;
						subNeighb[j]=-2;
					}
				}
				if (num==SubCellsVol[Depth]) {
					tNeighb[i]=subNeighb[i];
					subNeighb[i]=-4;
					NNeighbC[Pos]-=(SubCellsVol[Depth]-1);
					nSize++;
					for (int j=i+1; j<tSize; j++) if (subNeighb[j]==-2) subNeighb[j]=-3;
				}
				else {
					for (int j=i+1; j<tSize; j++) if (subNeighb[j]==-2) subNeighb[j]=-1;
					subNeighb[i]=-1;
				}
			}
		}

		long int *new_tNeighb=new long int[nSize];
		nSize=0;
		for (int i=0; i<tSize; i++) {
			if (subNeighb[i]==-4) {
				new_tNeighb[nSize++]=tNeighb[i];
				tNeighb[i]=-1;
			}
			else if (subNeighb[i]==-3) {
				tNeighb[i]=-1;
			}
		}
		cleanUpNeighb(Depth+1, nSize, new_tNeighb, Pos);
        delete []new_tNeighb;
        delete []subNeighb;
	}
	else {
		Neighb[Pos]= new long int[NNeighbC[Pos]];
		NNeighbC[Pos]=0;
	}
	
	for (int i=0; i<tSize; i++)
		if (tNeighb[i]>=0) {
			Neighb[Pos][NNeighbC[Pos]++]=tNeighb[i];
		}
};
