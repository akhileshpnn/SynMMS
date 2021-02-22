template <typename Comp_Type> void cellField<Comp_Type>::initNeighb(double *DiffRad) {
	NRecip   = new int[SpeciesSize];
	NNeighb  = new int[SpeciesSize];
	NNeighbC = new int[SpeciesSize];
	Neighb   = new long int*[SpeciesSize];

	int s, *rrange=new int[numSpecies];
	for (s=0; s<numSpecies; s++) {
		rrange[s]=(int)(DiffRad[s])+1;
	}
	int Volume=1;
	for (int i=0; i<DIM; i++) {
		int maxRange=rrange[0]*(anIso_pos[i][0]+anIso_neg[i][0]);
        for (s=1; s<numSpecies; s++) {
			maxRange=MAX(rrange[s]*(anIso_pos[i][0]+anIso_neg[i][0]),maxRange);
		}
	    if (maxRange==0) maxRange=1;
        else maxRange+=2;
        Volume*=maxRange;
	}

	int NN;

    long int *tempNeighb=new long int[Volume];

	for (s=0; s<SpeciesSize; s++) NRecip[s]=0;

	int state;
    long int Pos, SpecPos, sp1Pos;

    int speciesLoc[numSpecies]={255,255,255,255,255,255,255,255,255};
// Comp: cyt 255, ER 127, PMC 95, RE 63, PM 1
// States: cyt 1, ER 2, PM 3
//    int speciesLoc[2]={255,127};
    
    for (int sp=0; sp<numSpecies; sp++) {
        Pos=0;
		for (int z=0;z<size[2];z++) {
			for (int y=0;y<size[1];y++) {
				for (int x=0;x<size[0];x++) {
                  SpecPos=Pos+sp*CellsSize;


//cross-link RE and ER
//                    if (State[Pos] && (State[Pos]==speciesLoc[sp] || speciesLoc[sp]==255 || ((State[Pos]==127 || State[Pos]==63) && speciesLoc[sp]==127))) 
                    if (State[Pos] && (State[Pos]==speciesLoc[sp] || speciesLoc[sp]==255 || (State[Pos]==1 && speciesLoc[sp]==127) || (State[Pos]==63 && speciesLoc[sp]==127) || (State[Pos]==95 && speciesLoc[sp]==127))) {
                        double DR=DiffRad[sp];
    // hack: Randomize Diffusion Radii (with radial gradient) to generate inhomogeneities;
//                        if (sp==1) DR*=rand()*rand_norm*( exp(-sqrt((double)(SQUARE(x-160)+SQUARE(y-160)))/(80)) );
//                          if (sp==1 && State[Pos]!=1) DR=(.75+.5*rand()*rand_norm)*DR;

                        
                          if (DR>0) {
                            NN=0;
                            for (int dx=-anIso_neg[0][0]*rrange[sp];dx<=anIso_pos[0][0]*rrange[sp];dx++) 
                                for (int dy=-anIso_neg[1][0]*rrange[sp];dy<=anIso_pos[1][0]*rrange[sp];dy++) 
                                    for (int dz=-anIso_neg[2][0]*rrange[sp];dz<=anIso_pos[2][0]*rrange[sp];dz++) {
    //noflux randbedingungen wenn dann hier!
                                        bool innen=(x+dx>=0 && x+dx<size[0] && y+dy>=0 && y+dy<size[1] && z+dz>=0 && z+dz<size[2]);
    									if (innen) {    
    //                                    if (true) 
                                            int Index=(size[0]+x+dx)%size[0]+((size[1]+y+dy)%size[1])*size[0]+((size[2]+z+dz)%size[2])*size[0]*size[1];
    //transport across gradient of compartment-numbert from higher to lower                                        
    //                                        if (State[Index] && State[Index]<=speciesLoc[sp] && (State[Index]>=State[Pos] || speciesLoc[sp]==255)) 
    //                                        if (State[Index] && State[Index]<=speciesLoc[sp] && (State[Index]==State[Pos] || speciesLoc[sp]==255 || (State[Index]+State[Pos]==190 && speciesLoc[sp]==127))) 
                                            bool outNuc=State[Index]!=253 && State[Pos]!=253;
                                            bool nucPore=(State[Index]==254 || State[Pos]==254) && (State[Index]==253 || State[Pos]==253);
                                            bool lPDEnoPore= !(sp==4 && ((State[Index]!=254 && State[Pos]==254) || (State[Index]==254 && State[Pos]!=254)));
                                            if (State[Index] && (State[Index]==State[Pos] || State[Index]+State[Pos]==222 || (speciesLoc[sp]==255 && outNuc) || nucPore ) && lPDEnoPore) {
                                                vec3d Direction=Center[Index]+vec3d(dx,dy,dz)-Center[Pos];
                                                double Dist=anIsoCheck(Direction,0);
                                                if (Dist>=0) {
    //link voxels within same compartment 										
    //												if (State[Pos]==State[Index]) 
                                                    if (true) {
                                                        if (Dist<=DR) {
                                                            tempNeighb[NN++]=Index+sp*CellsSize;
                                                            NRecip[Index+sp*CellsSize]++;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                            NNeighb[SpecPos]=NN;
                            NNeighbC[SpecPos]=NN;
                            if (SubSize>1) cleanUpNeighb(1,NN,tempNeighb,SpecPos);
                            else {
                                if (NN>0) 
                                    Neighb[SpecPos]=new long int[NN];
                                else if (State[Pos]==127 || State[Pos])
                                    mexPrintf("No Neighbours! %i",Pos);
                                for (int i=0; i<NN; i++)
                                    Neighb[SpecPos][i]=tempNeighb[i];
                            }
                        }
                        else {
                            NNeighb[SpecPos]=1;
                            NNeighbC[SpecPos]=1;
                            Neighb[SpecPos]=new long int[1];
                            Neighb[SpecPos][0]=SpecPos;
                        }
                    }
					else {
						NNeighb[SpecPos]=0;
						NNeighbC[SpecPos]=0;
					}
                  Pos++;
				}
			}
//cout << "\r                              \rCalculating Neighborhood ... " << (int)((100*(z+sp*size[2]))/(size[2]*numSpecies)) << '%' << flush; 
		}
	}

	for (s=0; s<numSpecies; s++) maxRecip[s]=0;

	for (int sp=0; sp<numSpecies; sp++) 
		for (s=0; s<CellsSize; s++) {
			if (State[s]>0) maxRecip[sp]=MAX(maxRecip[sp],NRecip[s+sp*CellsSize]);
		}

    delete []rrange;
    delete []tempNeighb;
//cout << "\r                              \rCalculating Neighborhood ... 100%" << endl; 
//for (s=0; s<numStates; s++) cout << maxRecip[s]  << " ";
//cout << endl;
//Hardcode membrane-equalization
//maxRecip[0]=MAX(maxRecip[0],maxRecip[2]);
//maxRecip[2]=maxRecip[0];
};
