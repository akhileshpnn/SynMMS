template <typename Comp_Type> void cellField<Comp_Type>::deallocClass() {
    delete []Cells;
    delete []State;
    delete []Center;

    delete []NRecip;
    delete []NNeighb;
    for (int i=0;i<SpeciesSize;i++) if (NNeighbC[i]>0) delete[] Neighb[i];
    delete []Neighb;
    delete []NNeighbC;

    for (int i=0; i<DIM; i++) {
	  delete []anIso_pos[i];
	  delete []anIso_neg[i];
    }

    delete []maxRecip;
    for (int i=0;i<SubSize;++i) delete[] SubSizes[i];
    delete []SubSizes;
    delete []SubCellsOffset;
    delete []SubCellsVol;
    delete []SubCellsConstits;
    delete []SubCells;

#ifndef _MATLAB_MEX
//    delete []DiffRad;
#endif                

}
