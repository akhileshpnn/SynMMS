template <typename Comp_Type> void cellField<Comp_Type>::clearCenter() {
    long int OffSet=size[0]*size[1];
    long int Pos=0,NPos;
    for (int z=0;z<size[2];z++) {
        for (int y=0;y<size[1];y++) {
            for (int x=0;x<size[0];x++) {
                if (State[Pos]==127) {
                    bool innen=true;
                    for (int d=-1; d<=1 && innen; d+=2) {
                        if (size[0]>1) {
                            if (x+d<0 || x+d>=size[0] ) innen=false;
                            else NPos=Pos+d;
                            if (innen) if (State[NPos]<127 || State[NPos]>128) innen=false;
                        }
                        if (size[1]>1) {
                            if (y+d<0 || y+d>=size[1] ) innen=false;
                            else NPos=Pos+d*size[0];
                            if (innen) if (State[NPos]<127 || State[NPos]>128) innen=false;
                        }
                        if (size[2]>1) {
                            if (z+d<0 || z+d>=size[2] ) innen=false;
                            else NPos=Pos+d*OffSet;
                            if (innen) if (State[NPos]<127 || State[NPos]>128) innen=false;
                        }
                    }
                    if (innen) State[Pos]=128;
                }
                Pos++;
            }
        }
    }
    for (Pos=0;Pos<CellsSize;Pos++) {
        if (State[Pos]==128) State[Pos]=0;
        Cells[Pos]=State[Pos];
    }
}
