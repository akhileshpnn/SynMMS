template <typename Comp_Type> void cellField<Comp_Type>::interact() {
      Comp_Type Stat,Stat_p,Tub,Stat_p_Tub01,Stat_Tub01,Stat_p_Tub10,Stat_Tub10,Stat_p_Tub11,Stat_Tub11;

    double* p=static_cast<double*>(mxGetData(mxGetField(Params, 0, "IntPar")));
    
//Species:  0--cytosolic; 1--cytosolic; 2--membrane; 3--membrane; 4--cytosolic; 5--cytosolic;
//	        6--membrane;  7--membrane;  8--membrane;
//Locations:  0--outside; 1--PM; 63--recycling endosome; 95--perinuclear membranes; 127--endomembranes; 255--cytosol
//Params:   0--PPTase, 1--low kinase, 2--high kinase, 
//          3--diss Stat_Tub01_p, 4--diss Stat_Tub01, 5--ass Stat01_p, 6--ass Stat01
//          7--diss Stat_Tub10_p, 8--diss Stat_Tub10, 9--ass Stat10_p,10--ass Stat10
    
    for (int i=0; i<CellsSize; ++i) {
       Comp_Type &Tub_old           =Cells[i];                    Tub=           Tub_old;
       Comp_Type &Stat_old          =Cells[i+  CellsSize];       Stat=          Stat_old;
       Comp_Type &Stat_p_old        =Cells[i+2*CellsSize];     Stat_p=        Stat_p_old;
       Comp_Type &Stat_Tub01_old    =Cells[i+3*CellsSize];   Stat_Tub01=  Stat_Tub01_old;       
       Comp_Type &Stat_p_Tub01_old  =Cells[i+4*CellsSize]; Stat_p_Tub01=Stat_p_Tub01_old;
       Comp_Type &Stat_Tub10_old    =Cells[i+5*CellsSize];   Stat_Tub10=  Stat_Tub10_old;       
       Comp_Type &Stat_p_Tub10_old  =Cells[i+6*CellsSize]; Stat_p_Tub10=Stat_p_Tub10_old;
       Comp_Type &Stat_Tub11_old    =Cells[i+7*CellsSize];   Stat_Tub11=  Stat_Tub11_old;       
       Comp_Type &Stat_p_Tub11_old  =Cells[i+8*CellsSize]; Stat_p_Tub11=Stat_p_Tub11_old;

             
       //everywhere
       //Stathmin(p) binds free Tubulin and Stathmin-Tubulin(p) dissociates -- 3 species!
       asDissConc((p[ 6]*Stat*Tub-p[4]*Stat_Tub01),Stat_old,Tub_old,Stat_Tub01_old);
       asDissConc((p[10]*Stat*Tub-p[8]*Stat_Tub10),Stat_old,Tub_old,Stat_Tub10_old);

       asDissConc((p[10]*Stat_Tub01*Tub-p[8]*Stat_Tub11),Stat_Tub01_old,Tub_old,Stat_Tub11_old);
       asDissConc((p[ 6]*Stat_Tub10*Tub-p[4]*Stat_Tub11),Stat_Tub10_old,Tub_old,Stat_Tub11_old);
       
       asDissConc((p[5]*Stat_p*Tub-p[3]*Stat_p_Tub01),Stat_p_old,Tub_old,Stat_p_Tub01_old);
       asDissConc((p[9]*Stat_p*Tub-p[7]*Stat_p_Tub10),Stat_p_old,Tub_old,Stat_p_Tub10_old);

       asDissConc((p[9]*Stat_p_Tub01*Tub-p[7]*Stat_p_Tub11),Stat_p_Tub01_old,Tub_old,Stat_p_Tub11_old);
       asDissConc((p[5]*Stat_p_Tub10*Tub-p[3]*Stat_p_Tub11),Stat_p_Tub10_old,Tub_old,Stat_p_Tub11_old);
       
       //only in cyt:
       if (State[i]==255) {
       //Stathmin_p gets (de-)phosphorylated
         exchangeConc(p[0]*Stat_p      -p[1]*Stat,      Stat_p_old,      Stat_old);
         exchangeConc(p[0]*Stat_p_Tub01-p[1]*Stat_Tub01,Stat_p_Tub01_old,Stat_Tub01_old);
         exchangeConc(p[0]*Stat_p_Tub10-p[1]*Stat_Tub10,Stat_p_Tub10_old,Stat_Tub10_old);
         exchangeConc(p[0]*Stat_p_Tub11-p[1]*Stat_Tub11,Stat_p_Tub11_old,Stat_Tub11_old);
       }

       //in all membranes
       if (State[i]>0 && State[i]<255) {
       }
       //only at the PM
       if (State[i]==1) {
       //Stathmin_p gets (de-)phosphorylated
           exchangeConc(p[0]*Stat_p      -(p[2]+p[1])*Stat,      Stat_p_old,      Stat_old);
           exchangeConc(p[0]*Stat_p_Tub01-(p[2]+p[1])*Stat_Tub01,Stat_p_Tub01_old,Stat_Tub01_old);
           exchangeConc(p[0]*Stat_p_Tub10-(p[2]+p[1])*Stat_Tub10,Stat_p_Tub10_old,Stat_Tub10_old);
           exchangeConc(p[0]*Stat_p_Tub11-(p[2]+p[1])*Stat_Tub11,Stat_p_Tub11_old,Stat_Tub11_old);
       //Stathmin_p only gets phosphorylated
       //exchangeConc(-p[2]*Stat,Stat_p_old,Stat_old);
       //exchangeConc(-p[2]*Stat_Tub,Stat_p_Tub_old,Stat_Tub_old);
       }
       //only in ER:
       if (State[i]==127) {
       }
       //only in PMC:
       if (State[i]==95) {
       }
       //only in lysosomes:
       if (State[i]==63) {
       }
 
     }
}

template <typename Comp_Type> void cellField<Comp_Type>::exchangeConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere) {
  if (fromWhere<howMuch || toWhere<-howMuch) {
      if (howMuch>0) howMuch=fromWhere;
      else howMuch=-toWhere;
  }
  fromWhere-=howMuch;toWhere+=howMuch;
}

template <typename Comp_Type> void cellField<Comp_Type>::oneWayConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere) {
  if (howMuch>fromWhere) howMuch=fromWhere;
  toWhere+=howMuch;
  fromWhere-=howMuch; 
}

template <typename Comp_Type> void cellField<Comp_Type>::dissociateConc(double howMuch, Comp_Type &fromWhere, Comp_Type &toWhere, Comp_Type &andWhere) {
  if (howMuch>fromWhere) howMuch=fromWhere;
  toWhere+=howMuch;
  andWhere+=howMuch;
  fromWhere-=howMuch; 
}

template <typename Comp_Type> void cellField<Comp_Type>::associateConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andWhere, Comp_Type &toWhere) {
  if (howMuch>fromWhere) howMuch=fromWhere;
  if (howMuch>andWhere) howMuch=andWhere;
  toWhere+=howMuch;
  fromWhere-=howMuch; 
  andWhere-=howMuch;
}

template <typename Comp_Type> void cellField<Comp_Type>::asDissConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andWhere, Comp_Type &toWhere) {
// A+B <-> AB
  if (fromWhere<howMuch || andWhere<howMuch || toWhere<-howMuch) {
      if (howMuch>0) howMuch=MIN(fromWhere,andWhere);
      else howMuch=-toWhere;
  }
  toWhere+=howMuch;
  andWhere-=howMuch;
  fromWhere-=howMuch; 
}


template <typename Comp_Type> void cellField<Comp_Type>::asymAsDissConc(double howMuch, Comp_Type &fromWhere, Comp_Type &andFromWhere, Comp_Type &toWhere, Comp_Type &andToWhere) {
// A+B -> AB -> A+B*
  if (fromWhere<howMuch || andFromWhere<howMuch || toWhere<-howMuch) {
      if (howMuch>0) howMuch=MIN(fromWhere,andFromWhere);
      else howMuch=-toWhere;
  }
  toWhere+=howMuch;
  if (howMuch<0) andToWhere-=howMuch;
  if (howMuch>0) andFromWhere-=howMuch;
  fromWhere-=howMuch; 
}

