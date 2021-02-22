template <typename Comp_Type> double cellField<Comp_Type>::anIsoCheck(vec3d direction, int state) {
		double direction_aniso[3]={0.,0.,0.};

		if (direction.X()>0) {
			if (anIso_pos[0][state]>0) direction_aniso[0]=direction.X()/anIso_pos[0][state];
			else return -1;
		}
		else if (direction.X()<0) {
			if (anIso_neg[0][state]>0) direction_aniso[0]=direction.X()/anIso_neg[0][state];
			else return -1;
		}
		if (direction.Y()>0) {
			if (anIso_pos[1][state]>0) direction_aniso[1]=direction.Y()/anIso_pos[1][state];
			else return -1;
		}
		else if (direction.Y()<0) {
			if (anIso_neg[1][state]>0) direction_aniso[1]=direction.Y()/anIso_neg[1][state];
			else return -1;
		}
		if (direction.Z()>0) {
			if (anIso_pos[2][state]>0) direction_aniso[2]=direction.Z()/anIso_pos[2][state];
			else return -1;
		}
		else if (direction.Z()<0) {
			if (anIso_neg[2][state]>0) direction_aniso[2]=direction.Z()/anIso_neg[2][state];
			else return -1;
		}
		return sqrt(SQUARE(direction_aniso[0])+SQUARE(direction_aniso[1])+SQUARE(direction_aniso[2]));
};
