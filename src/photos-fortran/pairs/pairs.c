// after investigations PHORO3 of PhotosUtilities.cxx will be used instead
// but it must be checked first if it works

void ROTOD3_(double ANGLE,double PVEC[4],double QVEC[4]){
  int j=1;  // convention of indices of Riemann space must be preserved.
  double CS,SN;
  CS=cos(ANGLE)*PVEC[1-j]-sin(ANGLE)*PVEC[2-j];
  SN=sin(ANGLE)*PVEC[1-j]+cos(ANGLE)*PVEC[2-j];
  QVEC[1-j]=CS;
  QVEC[2-j]=SN;
  QVEC[3-j]=PVEC[3-j];
  QVEC[4-j]=PVEC[4-j];
}
