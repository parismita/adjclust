//heap functions in C

#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MINCL1 1
#define MAXCL1 2
#define MINCL2 3
#define MAXCL2 4
#define LAB1 5
#define LAB2 6
#define POSL 7
#define POSR 8
#define SII 9
#define SJJ 10
#define SIJ 11
#define VALID 12

#define DIS(i,j,p) dissim[p*(j-1)+i-1]
#define CHAIN(i,j) chainedL[(12*(j-1)+i-1)]
#define MERGE(i,j,k) merge[((*p-1)*(j-1)+i-1)]

#define MIN(a,b) ((a) < (b) ? a : b)

//@export
SEXP percDown(SEXP Rpositions, SEXP Rdistances, SEXP Rl, SEXP Rpos){
  int mc, right, left;
  int *positions, *l, *pos;
  double *distances, tmp, val;
  Rpositions = PROTECT(coerceVector(Rpositions, INTSXP));
  positions = INTEGER(Rpositions);
  distances = REAL(Rdistances);
  l = INTEGER(Rl);
  pos = INTEGER(Rpos);
  *pos = *pos - 1;
  val = distances[positions[*pos]-1];

  while((2**pos+1) < *l) {
    //i=i+1;
    if ((2**pos+2) == *l){
      left = 2**pos+1;
      if (val > distances[positions[left]-1]) {  // positions in C start from 0 !!!
	// swap positions
	tmp = positions[*pos];
	positions[*pos] = positions[left];
	positions[left] = tmp;
	// update pos
	*pos = left;
      }
      else
	*pos = *l;
    }
    else {
      left = 2**pos+1;
      right = 2**pos+2;
      mc = right;
      if (distances[positions[left]-1] < distances[positions[right]-1]) // positions in C start from 0 !!!
	mc = left;

      if (val > distances[positions[mc]-1]) {  // positions in C start from 0 !!!
	// swap positions
	tmp = positions[*pos];
	positions[*pos] = positions[mc];
	positions[mc] = tmp;
	// update pos
	*pos = mc;
      }
      else
	*pos = *l;
    }
  }
  UNPROTECT (1) ;
  // return R_NilValue;
  return(Rpositions);
}

int* deleteMin_C(int *positions, double *distances, int l){
  int mc, right, left;
  int pos;
  double tmp, val;

  // pre-processing
  positions[0] = positions[l-1];
  l = l-1;
  pos = 1;

  pos = pos - 1;
  val = distances[positions[pos]-1];

  while((2*pos+1) < l) {
    // i=i+1;
    if ((2*pos+2) == l){
      left = 2*pos+1;
      if (val > distances[positions[left]-1]) {  // positions in C start from 0 !!!
	// swap positions
	tmp = positions[pos];
	positions[pos] = positions[left];
	positions[left] = tmp;
	// update pos
	pos = left;
      }
      else
	pos = l;
    }
    else {
      left = 2*pos+1;
      right = 2*pos+2;
      mc = right;
      if (distances[positions[left]-1] < distances[positions[right]-1]) // positions in C start from 0 !!!
	mc = left;

      if (val > distances[positions[mc]-1]) {  // positions in C start from 0 !!!
	// swap positions
	tmp = positions[pos];
	positions[pos] = positions[mc];
	positions[mc] = tmp;
	// update pos
	pos = mc;
      }
      else
	pos = l;
    }
  }
  // printf("%d\n", i);
  return(positions);
}

int* insertHeap_C(int *positions, double *distances, int l, int key){
  int pos;
  int parent;
  double tmp, val;

  // pre-processing
  positions[l] = key;
  pos = l+1;
  //*l = *l+1;

  parent = pos/2;
  val = distances[positions[pos-1]-1];

  while((pos-1) > 0){
    if ( distances[positions[parent-1]-1] > val){
      // i = i+1;
      // swap positions
      tmp = positions[parent-1];
      positions[parent-1] = positions[pos-1];
      positions[pos-1] = tmp;
      pos = pos/2;
      parent = pos/2;
    } else {
      pos=1;
    }
  }
  // printf("%d\n", i);

  return(positions);
}

// nrow(chainedL)==9
// sense==1/-1
// sense==1 --> right
int* neighborCl_C(int sense, int posMin, double *chainedL){
  int *res = malloc(4*sizeof(int));  // vector of size 4
  if ((sense==1) & (CHAIN(POSR, posMin)<0)){
    res[0] = -1;
  } else if ((sense==1) & (CHAIN(POSR, posMin)>0)) {
    res[0] = CHAIN(POSR, posMin);
    res[1] = CHAIN(MINCL2, res[0]);
    res[2] = CHAIN(MAXCL2, res[0]);
    res[3] = CHAIN(LAB2, res[0]);
  } else if ((sense==-1) & (CHAIN(POSL, posMin)<0)){
    res[0] = -1;
  } else {
    res[0] = CHAIN(POSL, posMin);
    res[1] = CHAIN(MINCL1, res[0]);
    res[2] = CHAIN(MAXCL1, res[0]);
    res[3] = CHAIN(LAB1, res[0]);

  }
  return res;
}

int neiNeighborPos_C(int sense, int posMin, double *chainedL){
  int pos;
  if (sense==1){
    if (CHAIN(POSR, posMin)<0){
      pos = -1;
    } else if (CHAIN(POSR, (int)CHAIN(POSR, posMin))<0) {
      pos = -1;
    }  else {
      pos = CHAIN(POSR, (int)CHAIN(POSR, posMin));
    }
  } else {
    if (CHAIN(POSL, posMin)<0){
      pos = -1;
    } else if (CHAIN(POSL, (int)CHAIN(POSL, posMin))<0) {
      pos = -1;
    }  else {
      pos = CHAIN(POSL, (int)CHAIN(POSL, posMin));
    }
  }
  return(pos);
}

double distance_C(int mini, int maxi, int minj, int maxj, double *dissim, int h, int p){
  double res = INT_MAX;
  int i=0;
  int ni, nj, mI, mJ, mIJ;

  ni = maxi - mini + 1;
  nj = maxj - minj + 1;
  mI = ni - 1;
  mJ = nj - 1;
  mIJ = maxj - mini;

  // given with 2 clusters we have to find min dist.
//max clust items = p, max diff in clust index = h
//time complexity worst case can be reduced to O(ph) as the distancematrix is not sorted and sorting would take O(plogh) time complexity
// hence sequencial method preferred

//traverse dissim from (mini,minj - mini) to (maxi,maxj-maxi) or (maxi,h) and find min dissim.
//time comp O(p*h) worst case
  for(int i=mini;i<=maxi;i++)
  {
    for(int j=minj-i;(j<=maxj-i)&&(j<=h);j++)
    {   
    res = MIN(DIS(i,j,p),res);
    }
  }


return res;
}

SEXP cWardHeaps(SEXP Rdissim, SEXP Rh, SEXP Rp, SEXP RchainedL, SEXP Rpositions, SEXP Rdistances, SEXP RlHeap, SEXP Rmerge, SEXP Rgains, SEXP RblMin){

  int *h, *p, *positions, *lHeap, *merge, *neiL, *neiR, *blMin;
  int posMin, neineiL, neineiR, k;
  double *distances, *chainedL, *gains, d1, d2, dLast, newDR, newDL, clMin_11, clMin_12, clMin_21, clMin_22;
  int jj, step, stepInv;
  double *dissim;

  Rpositions = PROTECT(coerceVector(Rpositions, INTSXP));
  Rmerge = PROTECT(coerceVector(Rmerge, INTSXP));  // matrix
  Rgains = PROTECT(coerceVector(Rgains, REALSXP));  // vector
  Rdistances = PROTECT(coerceVector(Rdistances, REALSXP));
  RchainedL = PROTECT(coerceVector(RchainedL, REALSXP));
  dissim = REAL(Rdissim);
  h = INTEGER(Rh);
  p = INTEGER(Rp);
  chainedL = REAL(RchainedL);
  positions = INTEGER(Rpositions);
  distances = REAL(Rdistances);
  lHeap = INTEGER(RlHeap);
  merge = INTEGER(Rmerge);
  gains = REAL(Rgains);
  blMin = INTEGER(RblMin);

  k = *p - *blMin;

  jj = *p;

  for ( step=1; step < (*p-1); step=step+1 ){
    while(CHAIN(VALID, positions[0])==0){
      // printf("%f\n", distances[positions[0]-1]);
      // printf("%d\n", *lHeap);
      positions = deleteMin_C(positions, distances, *lHeap);
      *lHeap = *lHeap -1;
    }
    posMin = positions[0];
    clMin_11 = CHAIN(MINCL1, posMin);
    clMin_12 = CHAIN(MAXCL1, posMin);
    clMin_21 = CHAIN(MINCL2, posMin);
    clMin_22 = CHAIN(MAXCL2, posMin);
    gains[step-1] = distances[positions[0]-1];

    //remove head
    positions = deleteMin_C(positions, distances, *lHeap);
    *lHeap = *lHeap - 1;

    neiL = neighborCl_C(-1, posMin, chainedL);
    neiR = neighborCl_C(1, posMin, chainedL);
    neineiR = neiNeighborPos_C(1, posMin, chainedL);
    neineiL = neiNeighborPos_C(-1, posMin, chainedL);


    if (clMin_11==1){
      d1 = distance_C(clMin_11, clMin_22, neiR[1], neiR[2], dissim, *h, *p);
      newDR = d1;
      CHAIN(MINCL1, jj) = CHAIN(MINCL1, posMin);
      CHAIN(MAXCL1, jj) = CHAIN(MAXCL2, posMin);
      CHAIN(MINCL2, jj) = neiR[1];
      CHAIN(MAXCL2, jj) = neiR[2];
      CHAIN(LAB1, jj) = step;
      CHAIN(LAB2, jj) = neiR[3];
      CHAIN(POSL, jj) = -1;
      CHAIN(POSR, jj) = neineiR;
      CHAIN(VALID, jj) = 1;
      if (neineiR>0){
    	CHAIN(POSL, neineiR) = jj;
      }
      CHAIN(VALID, posMin) = 0;
      CHAIN(VALID, neiR[0]) = 0;
      distances[jj-1] = newDR;
      positions = insertHeap_C(positions, distances, *lHeap, jj);
      *lHeap = *lHeap + 1;
      jj = jj+1;
    }
    else if (clMin_22==*p){
      d2 = distance_C(neiL[1], neiL[2], clMin_11, clMin_22, dissim, *h, *p);
      newDL = d2 ;
      CHAIN(MINCL1, jj) = neiL[1];
      CHAIN(MAXCL1, jj) = neiL[2];
      CHAIN(MINCL2, jj) = CHAIN(MINCL1, posMin);
      CHAIN(MAXCL2, jj) = CHAIN(MAXCL2, posMin);
      CHAIN(LAB1, jj) = neiL[3];
      CHAIN(LAB2, jj) = step;
      CHAIN(POSL, jj) = neineiL;
      CHAIN(POSR, jj) = -1;
      CHAIN(VALID, jj) = 1;
      if (neineiL>0){
    	CHAIN(POSR, neineiL) = jj;
      }
      CHAIN(VALID, posMin) = 0;
      CHAIN(VALID, neiL[0]) = 0;
      distances[jj-1] = newDL;
      positions = insertHeap_C(positions, distances, *lHeap, jj);
      *lHeap = *lHeap + 1;
      jj = jj+1;
    }
    else {
      d2 = distance_C(clMin_11, clMin_22, neiR[1], neiR[2], dissim, *h, *p);
      d1 = distance_C(neiL[1], neiL[2], clMin_11, clMin_22, dissim, *h, *p);
      newDR = d2;
      newDL = d1;
      CHAIN(MINCL1, jj) = neiL[1];
      CHAIN(MAXCL1, jj) = neiL[2];
      CHAIN(MINCL2, jj) = CHAIN(MINCL1, posMin);
      CHAIN(MAXCL2, jj) = CHAIN(MAXCL2, posMin);
      CHAIN(LAB1, jj) = neiL[3];
      CHAIN(LAB2, jj) = step;
      CHAIN(POSL, jj) = neineiL;
      CHAIN(POSR, jj) = jj+1;

      CHAIN(VALID, jj) = 1;
      CHAIN(MINCL1, jj+1) = CHAIN(MINCL1, posMin);
      CHAIN(MAXCL1, jj+1) = CHAIN(MAXCL2, posMin);
      CHAIN(MINCL2, jj+1) = neiR[1];
      CHAIN(MAXCL2, jj+1) = neiR[2];
      CHAIN(LAB1, jj+1) = step;
      CHAIN(LAB2, jj+1) = neiR[3];
      CHAIN(POSL, jj+1) = jj;
      CHAIN(POSR, jj+1) = neineiR;

      CHAIN(VALID, jj+1) = 1;
      if (neineiL>0){
     	CHAIN(POSR, neineiL) = jj;
      }
     if (neineiR>0){
       CHAIN(POSL, neineiR) = jj+1;
     }
     CHAIN(VALID, posMin) = 0;
     CHAIN(VALID, neiL[0]) = 0;
     CHAIN(VALID, neiR[0]) = 0;
     distances[jj-1] = newDL;
     positions = insertHeap_C(positions, distances, *lHeap, jj);
     *lHeap = *lHeap + 1;
     distances[jj] = newDR;
     positions = insertHeap_C(positions, distances, *lHeap, jj+1);
     *lHeap = *lHeap + 1;
     jj = jj+2;
    }

    stepInv = *p - step;
    MERGE(step, 1, k) = CHAIN(LAB1, posMin);
    MERGE(step, 2, k) = CHAIN(LAB2, posMin);
  }

  // merging the remaining two classes
  step = *p-1;
  posMin = jj-1;

  dLast = distance_C(CHAIN(MINCL1, posMin), CHAIN(MAXCL1, posMin), CHAIN(MINCL2, posMin), CHAIN(MAXCL2, posMin), dissim, *h, *p);


  gains[step-1] = dLast;
  MERGE(step, 1, k) = CHAIN(LAB1, jj-1);
  MERGE(step, 2, k) = CHAIN(LAB2, jj-1);

  UNPROTECT(6);
  return(Rmerge);
}

