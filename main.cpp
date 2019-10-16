#include "vars.h"

int numx;
int numy;
int numz;
int nump;
int numhalo;
int numbound;
double h;

double timetot;
double timepre;

int batchsize;
int numbatch;
int numiter;

int *boundmap;
int *halodispl;
int *haloind;
double *haloval;

int *displ;
int *ind;
int numind;
double *val;
double *diag;
int *batchdispl;

double *posx;
double *posy;
double *posz;
double *poshalox;
double *poshaloy;
double *poshaloz;

double *ref;

int main(int argc, char** argv){

  int numproc;
  int myid;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  MPI_Barrier(MPI_COMM_WORLD);
  timetot = MPI_Wtime();

  int numthreads;
  #pragma omp parallel
  if(omp_get_thread_num()==0)numthreads=omp_get_num_threads();

  if(myid==0){
    printf("NUMBER OF PROCESSES    : %d\n",numproc);
    printf("NUMBER OF THRD./PROC.  : %d\n",numthreads);
  }

  //SCANNING GEOMETRY DATA
  char *chartemp;
  chartemp = getenv("NUMX");
  numx = atoi(chartemp);
  chartemp = getenv("NUMY");
  numy = atoi(chartemp);
  chartemp = getenv("NUMZ");
  numz = atoi(chartemp);
  chartemp = getenv("DISC");
  h = strtod(chartemp,NULL);
  chartemp = getenv("NUMITER");
  numiter = atoi(chartemp);
  chartemp = getenv("BATCHSIZE");
  batchsize = atoi(chartemp);

  //INITIAL DATA
  nump = numx*numy*numz;
  numbatch = nump/batchsize;
  if(nump%batchsize)numbatch++;
  numhalo = numx*numy*2+numx*numz*2+numy*numz*2+numx*4+numy*4+numz*4+8;
  numbound = numx*numy*2+numx*numz*2+numy*numz*2-numx*4-numy*4-numz*4+8;

  //PRINT PARAMETERS
  if(myid==0){
    printf("GRID (X,Y,Z)            : %dx%dx%d\n",numx,numy,numz);
    printf("NUMBER OF POINTS        : %d\n",nump);
    printf("NUMBER OF HALO          : %d\n",numhalo);
    printf("NUMBER OF BOUNDARIES    : %d\n",numbound);
    printf("DISCRETIZATION          : %.3e\n",h);
    printf("DIMENSION (X,Y,Z)       : %.3ex%.3ex%.3e\n",numx*h,numy*h,numz*h);
    printf("\n");
    printf("NUMBER OF INTERACTIONS  : %d\n",nump*27);
    printf("GRID MEMORY             : %f MB\n",nump/1024.0/1024.0*sizeof(double));
    printf("MATRIX MEMORY           : %f GB\n",nump/1024.0/1024.0/1024.0*(sizeof(double)+sizeof(int))*27.0);
    printf("HALO MEMORY             : %f MB\n",numhalo/1024.0/1024.0*sizeof(double));
    printf("HALO MATRIX MEMORY      : %f MB\n",numhalo/1024.0/1024.0*(sizeof(double)+sizeof(int))*9.0);
    printf("\n");
    printf("NUMBER OF ITERATIONS : %d\n",numiter);
    printf("BATCH SIZE           : %d\n",batchsize);
    printf("BUMBER OF BATHCES    : %d\n",numbatch);
  }
  double matsize = nump/1024.0/1024.0/1024.0*(sizeof(double)+sizeof(int))*27.0;
  MPI_Barrier(MPI_COMM_WORLD);
  timepre = MPI_Wtime();
  
  posx = new double[nump];
  posy = new double[nump];
  posz = new double[nump];
  poshalox = new double[numhalo];
  poshaloy = new double[numhalo];
  poshaloz = new double[numhalo];
  setpositions();

  //LOOKUP MATRIX
  int *lookind = new int[27];
  double *lookval = new double[27];
  {
    /*double w0 = -6;
    double wf = 1;
    double we = 0;
    double wc = 0;
    double w = 1/(h*h);*/
    double w0 = -24;
    double wf = 2;
    double we = 1;
    double wc = 0;
    double w = 1/(6*h*h);
    /*double w0 = -128;
    double wf = 14;
    double we = 3;
    double wc = 1;
    double w = 1/(30*h*h);*/
   
    int count = 0;
    for(int kz = -1; kz < 2; kz++)
      for(int ky = -1; ky < 2; ky++)
        for(int kx = -1; kx < 2; kx++){
          int diff = kz*numx*numy+ky*numx+kx;
          lookind[count] = diff;
          lookval[count] = 0;
          //SELF
          if(kx==0 && ky==0 && kz==0)
            lookval[count] = w0;
          //FACE
          else if((ky==0&&kz==0)||
                  (kx==0&&kz==0)||
                  (kx==0&&ky==0))
            lookval[count] = wf;
          //EDGE
          else if(kx==0 || ky==0 || kz==0)
            lookval[count] = we;
          //CORNER
          else
            lookval[count] = wc;
          count++;
        }
    if(myid == 0){
      printf("\n");
      printf("PRINT 27-POINT LOOKUP TABLE\n");
      for(int n = 0; n < 27; n++)
        printf("%d: %d \t %.0f\n",n,lookind[n],lookval[n]);
      double sum = 0;
      for(int n = 0; n < 27; n++)
        sum += lookval[n];
      printf("sum: %e\n",sum);
    }
    for(int n = 0; n < 27; n++)
      lookval[n] = lookval[n]*w;
  }

  //FILL BOUNDARY INDICES & VALUES
  boundmap = new int[numbound];
  halodispl = new int[numbound+1];
  {
    int *boundnz = new int[numbound];
    int count = 0;
    for(int m = 0; m < nump; m++)
      if(posx[m] < h || posy[m] < h || posz[m] < h ||
         posx[m] > (numx-1)*h || posy[m] > (numy-1)*h || posz[m] > (numz-1)*h){
        boundmap[count] = m;
        boundnz[count] = 0;
        for(int n = 0; n < numhalo; n++){
          double dist = sqrt(pow(posx[m]-poshalox[n],2)+pow(posy[m]-poshaloy[n],2)+pow(posz[m]-poshaloz[n],2));
          if(dist<1.8*h)
            boundnz[count]++;
        }
        count++;
      }
    halodispl[0] = 0;
    for(int n = 1; n < numbound+1; n++)
      halodispl[n] = halodispl[n-1]+boundnz[n-1];
    int boundtot = halodispl[numbound];
    if(myid==0)printf("\nBOUNDARY-HALO INTERACTIONS: %d\n",boundtot);
    haloind = new int[boundtot];
    haloval = new double[boundtot];
    for(int m = 0; m < numbound; m++){
      boundnz[m] = 0;
      for(int n = 0; n < numhalo; n++){
        double dist = sqrt(pow(posx[boundmap[m]]-poshalox[n],2)+pow(posy[boundmap[m]]-poshaloy[n],2)+pow(posz[boundmap[m]]-poshaloz[n],2));
        if(dist<1.8*h){
          haloind[halodispl[m]+boundnz[m]] = n;
          int xdif = round((poshalox[n]-posx[boundmap[m]])/h);
          int ydif = round((poshaloy[n]-posy[boundmap[m]])/h);
          int zdif = round((poshaloz[n]-posz[boundmap[m]])/h);
          int ind = (xdif+1)+(ydif+1)*3+(zdif+1)*9;
          haloval[halodispl[m]+boundnz[m]] = lookval[ind];
          boundnz[m]++;
        }
      }
    }
    /*if(myid==0){
      for(int m = 0; m < numbound; m++){
        printf("bound %d: %d <- ",m,boundmap[m]);
        for(int n = halodispl[m]; n < halodispl[m+1]; n++)
          printf("%d ",haloind[n]);
        printf("\n");
      }
      printf("\n");
      for(int m = 0; m < numbound; m++){
        printf("bound %d: %d <- ",m,boundmap[m]);
        for(int n = halodispl[m]; n < halodispl[m+1]; n++)
          printf("%.0f ",haloval[n]);
        printf("\n");
      }
    }*/
  }
  //FILL INDICES & VALUES
  {
    int *numnz = new int[nump];
    for(int m = 0; m < nump; m++){
      numnz[m] = 0;
      for(int n = 0; n < 27; n++){
        int k = m+lookind[n];
        if(k>-1 && k < nump){
          double dist = sqrt(pow(posx[m]-posx[k],2)+pow(posy[m]-posy[k],2)+pow(posz[m]-posz[k],2));
          if(dist<1.8*h)
            numnz[m]++;
        }
      }
    }
    displ = new int[nump+1];
    displ[0] = 0;
    for(int m = 1; m < nump+1; m++)
      displ[m] = displ[m-1]+numnz[m-1];
    numind = displ[nump];
    if(myid==0)printf("\nNUMBER OF NON-ZEROES: %d\n",numind);
    ind = new int[numind];
    val = new double[numind];
    for(int m = 0; m < nump; m++){
      numnz[m] = 0;
      for(int n = 0; n < 27; n++){
        int k = m+lookind[n];
        if(k>-1 && k < nump){
          double dist = sqrt(pow(posx[m]-posx[k],2)+pow(posy[m]-posy[k],2)+pow(posz[m]-posz[k],2));
          if(dist<1.8*h){
            ind[displ[m]+numnz[m]]=k;
            val[displ[m]+numnz[m]]=lookval[n];
            numnz[m]++;
          }
        }
      }
    }
    /*if(myid == 0){
      for(int m = 0; m < nump; m++){
        printf("inter %d <- ",m);
        for(int n = displ[m]; n < displ[m+1]; n++)
          printf("%d ",ind[n]);
        printf("\n");
      }
      printf("\n");
      for(int m = 0; m < nump; m++){
        printf("inter %d <- ",m);
        for(int n = displ[m]; n < displ[m+1]; n++)
          printf("%.0f ",val[n]);
        printf("\n");
      }
      printf("\n");
      int *mat = new int[nump*nump];
      for(int n = 0; n < nump*nump; n++)
        mat[n] = 0;
      for(int m = 0; m < nump; m++)
        for(int n = displ[m]; n < displ[m+1]; n++)
          mat[m*nump+ind[n]] = val[n];
      for(int m = 0; m < nump; m++){
        for(int n = 0; n < nump; n++){
          printf("%d ",abs((int)mat[m*nump+n]));
          if(n%3==2)printf("\n");
          if(n%9==8)printf("\n");
        }
        printf("\n\n");
      }
    }*/
  }
  diag = new double[nump];
  for(int m = 0; m < nump; m++)
    for(int n = displ[m]; n < displ[m+1]; n++)
      if(ind[n]==m){
        diag[m] = val[n];
        break;
      }
  MPI_Barrier(MPI_COMM_WORLD);
  timepre = MPI_Wtime()-timepre;
  if(myid==0)printf("\nPREPROCESSING: %e s\n",timepre);

  //BATCH REORDERING
  int *batchmap;
  int *invbatchmap;
  {
    //BREADTH-FIRST SEARCH
    int *oldfront = new int[nump];
    int *newfront = new int[nump];
    int *tagfront = new int[nump];
    int *tagbatch = new int[nump];
    for(int n = 0; n < nump; n++)
      tagbatch[n] = -1;
    for(int n = 0; n < nump; n++)
      tagfront[n] = -1;
    int numbatch = 0;
    int numfront = 0;
    int seed = 0;
    while(true){
      oldfront[0] = seed;
      tagbatch[seed] = numbatch;
      tagfront[seed] = numfront;
      int numoldfront = 1;
      int numnewfront = 0;
      int sizetemp = 1;
      //printf("seed: %d\n",seed);
      while(true){
        numnewfront = 0;
        numfront++;
        for(int m = 0; m < numoldfront; m++){
          int old = oldfront[m];
          for(int n = displ[old]; n < displ[old+1]; n++){
            int con = ind[n];
            if(tagfront[con]==-1 && sizetemp < batchsize){
              newfront[numnewfront] = con;
              tagfront[con] = numfront;
              tagbatch[con] = numbatch;
              numnewfront++;
              sizetemp++;
            }
          }
        }
        if(numnewfront==0)break;

        numoldfront = numnewfront;
        int *temp = oldfront;
        oldfront = newfront;
        newfront = temp;
      }
      numbatch++;
      int newseed = -1;
      for(int n = seed; n < nump; n++)
        if(tagbatch[n]==-1){
          newseed = n;
          break;
        }
      if(newseed==-1)
        break;
      else
        seed = newseed;
    }
    delete[] oldfront;
    delete[] newfront;
    delete[] tagfront;
    if(myid==0)printf("number of batch / front: %d %d\n",numbatch,numfront);

    //MAPPING
    batchmap = new int[nump];
    int *batchnz = new int[numbatch];
    batchdispl = new int[numbatch+1];
    for(int n = 0; n < numbatch; n++)
      batchnz[n] = 0;
    for(int n = 0; n < nump; n++)
      batchnz[tagbatch[n]]++;
    batchdispl[0] = 0;
    for(int n = 1; n < numbatch+1; n++)
      batchdispl[n] = batchdispl[n-1]+batchnz[n-1];
    for(int n = 0; n < numbatch; n++)
      batchnz[n] = 0;
    for(int n = 0; n < nump; n++){
      batchmap[batchdispl[tagbatch[n]]+batchnz[tagbatch[n]]] = n;
      batchnz[tagbatch[n]]++;
    }
    invbatchmap = new int[nump];
    for(int n = 0; n < nump; n++)
      invbatchmap[batchmap[n]] = n;

    //REORDERING
    int *newnumnz = new int[nump];
    int *newdispl = new int[nump+1];
    int *newind = new int[numind];
    double *newval = new double[numind];
    for(int n = 0; n < nump; n++)
      newnumnz[n] = displ[batchmap[n]+1]-displ[batchmap[n]];
    newdispl[0] = 0;
    for(int n = 1; n < nump+1; n++)
      newdispl[n] = newdispl[n-1]+newnumnz[n-1];
    for(int m = 0; m < nump; m++)
      for(int n = 0; n < newnumnz[m]; n++){
        newind[newdispl[m]+n] = invbatchmap[ind[displ[batchmap[m]]+n]];
        newval[newdispl[m]+n] = val[displ[batchmap[m]]+n];
      }
    delete[] batchnz;
    delete[] batchdispl;
    delete[] displ;
    delete[] ind;
    delete[] val;
    displ = newdispl;
    ind = newind;
    val = newval;
    for(int m = 0; m < nump; m++)
      for(int n = displ[m]; n < displ[m+1]; n++)
        if(ind[n]==m){
          diag[m] = val[n];
        }
  }
  /*//BUFFERING
  int *mapnz;
  int *mapdispl;
  int *map;
  unsigned short *sndex;
  //int *sndex;
  {
    mapnz = new int[numbatch];
    int *batchtag = new int[nump];
    for(int n = 0; n < nump; n++)
      batchtag[n] = -1;
    for(int batch = 0; batch < numbatch; batch++){
      mapnz[batch] = 0;
      for(int m = batchdispl[batch]; m < batchdispl[batch+1]; m++)
        for(int n = displ[m]; n < displ[m+1]; n++)
          if(ind[n]/batchsize!=batch && batchtag[ind[n]]!=batch){
            batchtag[ind[n]]=batch;
            mapnz[batch]++;
          }
    }
    int maxnz = 0;
    for(int n = 0; n < numbatch; n++)
      if(mapnz[n]>maxnz)maxnz=mapnz[n];
    if(myid==0){
      printf("\n");
      printf("batchsize: %d (%f KB)\n",batchsize,batchsize*sizeof(double)/1024.0);
      printf("maxnz: %d max buffer size: %d (%f KB)\n",maxnz,batchsize+maxnz,(batchsize+maxnz)*sizeof(double)/1024.0);
    }
    //if((batchsize+maxnz)>2*32768)return 0;
    mapdispl = new int[numbatch+1];
    mapdispl[0] = 0;
    for(int n = 1; n < numbatch+1; n++)
      mapdispl[n] = mapdispl[n-1]+mapnz[n-1];
    int maptot = mapdispl[numbatch];
    if(myid==0){
      printf("maptot: %d/%d (%.2f)\n",maptot,nump,maptot/(float)nump);
    }
    map = new int[maptot];
    for(int n = 0; n < nump; n++)
      batchtag[n] = -1;
    for(int batch = 0; batch < numbatch; batch++){
      mapnz[batch] = 0;
      for(int m = batchdispl[batch]; m < batchdispl[m+1]; m++)
        for(int n = displ[m]; n < displ[m+1]; n++)
          if(ind[n]/batchsize!=batch && batchtag[ind[n]]!=batch){
            batchtag[ind[n]]=batch;
            map[mapdispl[batch]+mapnz[batch]]=ind[n];
            mapnz[batch]++;
          }
    }
    delete[] mapnz;
    //for(int m = 0; m < numbatch; m++){
    //  printf("batch %d ****\n",m);
    //  for(int n = mapdispl[m]; n < mapdispl[m+1]; n++)
    //    printf("  %d\n",map[n]);
    //}
    //sndex = new unsigned short[numindex];
    sndex = new unsigned short[numind];
    for(int batch = 0; batch < numbatch; batch++){
      int count = 0;
      for(int m = batchdispl[batch]; m < batchdispl[batch+1]; m++){
        batchtag[m] = count;
        count++;
      }
      for(int n = mapdispl[batch]; n < mapdispl[batch+1]; n++){
        batchtag[map[n]]=count;
        count++;
      }
      for(int m = batchdispl[batch]; m < batchdispl[batch+1]; m++)
        for(int n = displ[m]; n < displ[m+1]; n++)
          sndex[n] = batchtag[ind[n]];
    }
  }*/

  //INITALIZE
  double *halo = new double[numhalo];
  double *f = new double[nump];
  double *x = new double[nump];
  double *b = new double[nump];

  /*for(int n = 0; n < nump; n++)
    x[n] = sin(posx[n]+posy[n]+posz[n])*cos(posx[n]+posy[n]+posz[n]);
  for(int n = 0; n < numhalo; n++)
    halo[n] = sin(poshalox[n]+poshaloy[n]+poshaloz[n])*cos(poshalox[n]+poshaloy[n]+poshaloz[n]);
  for(int n = 0; n < nump; n++){
    f[n] = -6*sin(2*(posx[n]+posy[n]+posz[n]));
    f[n] += (h*h)/12.0*72*sin(2*(posx[n]+posy[n]+posz[n]));
    //double a = -96*sin(2*(posx[n]+posy[n]+posz[n]));
    //f[n] += (h*h*h*h)/120*15*a;
  }

  double time = MPI_Wtime();
  for(int m = 0; m < nump; m++){
    double reduce = 0;
    for(int n = displ[m]; n < displ[m+1]; n++)
      reduce += val[n]*x[ind[n]];
    b[m] = reduce;
  }
  if(myid==0)printf("spmv time: %e\n",MPI_Wtime()-time);
  time = MPI_Wtime();
  for(int m = 0; m < numbound; m++){
    double reduce = 0;
    for(int n = halodispl[m]; n < halodispl[m+1]; n++)
      reduce += haloval[n]*halo[haloind[n]];
    b[boundmap[m]] += reduce;
  }
  if(myid==0)printf("halo time: %e\n",MPI_Wtime()-time);

  double l2 = 0;
  for(int n = 0; n < nump; n++)
    l2 += norm(b[n]-f[n]);
  l2 = sqrt(l2/nump);
  double l0 = abs(b[0]-f[0]);
  for(int n = 0; n < nump; n++)
    if(l0 < abs(b[n]-f[n]))
      l0 = abs(b[n]-f[n]);
  if(myid==0)printf("L0: %e L2: %e\n",l0,l2);*/

  //FORCING FUNCTION
  for(int n = 0; n < nump; n++){
    f[n] = -6*sin(2*(posx[n]+posy[n]+posz[n]));
    f[n] += (h*h)/12.0*72*sin(2*(posx[n]+posy[n]+posz[n]));
  }
  //BOUNDARY CONDITION
  for(int n = 0; n < numhalo; n++)
    halo[n] = sin(poshalox[n]+poshaloy[n]+poshaloz[n])*cos(poshalox[n]+poshaloy[n]+poshaloz[n]);
  //CONTRIBUTION FROM HALO
  for(int m = 0; m < numbound; m++){
    double reduce = 0;
    for(int n = halodispl[m]; n < halodispl[m+1]; n++)
      reduce += haloval[n]*halo[haloind[n]];
    f[boundmap[m]] -= reduce;
  }
  //CONTRIBUTION FROM RHS
  for(int n = 0; n < nump; n++)
    b[invbatchmap[n]] = f[n];
  ref = new double[nump];
  for(int n = 0; n < nump; n++)
    ref[invbatchmap[n]] = sin(posx[n]+posy[n]+posz[n])*cos(posx[n]+posy[n]+posz[n]);

  double *temp = new double[nump];
  if(myid==0){
    FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/f.bin","wb");
    fwrite(ref,sizeof(double),nump,filef);
    fclose(filef);
  }
  if(myid==0){
    FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/b.bin","wb");
    fwrite(b,sizeof(double),nump,filef);
    fclose(filef);
  }
  //JACOBI METHOD
  jacobi(x,b);
  if(myid==0){
    //for(int n = 0; n < nump; n++)
    //  temp[batchmap[n]] = x[n];
    FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/jacobi.bin","wb");
    //fwrite(temp,sizeof(double),nump,filef);
    fwrite(x,sizeof(double),nump,filef);
    fclose(filef);
  }
  //GAUSS-SEIDEL METHOD
  gaussseidel(x,b);
  if(myid==0){
    //for(int n = 0; n < nump; n++)
    //  temp[batchmap[n]] = x[n];
    FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/gaussseidel.bin","wb");
    //fwrite(temp,sizeof(double),nump,filef);
    fwrite(x,sizeof(double),nump,filef);
    fclose(filef);
  }
  /*//GRADIENT-DESCENT METHOD
  gradientdescent(x,b);
  if(myid==0){
    //for(int n = 0; n < nump; n++)
    //  temp[batchmap[n]] = x[n];
    FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/gradientdescent.bin","wb");
    //fwrite(temp,sizeof(double),nump,filef);
    fwrite(x,sizeof(double),nump,filef);
    fclose(filef);
  }
  //CONJUGATE-GRADIENT METHOD
  conjugategradient(x,b);
  if(myid==0){
    //for(int n = 0; n < nump; n++)
    //  temp[batchmap[n]] = x[n];
    FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/conjugategradient.bin","wb");
    //fwrite(temp,sizeof(double),nump,filef);
    fwrite(x,sizeof(double),nump,filef);
    fclose(filef);
  }*/

  MPI_Barrier(MPI_COMM_WORLD);
  timetot = MPI_Wtime()-timetot;
  if(myid==0)printf("TOTAL TIME: %e s\n",timetot);

  MPI_Finalize();
}
