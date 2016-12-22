/* C modules for I/O etc. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <unistd.h>


#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;


#ifdef UPCASE
void BWRITE_HFET 
(int *n, double *f, double *fd, double *ft, double *fdd, double *ftt,
 double *fdt, double *fddt, double *fdtt, double *fddtt, double *dpdf, double *dpdfd,
 double *dpdft, double *dpdfdt, double *ef, double *efd, double *eft, double *efdt, 
 double *xf, double *xfd, double *xft, double *xfdt)
#elif IBM
void bwrite_hfet 
(int *n, double *f, double *fd, double *ft, double *fdd, double *ftt,
 double *fdt, double *fddt, double *fdtt, double *fddtt, double *dpdf, double *dpdfd,
 double *dpdft, double *dpdfdt, double *ef, double *efd, double *eft, double *efdt, 
 double *xf, double *xfd, double *xft, double *xfdt)
#else
void bwrite_hfet_
(int *n, double *f, double *fd, double *ft, double *fdd, double *ftt,
 double *fdt, double *fddt, double *fdtt, double *fddtt, double *dpdf, double *dpdfd,
 double *dpdft, double *dpdfdt, double *ef, double *efd, double *eft, double *efdt, 
 double *xf, double *xfd, double *xft, double *xfdt)
#endif
{
  FILE *fp;
  int nr;


#ifdef SAFE
  if (*n<0)
    {printf("bwrite_hfet_() :: n must be positive\n"); abort();}
#endif

  if (!(fp=fopen("helm_table.bdat","wb")))
    {printf("bwrite_hfet_() :: fopen failure!\n"); abort();}

  nr=fwrite(f,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on f!\n"); abort();}
  nr=fwrite(fd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on fd!\n"); abort();}
  nr=fwrite(ft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on ft!\n"); abort();}
  nr=fwrite(fdd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on fdd!\n"); abort();}
  nr=fwrite(ftt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on ftt!\n"); abort();}
  nr=fwrite(fdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on fdt!\n"); abort();}
  nr=fwrite(fddt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on fddt!\n"); abort();}
  nr=fwrite(fdtt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on fdtt!\n"); abort();}
  nr=fwrite(fddtt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on fddtt!\n"); abort();}
  nr=fwrite(dpdf,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on dpdf!\n"); abort();}
  nr=fwrite(dpdfd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on dpdfd!\n"); abort();}
  nr=fwrite(dpdft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on dpdft!\n"); abort();}
  nr=fwrite(dpdfdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on dpdfdt!\n"); abort();}
  nr=fwrite(ef,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on ef!\n"); abort();}
  nr=fwrite(efd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on efd!\n"); abort();}
  nr=fwrite(eft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on eft!\n"); abort();}
  nr=fwrite(efdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on efdt!\n"); abort();}
  nr=fwrite(xf,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on xf!\n"); abort();}
  nr=fwrite(xfd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on xfd!\n"); abort();}
  nr=fwrite(xft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on xft!\n"); abort();}
  nr=fwrite(xfdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bwrite_hfet() :: failed write on xfdt!\n"); abort();}

  if (fclose(fp))
    {printf("bwrite_hfet() :: couldn't close file!\n"); abort();}
}


#ifdef UPCASE
void BREAD_HFET 
(int *n, double *f, double *fd, double *ft, double *fdd, double *ftt,
 double *fdt, double *fddt, double *fdtt, double *fddtt, double *dpdf, double *dpdfd,
 double *dpdft, double *dpdfdt, double *ef, double *efd, double *eft, double *efdt, 
 double *xf, double *xfd, double *xft, double *xfdt)
#elif IBM
void bread_hfet 
(int *n, double *f, double *fd, double *ft, double *fdd, double *ftt,
 double *fdt, double *fddt, double *fdtt, double *fddtt, double *dpdf, double *dpdfd,
 double *dpdft, double *dpdfdt, double *ef, double *efd, double *eft, double *efdt, 
 double *xf, double *xfd, double *xft, double *xfdt)
#else
void bread_hfet_
(int *n, double *f, double *fd, double *ft, double *fdd, double *ftt,
 double *fdt, double *fddt, double *fdtt, double *fddtt, double *dpdf, double *dpdfd,
 double *dpdft, double *dpdfdt, double *ef, double *efd, double *eft, double *efdt, 
 double *xf, double *xfd, double *xft, double *xfdt)
#endif
{
  FILE *fp;
  int nr;


#ifdef SAFE
  if (*n<0)
    {printf("bread_hfet_() :: n must be positive\n"); abort();}
#endif

  if (!(fp=fopen("helm_table.bdat","rb")))
    {printf("bread_hfet_() :: fopen failure!\n"); abort();}

  nr=fread(f,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on f!\n"); abort();}
  nr=fread(fd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on fd!\n"); abort();}
  nr=fread(ft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on ft!\n"); abort();}
  nr=fread(fdd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on fdd!\n"); abort();}
  nr=fread(ftt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on ftt!\n"); abort();}
  nr=fread(fdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on fdt!\n"); abort();}
  nr=fread(fddt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on fddt!\n"); abort();}
  nr=fread(fdtt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on fdtt!\n"); abort();}
  nr=fread(fddtt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on fddtt!\n"); abort();}
  nr=fread(dpdf,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on dpdf!\n"); abort();}
  nr=fread(dpdfd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on dpdfd!\n"); abort();}
  nr=fread(dpdft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on dpdft!\n"); abort();}
  nr=fread(dpdfdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on dpdfdt!\n"); abort();}
  nr=fread(ef,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on ef!\n"); abort();}
  nr=fread(efd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on efd!\n"); abort();}
  nr=fread(eft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on eft!\n"); abort();}
  nr=fread(efdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on efdt!\n"); abort();}
  nr=fread(xf,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on xf!\n"); abort();}
  nr=fread(xfd,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on xfd!\n"); abort();}
  nr=fread(xft,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on xft!\n"); abort();}
  nr=fread(xfdt,sizeof(double),*n,fp);
  if (nr!=*n)
    {printf("bread_hfet() :: failed read on xfdt!\n"); abort();}

  if (fclose(fp))
    {printf("bread_hfet() :: couldn't close file!\n"); abort();}
}



#ifdef UPCASE
void
BYTE_REVERSE (float *buf, int *nn)
#elif IBM
void
byte_reverse (float *buf, int *nn)
#else
void
byte_reverse_(float *buf, int *nn)
#endif
{
  int n;
  char temp, *ptr;


#ifdef SAFE
  if (*nn<0)
    {printf("byte_reverse() :: n must be positive\n"); abort();}
#endif
  
  for (ptr=(char *)buf,n=*nn; n--; ptr+=4)
    {
      SWAP(ptr[0],ptr[3])
      SWAP(ptr[1],ptr[2])
     }
}
