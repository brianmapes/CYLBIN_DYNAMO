
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>

#include "rsl.h"

/* NOTE: 1. When using this program on a different machine, the only thing 
            needs to be changed is the output directory, which is stored in
            the variable cdir. Please find it in below and change it. 
         2. This program reads reflectivity (Z), velocity (Vr), and spectral 
            width (w). It works when there is not w value in the orginal data,
            but not work when either Z or Vr value is not available.
*/


/* Define all bin width and No of bins */

  #define AZ_BIN_WIDTH 15    /* degrees */
  #define NABINS 24

  #define R_BIN_WIDTH  8000 /* meters */
  #define NRBINS 12

  #define Z_BIN_WIDTH  500  /* meters */
  #define NZBINS 36

  #define DBZ_BIN_WIDTH  1   
  #define NDBZBINS 66        /* bad and 1 (or less) to 65 */

  #define SDZ_BIN_WIDTH  1  
  #define NSDZBINS 36       
  #define LENGTH_SDZ  5000   /* meters */

  #define VR_BIN_WIDTH  4    /* m/s */
  #define NVRBINS  20

  #define WID_BIN_WIDTH 0.5    /* m/s */
  #define NWIDBINS 40

  #define EARTH_R43  8.5e6  /* meters  4/3 earth redius*/
  #define NSIG_NO_ECHO -32


void main(int argc, char **argv)
{
/* Assign output subdirectory */
  char cdir[ ]="output/";

/* Output arrays */
  static unsigned long Zhist [NZBINS] [NRBINS] [NABINS] [NDBZBINS];
  static unsigned long SDZhist [NRBINS] [NABINS] [NDBZBINS] [NSDZBINS];
  static unsigned long Vhist [NZBINS] [NRBINS] [NABINS] [NVRBINS];
  static unsigned long Whist [NZBINS] [NRBINS] [NABINS] [NWIDBINS];
  static float Vsum  [NZBINS] [NRBINS] [NABINS];
  static float header [10];

/* Input and output file names */
  char ctemp[20],radar_name[5];
  int iyear,imon,iday,ihour,imin;
  static long int time;
  char filename0[80];
  char filename1[80],filename2[80],filename3[80],filename4[80],filename5[80];
  char filename6[80];
  FILE *outfile;
  char iptfile[100] = "/data/jlin/radar/data/RHB/RHB.19970808.0001.sig";
  int ifile;

/* Input data */
  Radar *radar;
  Volume *volumez;
  Ray *ray;
  int nsweeps;
  int isweep,iray,irange;

/* Coordinates */
  float elev, azimuth, range;
  double coselev,sinelev;
  float r,z;
  int i,j,k,itemp;

/* Physical values */
  float lon,lat;
  float pw;            /* Pulse width (micro-sec). */
  float prf;           /* Pulse repitition frequency, in Hz. */
  float wave;          /* Wavelength. Meters. */
  float nyq_vel;       /* Nyquist velocity. m/s */
  float rng_first_bin; /* Range to first gate.(meters) */
  float bin_space;     /* Data gate size (meters)*/
  float zray[2000];
  float avg, na, sd;
  int   id;
  float zvalue, vvalue, wvalue;
  int   zindex, sdzindex, vindex, windex;
  int   ifield_z0[4]={4,0,3,15};
  int   ifield_z, ifield;
  

/* read in filename */
  printf("\n input filename \n");
  scanf("%s",iptfile); 
  printf("\n%s \n",iptfile);
  printf("%f\n",BADVAL);

/* read data in this file */
  radar = RSL_anyformat_to_radar(iptfile, NULL); 

  ifield_z=-1;
  for(i=0; i<4; i++)
  {
    ifield=ifield_z0[i];
    if(radar->v[ifield] != NULL) ifield_z=ifield;
  }
  printf("ifield_z=%d\n",ifield_z);

  if(ifield_z < 0 || radar->v[1] == NULL)
  {
    printf("NULL data");
    exit(itemp);
  } 

  volumez = radar->v[ifield_z];
  nsweeps = volumez->h.nsweeps;
  header[3]=nsweeps;
  printf("nsweeps %i \n",nsweeps);
  if(nsweeps < 6) exit(itemp); 


/* get radar name and time, construct output file name */
  strcpy(ctemp,radar->h.radar_name);
  for (i=0; i<4; i++) {radar_name[i] = ctemp[i];}
  iyear = radar->h.year;
  imon = radar->h.month;
  iday = radar->h.day;
  ihour = radar->h.hour;
  imin = radar->h.minute;
  lon = radar->h.lond;
  lat = radar->h.latd;
  header[0]=imin;
  header[1]=lon;
  header[2]=lat;
  time = iyear*1000000+imon*10000+iday*100+ihour; 
  sprintf(filename0,"%s%s%ld",cdir,radar_name,time);
  sprintf(filename1,"%s.Zhist",filename0);
  sprintf(filename2,"%s.Whist",filename0);
  sprintf(filename3,"%s.Vhist",filename0);
  sprintf(filename4,"%s.Vsum",filename0);
  sprintf(filename5,"%s.SDZhist",filename0);
  sprintf(filename6,"%s.header",filename0);
  printf("%s\n",filename0);

/* if the output files already exist, read in the output arrays */
  outfile = fopen(filename1,"r");
  if(outfile != NULL) 
  {
    fread(Zhist,sizeof(Zhist),1,outfile);
    fclose(outfile);
    outfile = fopen(filename2,"r");
    fread(Whist,sizeof(Whist),1,outfile);
    fclose(outfile);
    outfile = fopen(filename3,"r");
    fread(Vhist,sizeof(Vhist),1,outfile);
    fclose(outfile);
    outfile = fopen(filename4,"r");
    fread(Vsum,sizeof(Vsum),1,outfile);
    fclose(outfile);
    outfile = fopen(filename5,"r");
    fread(SDZhist,sizeof(SDZhist),1,outfile);
    fclose(outfile);
  }

/* loop through the data and add to the output arrays */

  for (isweep=0; isweep<nsweeps; isweep++) 
  {
    if (volumez->sweep[isweep] == NULL) continue;

    printf("isweep %i nrays %i \n",isweep,volumez->sweep[isweep]->h.nrays);
    for (iray=0; iray<volumez->sweep[isweep]->h.nrays; iray++) 
    {
      ray = volumez->sweep[isweep]->ray[iray];
      if (ray == NULL) continue;

      elev = ray->h.elev;
      sinelev = sin(elev/180.0 *3.14159);
      coselev = cos(elev/180.0 *3.14159);
      rng_first_bin = ray->h.range_bin1;
      bin_space = ray->h.gate_size;
      pw = ray->h.pulse_width;
      prf = ray->h.prf;
      wave = ray->h.wavelength;
/*      nyq_vel = ray->h.nyq_vel; */
      nyq_vel = wave*prf/4.0;
      header[4]=rng_first_bin;
      header[5]=bin_space;
      header[6]=pw;
      header[7]=prf;
      header[8]=wave;
      header[9]=nyq_vel;

      azimuth = ray->h.azimuth;
      j = (int)(azimuth/AZ_BIN_WIDTH);
      if(j < 0 || j > NABINS-1) continue;

      ray = volumez->sweep[isweep]->ray[iray];
      for (irange=0; irange<ray->h.nbins; irange++) 
      {
        zvalue = ray->h.f(ray->range[irange]);
        if(zvalue > BADVAL-1 || zvalue < NSIG_NO_ECHO+2) zvalue = -9999.0;
        zray[irange] = zvalue;
      }

      id = LENGTH_SDZ/bin_space;
      for (irange=0; irange<ray->h.nbins; irange++) 
      {
        range = rng_first_bin + irange*bin_space;
        r = range*coselev;
        z = range*sinelev + r*r /EARTH_R43;

        i = (int)(r/R_BIN_WIDTH);
        k = (int)(z/Z_BIN_WIDTH);
        if(i < 0 || i > NRBINS-1) 
	{
          break;
	  continue;
	}
        if(k < 0 || k > NZBINS-1) 
	{
          break;
	  continue;
	}

        ray = volumez->sweep[isweep]->ray[iray];
        zvalue = ray->h.f(ray->range[irange]);
        ray = radar->v[1]->sweep[isweep]->ray[iray];
        vvalue = ray->h.f(ray->range[irange]);
        if(radar->v[2] == NULL)
          wvalue = 99.0;
        else
        {
          ray = radar->v[2]->sweep[isweep]->ray[iray];
          wvalue = ray->h.f(ray->range[irange]);
        }

/* now index values: Z */
          if(zvalue < 2) 
            zindex = 1;
          else 
          {
            zindex = zvalue;
            if (zvalue == BADVAL) zindex=0;
            if (zindex > NDBZBINS-1) zindex=NDBZBINS-1;
          }
/* now index values: V */
          vindex = (nyq_vel + vvalue)/VR_BIN_WIDTH +1;
          if (vvalue == BADVAL) vindex=0;
          if (vindex > NVRBINS-1) vindex=NVRBINS-1;
/* now index values: width */
          windex = (int)(wvalue/WID_BIN_WIDTH);
          if (wvalue == BADVAL) windex=0;
          if (windex > NWIDBINS-1) windex = NWIDBINS-1;

/* BAD index values, sum only good velocities */
          if( (vvalue <= -nyq_vel) || (zvalue < NSIG_NO_ECHO+2))
          {
            zindex = 0; 
            vindex = 0; 
            windex = NWIDBINS-1;
          }

          if(zindex < 0 || zindex > NDBZBINS-1 ||
             vindex < 0 || vindex > NVRBINS-1 ||
             windex < 0 || windex > NWIDBINS-1) 
            continue;
          else
          {
            Zhist[k][i][j][zindex]++;
            Vhist[k][i][j][vindex]++;
            Whist[k][i][j][windex]++;
            if(vindex >= 1) Vsum [k][i][j] += vvalue;
          }

/* now index values: SDZ only between 2-4 km */
          if(k>3 && k<9)
          {
            avg = 0.0;
            na = 0.0;
            for (itemp=irange-id; itemp<irange+id+1; itemp++)
            {
              if(zray[itemp] > -9990.0)
              {
                avg = avg + zray[itemp];
                na = na + 1.0;
              }
            }

            if(na > 5)
            {
              avg = avg/na;
              sd = 0.0;
              for (itemp=irange-id; itemp<irange+id+1; itemp++)
              {
                if(zray[itemp] > -9990.0)
                {
                  sd = sd + pow((zray[itemp]-avg),2);
                }
              }
              sd = sqrt(sd/(na-1.0));
            }
            else
              sd = -9999.0;

            if(sd > -9990.0)
            {
              sdzindex = (int)(sd/SDZ_BIN_WIDTH)+1;
              if(sdzindex > NSDZBINS-1) sdzindex = NSDZBINS-1;
            }
            else
              sdzindex = 0;

            SDZhist[i][j][zindex][sdzindex]++;
          }

      };
    };
  };

  printf("nyq_vel= %f \n",nyq_vel);
  outfile = fopen(filename6,"a");
  fwrite(header,sizeof(header),1,outfile);
  fclose(outfile);

  outfile = fopen(filename1,"w");
  fwrite(Zhist,sizeof(Zhist),1,outfile);
  fclose(outfile);
  outfile = fopen(filename5,"w");
  fwrite(SDZhist,sizeof(SDZhist),1,outfile);
  fclose(outfile);
  outfile = fopen(filename2,"w");
  fwrite(Whist,sizeof(Whist),1,outfile);
  fclose(outfile);
  outfile = fopen(filename3,"w");
  fwrite(Vhist,sizeof(Vhist),1,outfile);
  fclose(outfile);
  outfile = fopen(filename4,"w");
  fwrite(Vsum,sizeof(Vsum),1,outfile);
  fclose(outfile);

}
