// DUMPHALO.CPP
// by Andrew Pontzen
//
// Dump a halo from GADGET/FoF_Special output

#include <siman/base.hpp>
#include <siman/extra.hpp>

int main(int argc, char **argv)
{


  char basepath[1024];
  char fof_catalogue_path[1024];
  char fof_indexlist_path[1024];
 
  char output[1024];
  char snapname[1024];


  int halonum = 0;
  int snapshot = 0;
  int verbose =0;


  int *particle_list = NULL;
  int num_particles = 0;

  int n;

  if(argc < 6) {

    printf("dumphalo particle_data_path snap_basename snapshotnum halonumbermin halonumbermax [verbose] [live]\n");
    printf("  particle_data_path = path to particledata (other paths inferred)\n");
    printf("  snap_basename = snapshot basename (e.g. snap, snap_DM etc)\n");
    printf("  snapshot number   = gadget snapshot number >=0\n");
    printf("  halonumbermin/max   = halo number min/max \n");
    printf("\n\n");
    printf("Output (to stdout) consists of columns as follows:\n");
    printf("x/kpc\ty/kpc\tz/kpc\tvx/km s^-1\tvy\tvz\tmass/M_sol\ttype\n");
    printf("\nAll units are PHYSICAL\n\n");
    exit(0);
  }

  strcpy(basepath,argv[1]);
  strcpy(snapname,argv[2]);
 
  snapshot = atoi(argv[3]);
  halonum = atoi(argv[4]);

  bool live_from_disk = false;

  CGadgetFile gadget(basepath,snapname,snapshot);
  //CSimSnap *pSf = CBaseSimSnap::makeSlab(89906675+89915392+8718,10,10);
  
  CSimSnap *pSf= &gadget;
  
  fprintf(stderr,"main: Opening up FoF files.\n"); 

  CFoFFile fof(basepath,snapname,snapshot);

  // fprintf(stderr,"main: converting to physical units.\n");

  pSf->convertUnits(units::CUnit(units::len_kpc),units::CUnit(units::mass_Msol)*1.e10,units::CUnit(units::vel_kmPerS), units::CUnit(units::mass_Msol)*pow(units::CUnit(units::len_kpc),-3)*1.e10, pow(units::CUnit(units::vel_kmPerS),2) );

  // Now we're ready

  
  for(n=atoi(argv[4]);n<=atoi(argv[5]);n++) {
    float cx,cy,cz;
    auto_ptr<CSimSnap> currenthalo_nc = fof.getGroup(*pSf,n);
    CGeometry currenthalo(currenthalo_nc.get());
    
    fof.getGroupCentre(n,&cx,&cy,&cz);

    currenthalo.reCentre(cx,cy,cz);
    currenthalo.wrapping = true;
    currenthalo.apply();

    int nlines = currenthalo_nc->getNumParticles();

   
    CParticle *particles;

    
    sprintf(output,"%sextract_ap/%d.gadget",basepath,n);
    
    currenthalo_nc->write((string)output, CSimSnap::gadget);

    /*
    FILE *outputfile;
    if(!(outputfile = fopen(output,"w"))) {
      
      fprintf(stderr,"main: Failed opening output %s\n",output);
      exit(0);
    }
    
    printf("main: Now outputing %s (expect %d lines)\n",output,nlines);
   

    int pos;
    char writestring[2048];
   
    fprintf(outputfile,"(gadg snapshot: %s\nRedshift: %f   Boxsize (phys): %f   h=%f\n\nx\ty\tz\tvx\tvy\tvz\tmass\ttype\tTemp\n",basepath,(*pSf).getRedshift(),(*pSf).getBoxSize(),(*pSf).getHubble());

    const float BoxSize = (*pSf).getBoxSize();

    for(pos=0;pos<nlines;pos++) {

      CParticle *p = currenthalo.getParticle(pos);

      int type =0;
      switch(p->type) {
      case CParticle::gas:
	type = 0;
	break;
      case CParticle::dm:
	type = 1;
	break;
      case CParticle::star:
	type = 4;
	break;
      }

      fprintf(outputfile,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%f\n",p->x,p->y,p->z,p->vx,p->vy,p->vz,p->mass,type,p->temp);
      
      currenthalo.releaseParticle(p);

    }

    */

  }
  fprintf(stderr, "main: All done.\n");
}








  











