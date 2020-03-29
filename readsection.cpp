#include "readsection.h"
// im must be free by the user.  axis=[0,1,2]->[x,y,z] , slice is the distance from the volume center.
bool readMRCSection(float **im, int &dimx, int &dimy, const char* filename, int axis, int slice)
{
    MRC mrc;
    if(mrc.open(filename, "rb")<=0)
    {
        printf("Warning: Failed to read %s\n",filename);
        return false;
    }
    int nx=mrc.getNx();
    int ny=mrc.getNy();
    int nz=mrc.getNz();

    float *buf=NULL;
    int ns,i,j;
    switch(axis)
    {
        case 0: //along x
            dimx=ny;
            dimy=nz;
            //ns=slice+nx/2;
            ns = slice;
            if(ns<0) ns=0;
            if(ns>=nx) ns=nx-1;
            buf=new float[dimx*dimy];
            for(j=0;j<dimy;j++)
                for(i=0;i<dimx;i++)
                {
                    mrc.readPixel(&buf[i+j*dimx], j, i, ns);
                }
            break;

        case 1: //along y
            dimx=nx;
            dimy=nz;
            //ns=slice+ny/2;
            ns=slice;
            if(ns<0) ns=0;
            if(ns>=nx) ns=nx-1;
            buf=new float[dimx*dimy];
            for(j=0;j<dimy;j++)
                for(i=0;i<dimx;i++)
                {
                    mrc.readPixel(&buf[i+j*dimx], j, ns, i);
                }
            break;

        case 2: //along z
            dimx=nx;
            dimy=ny;
            //ns=slice+nz/2;
            ns=slice;
            if(ns<0) ns=0;
            if(ns>=nx) ns=nx-1;
            buf=new float[dimx*dimy];
            mrc.read2DIm_32bit(buf,ns);
            break;
    }
    mrc.close();

    *im=buf;

    return true;
}

// im must be free by the user.  axis=[0,1,2]->[x,y,z] , slice starts from zero along the axis direction.
bool readMRCSection2(float **im, int &dimx, int &dimy, MRC &mrc, int axis, int slice)
{
    int nx=mrc.getNx();
    int ny=mrc.getNy();
    int nz=mrc.getNz();

    float *buf=NULL;
    int ns,i,j;
    switch(axis)
    {
        case 0: //along x
            dimx=ny;
            dimy=nz;
            //ns=slice+nx/2;
            ns = slice;
            if(ns<0) ns=0;
            if(ns>=nx) ns=nx-1;
            buf=new float[dimx*dimy];
            for(j=0;j<dimy;j++)
                for(i=0;i<dimx;i++)
                {
                    mrc.readPixel(&buf[i+j*dimx], j, i, ns);
                }
            break;

        case 1: //along y
            dimx=nx;
            dimy=nz;
            //ns=slice+ny/2;
            ns=slice;
            if(ns<0) ns=0;
            if(ns>=nx) ns=nx-1;
            buf=new float[dimx*dimy];
            for(j=0;j<dimy;j++)
                for(i=0;i<dimx;i++)
                {
                    mrc.readPixel(&buf[i+j*dimx], j, ns, i);
                }
            break;

        case 2: //along z
            dimx=nx;
            dimy=ny;
            //ns=slice+nz/2;
            ns=slice;
            if(ns<0) ns=0;
            if(ns>=nx) ns=nx-1;
            buf=new float[dimx*dimy];
            mrc.read2DIm_32bit(buf,ns);
            break;
    }

    *im=buf;

    return true;
}

