import netCDF4 as nc
import numpy as np
import os
import matplotlib.pyplot as plt

### class gen_netcdf ###{{{
class gen_netcdf:
    """Create netCDF files with attributes.

    Functions:
        create(ncname,ncdimlist,ncdimsize)
            Create a netCDF file and add dimensions.

        addglobatt(ncname,ncglobatt)
            Add global atributes (if they are needed) to the netCDF file.

        addvar(ncname,ncvarlist,ncvardata)
            Add avariable to the netCDF file, to add more than one variable use a loop (see example).

    Inputs:
        ncname: Name of netCDF file.
        ncdimlist: List with the name of dimensions variables.
        ncdimsize: List with the dimension length.
        ncvarlist: List of variables that can include the attributes.
        ncglobatt: List with global attributes.

    Example:
        ncname="example.nc"
        ncdimlist="y","x"
        ncdimsize=827,494

        ncvarlist=[['lon','f',['y','x'],'long_name','Longitude','units','degree west'],
                   ['lat','f',['y','x'],'units','degree north'],
                   ['mask','d',['y','x'],'units','mask']]

        Each row of the ncvarlist is composed by the variable name, format, dimensions and the attributes (attname,att,attname02,att02...) if they are neded.

        ncglobatt=[['title','Example'],['creator_name','M. Larranaga'],['creator_email','larranaga.marco@gmail.com']]

        create(ncname,ncdimlist,ncdimsize)

        addglobatt(ncname,ncglobatt)

        for j in range(0,len(ncvarlist)):
            addvar(ncname,ncvarlist[j],ncvarlist[j][0])
    """

    ### create ###{{{
    def create(ncname,ncdimlist,ncdimsize):
        """create(ncname,ncdimlist,ncdimsize)"""
    
        if os.path.exists(ncname):
            os.remove(ncname)
    
        ncfile=nc.Dataset(ncname,'w',format='NETCDF4')
    
        for j in range(0,len(ncdimlist)):
    
            ncfile.createDimension(ncdimlist[j],ncdimsize[j])
    
        ncfile.close()
    #}}}
    
    ### addglobatt ###{{{
    def addglobatt(ncname,ncglobatt):
        """addglobatt(ncname,ncglobatt)"""
    
        ncfile=nc.Dataset(ncname,'a')
    
        for j in range(0,len(ncglobatt)):
    
            ncfile.setncattr(ncglobatt[j][0],ncglobatt[j][1])
    
        ncfile.close()
    #}}}
    
    ### addvar ###{{{
    def addvar(ncname,ncvarlist,ncvardata):
        """addvar(ncname,ncvarlist,ncvardata)"""
    
        ncfile=nc.Dataset(ncname,'a')
    
        ncvar=ncfile.createVariable(ncvarlist[0],ncvarlist[1],ncvarlist[2])
    
        if len(ncvarlist) > 3:
            for j in range(3,len(ncvarlist),2):
                ncvar.setncattr(ncvarlist[j],ncvarlist[j+1])
    
        ncfile.variables[ncvarlist[0]]=ncvardata
    
        ncfile.close()
    #}}}
#}}}

### class nemo_tools ###{{{
class nemo_tools:

    ### geocurrents ####{{{
    def geocurrents(h,lon,lat):
        """u_g,v_g=geocurrents(h,lon,lat)

        Function that compute geostrophic current by considering a beta plane in latitudes near to 0 degrees.

        This function uses 2dx(y) to derivate.

        Inputs:
            h: Sea surface height or sea level anomaly (m). h could be a 2d or 3d matrix.
            lon: 2d longitude matrix (degrees).
            lat: 2d latitude matrix (degrees).

        Outputs"
            u_g: Nstward geostrophic currents (m/s)
            v_g: Northward geostrophic currents (m/s)
        """
    
        g=9.81
        omega=2*np.pi/(24*60*60)
        f=2*omega*np.sin(lat*np.pi/180)
    
        R=6371000;
        r=R*np.cos(lat*np.pi/180);
    
        deglat=111117.5*np.cos(lat*np.pi/180);
    
        if len(h.shape) == 2:
    
            ny,nx=h.shape
    
            dx=r[:,1:-1]*(2*np.pi)*(lon[:,2:]-lon[:,:-2])/360
            dy=lat[2:,:]*deglat[2:,:] - lat[:-2,:]*deglat[:-2,:]
    
            detay=(h[2:,:]-h[:-2,:])/dy
            detax=(h[:,2:]-h[:,:-2])/dx
        
            beta=(f[2:,:]-f[:-2,:])/dy*lat[1:-1,:]*deglat[1:-1,:]
            beta[np.where((beta < 2.3e-11) & (lat[1:-1,:] > 0))]=np.NaN
            beta[np.where((beta > -2.3e-11) & (lat[1:-1,:] < 0))]=np.NaN
        
            a1=(1 - np.exp(-(lat/2.2)**2))
            a1[np.where(lat == 0)]=0
        
            u_f=-g/f[1:-1,:]*detay*a1[1:-1,:]
            v_f=g/f[:,1:-1]*detax*a1[:,1:-1]
        
            b1=np.exp(-(lat/2.2)**2)
            b1[np.where(lat == 0)]=0
        
            u_b=-g/beta*(h[2:,:] + h[:-2,:] -2*h[1:-1,:])/dy**2*b1[1:-1,:]
            v_b=g/beta[:,1:-1]*(h[1:-1,2:] + h[1:-1,:-2] - 2*h[1:-1,1:-1])/dx[1:-1,:]**2*b1[1:-1,1:-1]
        
            u_g=u_f[:,1:-1] + u_b[:,1:-1]
            v_g=v_f[1:-1,:] + v_b
        
            u_f=u_f[:,1:-1]
            u_b=u_b[:,1:-1]
            v_f=v_f[1:-1,:]
    
        elif len(h.shape) > 2:
    
            ny,nx,nz=h.shape
    
            dx=np.ones((ny,nx-2,nz))
            dy=np.ones((ny-2,nx,nz))
        
            dx[:,:,0]=r[:,1:-1]*(2*np.pi)*(lon[:,2:]-lon[:,:-2])/360
            dy[:,:,0]=lat[2:,:]*deglat[2:,:] - lat[:-2,:]*deglat[:-2,:]
    
            for j in range(1,nz):
    
                dx[:,:,j]=dx[:,:,0]
                dy[:,:,j]=dy[:,:,0]
        
            detay=(h[2:,:,:]-h[:-2,:,:])/dy
            detax=(h[:,2:,:]-h[:,:-2,:])/dx
        
            beta=(f[2:,:,:]-f[:-2,:,:])/dy*lat[1:-1,:,:]*deglat[1:-1,:,:]
            beta[np.where((beta < 2.3e-11) & (lat[1:-1,:,:] > 0))]=np.NaN
            beta[np.where((beta > -2.3e-11) & (lat[1:-1,:,:] < 0))]=np.NaN
        
            a1=(1 - np.exp(-(lat/2.2)**2))
            a1[np.where(lat == 0)]=0
        
            u_f=-g/f[1:-1,:,:]*detay*a1[1:-1,:]
            v_f=g/f[:,1:-1,:]*detax*a1[:,1:-1]
        
            b1=np.exp(-(lat/2.2)**2)
            b1[np.where(lat == 0)]=0
        
            u_b=-g/beta*(h[2:,:,:] + h[:-2,:,:] -2*h[1:-1,:,:])/dy**2*b1[1:-1,:,:]
            v_b=g/beta[:,1:-1,:]*(h[1:-1,2:,:] + h[1:-1,:-2,:] - 2*h[1:-1,1:-1,:])/dx[1:-1,:,:]**2*b1[1:-1,1:-1,:]
        
            u_g=u_f[:,1:-1,:] + u_b[:,1:-1,:]
            v_g=v_f[1:-1,:,:] + v_b
        
            u_f=u_f[:,1:-1,:]
            u_b=u_b[:,1:-1,:]
            v_f=v_f[1:-1,:,:]

        return [u_g,v_g];
    #}}}

#}}}

### class quick_plot ###{{{
class quick_plot:
    """This module allows you to do quick plots that can be configured by basics parameters.

    """
    ### pcolor ###{{{
    def pcolor(*args):
        """This function allows you to do quick pcolors that includes a colorbar. The functions can be used with or without including "x" and "y" coordinates. Quick pcolor can be setup by the following parameters:

            cmap: String that define a colormap to be used in the pcolor. Default value 'viridis'.
            vmin: Used to define the colormap range.
            vmax: Used to define the colormap range.
            axis: Switch between axis='auto' and axis='equal'
            xlim: Used to define the limits of the horizontal axis.
            ylim: Used to define the limits of the vercical axis.

            Inputs:
                The first argument must be composed by the variables to plot.

                    pcolor(var)
                    pcolor([x,y,var])

                var: Variable to plot
                x,y: Coordinates eg. longitude an latitude.

                If optional arguments are needed, those have to be included as lists:

                    pcolor(var,['cmap','seismic'],['axis','equal'])

                    or

                    pcolor([x,y,var,['cmap','seismic'],['axis','equal'])
        """

        param=['cmap','vmin','vmax','axis','xlim','ylim']

        if type(args[0]) != list:
            var=args[0]
            ny,nx=var.shape

            paramv=['viridis',np.min(var),np.max(var),'auto',[0,nx],[0,ny]]

        elif len(args[0]) == 3:
            x=args[0][0]
            y=args[0][1]
            var=args[0][2]

            paramv=['viridis',np.min(var),np.max(var),'auto',[np.min(x),np.max(x)],[np.min(y),np.max(y)]] 


        if len(args) > 1:

            paramin=[0]*(len(args)-1)
            paramvin=[0]*(len(args)-1)
            for j in range(1,len(args)):

                vargs=args[j]

                paramin[j-1]=vargs[0]
                paramvin[j-1]=vargs[1]


            for k in range(0,len(param)):

                ind=[paramin.index(i) for i in paramin if param[k] in i]

                if len(ind) == 1:

                    paramv[k]=paramvin[ind[0]]

        fig,ax=plt.subplots(1,1)
        if type(args[0]) != list:
            pclr=ax.pcolor(var,cmap=paramv[0],vmin=paramv[1],vmax=paramv[2])
        elif len(args[0]) == 3:
            pclr=ax.pcolor(x,y,var,cmap=paramv[0],vmin=paramv[1],vmax=paramv[2])

        fig.colorbar(pclr,ax=ax)

        ax.axis(paramv[3])
        ax.set_xlim(paramv[4])
        ax.set_ylim(paramv[5])

        fig.show()
    #}}}
#}}}

