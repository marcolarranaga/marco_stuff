import netCDF4 as nc
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import path
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import cartopy
import cartopy.crs as ccrs
from datetime import datetime
import sys

### marco_stuff ###{{{
class marco_stuff:
    """marco_stuff

    Functions:
        data=load(fname)
        whosdic(dic)
        mask=inpolygon(xpoly,ypoly,xgrid,ygrid)
    """

    ### load ### {{{
    def load(fname):
        dic={}
        data=np.load(fname)
        varnames=data['varnames']
        for j in range(0,len(varnames)):
            dic[varnames[j]]=data[varnames[j]]

        return dic
    #}}}

    ### whosdic ###{{{
    def whos(dic):
        varnames=[]
        info=[]
        varnameslen=[]
        vartype=[]
        vartypelen=[]
        sizename=[]
        for dicvarnames,dicinfo in dic.items():
            varnames.append(dicvarnames)
            varnameslen.append(len(dicvarnames))
            vartype.append(type(dic[dicvarnames]).__name__)
            vartypelen.append(len(type(dic[dicvarnames]).__name__))
            if type(dic[dicvarnames]).__name__ == 'ndarray':
                info.append(dic[dicvarnames].shape)
                sizename.append('shape:')
            elif type(dic[dicvarnames]).__name__ == 'list':
                info.append(len(dic[dicvarnames]))
                sizename.append('  len:')
            elif type(dic[dicvarnames]).__name__ == 'tuple':
                info.append(len(dic[dicvarnames]))
                sizename.append('  len:')
            else:
                info.append(' ')

        if max(varnameslen)+4 > 13:
            titstr1='Variables'.ljust(max(varnameslen)+4)
        else:
            titstr1='Variables'.ljust(13)
 
        if max(vartypelen)+4 > 11:
            titstr2='Type'.ljust(max(vartypelen)+4)
        else:
            titstr2='Type'.ljust(11)

        tit=titstr1+titstr2+'Data/info'
 
        print(tit)
        print('-'.ljust(len(tit),'-'))

        for j in range(0,len(varnames)):
            print(varnames[j].ljust(len(titstr1)-1),vartype[j].ljust(len(titstr2)-1),sizename[j],info[j])
    #}}}

    ### inpolygon ###{{{
    def inpolygon(xpoly,ypoly,xgrid,ygrid):

        ny,nx=xgrid.shape

        poly=[]
        for j in range(0,len(xpoly)): 
            poly.append((xpoly[j],ypoly[j])) 

        p=path.Path(poly)

        xgrid=xgrid.flatten()
        ygrid=ygrid.flatten()
        points=np.vstack((xgrid,ygrid)).T

        mask=p.contains_points(points)
        mask=mask.reshape((ny,nx))

        return mask
    #}}}

#}}}

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

    ### vort ###{{{
    def vort(lon,lat,u,v,**kwargs):

        paramin=[]
        paramvin=[]
        for key, value in kwargs.items():

            paramin.append(key)
            paramvin.append(value)

        g=9.81
        omega=2*np.pi/(24*60*60)
        f=2*omega*np.sin(lat*np.pi/180)
    
        R=6371000;
        r=R*np.cos(lat*np.pi/180);
    
        deglat=111117.5*np.cos(lat*np.pi/180);

        if paramvin[0] == 'simple':

            dx=r[:,1:]*(2*np.pi)*(lon[:,1:]-lon[:,:-1])/360
            dy=lat[1:,:]*deglat[1:,:] - lat[:-1,:]*deglat[:-1,:]

            dvdx=(v[:,1:] - v[:,:-1])/dx
            dudy=(u[1:,:] - u[:-1,:])/dy

            vort=np.ones((ny,nx))*np.NaN

            vort[0:-1,0:-1]=dvdx[0:-1,:] - dudy[:,0:-1]

        elif paramvin[0] == 'double':

            ny,nx=lon.shape

            dx=r[:,1:-1]*(2*np.pi)*(lon[:,2:]-lon[:,:-2])/360
            dy=lat[2:,:]*deglat[2:,:] - lat[:-2,:]*deglat[:-2,:]

            dvdx=(v[:,2:] - v[:,:-2])/dx
            dudy=(u[2:,:] - u[:-2,:])/dy

            vort=np.ones((ny,nx))*np.NaN

            vort[1:-1,1:-1]=dvdx[1:-1,:] - dudy[:,1:-1]

        return vort;
            #}}}

#}}}

### class quick_plot ###{{{
class quick_plot:
    """This module allows you to do quick plots that can be configured by basics parameters.

    """
    ### pcolor ###{{{
    def pcolor(args,**kwargs):
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

        if type(args) != list:
            var=args[:]*1
            ny,nx=var.shape

            paramv=['viridis',np.nanmin(var),np.nanmax(var),'auto',[0,nx],[0,ny]]

        elif len(args) == 3:
            x=args[0]
            y=args[1]
            var=args[2]

            paramv=['viridis',np.nanmin(var),np.nanmax(var),'auto',[np.min(x),np.max(x)],[np.min(y),np.max(y)]] 


#        if len(kwargs) > 0:
        paramin=[]
        paramvin=[]
        for key, value in kwargs.items():

            paramin.append(key)
            paramvin.append(value)


        if len(paramin) > 0:
            for k in range(0,len(param)):

                ind=[paramin.index(i) for i in paramin if param[k] in i]

                if len(ind) == 1:

                    paramv[k]=paramvin[ind[0]]

        fig,ax=plt.subplots(1,1)
        if type(args) != list:
            pclr=ax.pcolormesh(var,cmap=paramv[0],vmin=paramv[1],vmax=paramv[2])
        elif len(args) == 3:
            pclr=ax.pcolormesh(x,y,var,cmap=paramv[0],vmin=paramv[1],vmax=paramv[2])

        fig.colorbar(pclr,ax=ax)

        ax.axis(paramv[3])
        ax.set_xlim(paramv[4])
        ax.set_ylim(paramv[5])

        ind=[paramin.index(i) for i in paramin if 'title' in i]

        if len(kwargs) > 0:
            if len(ind) == 1:
                ax.set_title(paramvin[ind[0]])

        fig.show()
    #}}}

    ### fill_poly ###{{{
    def fill_poly(args,**kwargs):

        param=['ax'      ,'alpha','edgecolor','facecolor','label','linestyle','linewidth','zorder','transform']
        paramv=[plt.gca(),1      ,'k'        ,'b'        ,''     ,'-'        ,1          ,0       ,ccrs.PlateCarree()]
        x=args[0]*1
        y=args[1]*1

        paramin=[]
        paramvin=[]
        for key, value in kwargs.items():

            paramin.append(key)
            paramvin.append(value)


        if len(paramin) > 0:
            for k in range(0,len(param)):

                ind=[paramin.index(i) for i in paramin if param[k] in i]

                if len(ind) == 1:

                    paramv[k]=paramvin[ind[0]]


        ind=np.where(np.isnan(x) == 1)[0]

        for j in range(len(ind)-1):

            polygon=np.zeros((len(range(ind[j]+1,ind[j+1]-2)),2))
            polygon[:,0]=x[ind[j]+1:ind[j+1]-2]
            polygon[:,1]=y[ind[j]+1:ind[j+1]-2]

            poly=mpatches.Polygon(polygon,alpha=paramv[0],edgecolor=paramv[1],facecolor=paramv[2],
                    label=paramv[3],linestyle=paramv[4],linewidth=paramv[5],zorder=paramv[6],
                    transform=paramv[7])

            ax.add_patch(poly)

        return poly

#}}}

###}}}

### class colormap ###{{{
class colormap:

    ### parula ###{{{
    def parula(*args):

        # colormap #{{{
        parulaorig=[[0.2422,0.1504,0.6603],
                    [0.2504,0.1650,0.7076],
                    [0.2578,0.1818,0.7511],
                    [0.2647,0.1978,0.7952],
                    [0.2706,0.2147,0.8364],
                    [0.2751,0.2342,0.8710],
                    [0.2783,0.2559,0.8991],
                    [0.2803,0.2782,0.9221],
                    [0.2813,0.3006,0.9414],
                    [0.2810,0.3228,0.9579],
                    [0.2795,0.3447,0.9717],
                    [0.2760,0.3667,0.9829],
                    [0.2699,0.3892,0.9906],
                    [0.2602,0.4123,0.9952],
                    [0.2440,0.4358,0.9988],
                    [0.2206,0.4603,0.9973],
                    [0.1963,0.4847,0.9892],
                    [0.1834,0.5074,0.9798],
                    [0.1786,0.5289,0.9682],
                    [0.1764,0.5499,0.9520],
                    [0.1687,0.5703,0.9359],
                    [0.1540,0.5902,0.9218],
                    [0.1460,0.6091,0.9079],
                    [0.1380,0.6276,0.8973],
                    [0.1248,0.6459,0.8883],
                    [0.1113,0.6635,0.8763],
                    [0.0952,0.6798,0.8598],
                    [0.0689,0.6948,0.8394],
                    [0.0297,0.7082,0.8163],
                    [0.0036,0.7203,0.7917],
                    [0.0067,0.7312,0.7660],
                    [0.0433,0.7411,0.7394],
                    [0.0964,0.7500,0.7120],
                    [0.1408,0.7584,0.6842],
                    [0.1717,0.7670,0.6554],
                    [0.1938,0.7758,0.6251],
                    [0.2161,0.7843,0.5923],
                    [0.2470,0.7918,0.5567],
                    [0.2906,0.7973,0.5188],
                    [0.3406,0.8008,0.4789],
                    [0.3909,0.8029,0.4354],
                    [0.4456,0.8024,0.3909],
                    [0.5044,0.7993,0.3480],
                    [0.5616,0.7942,0.3045],
                    [0.6174,0.7876,0.2612],
                    [0.6720,0.7793,0.2227],
                    [0.7242,0.7698,0.1910],
                    [0.7738,0.7598,0.1646],
                    [0.8203,0.7498,0.1535],
                    [0.8634,0.7406,0.1596],
                    [0.9035,0.7330,0.1774],
                    [0.9393,0.7288,0.2100],
                    [0.9728,0.7298,0.2394],
                    [0.9956,0.7434,0.2371],
                    [0.9970,0.7659,0.2199],
                    [0.9952,0.7893,0.2028],
                    [0.9892,0.8136,0.1885],
                    [0.9786,0.8386,0.1766],
                    [0.9676,0.8639,0.1643],
                    [0.9610,0.8890,0.1537],
                    [0.9597,0.9135,0.1423],
                    [0.9628,0.9373,0.1265],
                    [0.9691,0.9606,0.1064],
                    [0.9769,0.9839,0.0805]]
                     #}}}
        parulaorig=np.array(parulaorig)
        rows=np.linspace(0,1,64)

        if len(args) == 0:
            parula=ListedColormap(parulaorig)

        elif len(args) == 1:

            N=args[0]
            rowsi=np.linspace(0,1,N)
            parula=np.ones((len(rowsi),3))*np.NaN

            parula[:,0]=np.interp(rowsi,rows,parulaorig[:,0])
            parula[:,1]=np.interp(rowsi,rows,parulaorig[:,1])
            parula[:,2]=np.interp(rowsi,rows,parulaorig[:,2])

            parula=ListedColormap(parula)

        return parula
    ###}}}

    ### parula_r ###{{{
    def parula_r(*args):

        # colormap #{{{
        parulaorig=[[0.2422,0.1504,0.6603],
                    [0.2504,0.1650,0.7076],
                    [0.2578,0.1818,0.7511],
                    [0.2647,0.1978,0.7952],
                    [0.2706,0.2147,0.8364],
                    [0.2751,0.2342,0.8710],
                    [0.2783,0.2559,0.8991],
                    [0.2803,0.2782,0.9221],
                    [0.2813,0.3006,0.9414],
                    [0.2810,0.3228,0.9579],
                    [0.2795,0.3447,0.9717],
                    [0.2760,0.3667,0.9829],
                    [0.2699,0.3892,0.9906],
                    [0.2602,0.4123,0.9952],
                    [0.2440,0.4358,0.9988],
                    [0.2206,0.4603,0.9973],
                    [0.1963,0.4847,0.9892],
                    [0.1834,0.5074,0.9798],
                    [0.1786,0.5289,0.9682],
                    [0.1764,0.5499,0.9520],
                    [0.1687,0.5703,0.9359],
                    [0.1540,0.5902,0.9218],
                    [0.1460,0.6091,0.9079],
                    [0.1380,0.6276,0.8973],
                    [0.1248,0.6459,0.8883],
                    [0.1113,0.6635,0.8763],
                    [0.0952,0.6798,0.8598],
                    [0.0689,0.6948,0.8394],
                    [0.0297,0.7082,0.8163],
                    [0.0036,0.7203,0.7917],
                    [0.0067,0.7312,0.7660],
                    [0.0433,0.7411,0.7394],
                    [0.0964,0.7500,0.7120],
                    [0.1408,0.7584,0.6842],
                    [0.1717,0.7670,0.6554],
                    [0.1938,0.7758,0.6251],
                    [0.2161,0.7843,0.5923],
                    [0.2470,0.7918,0.5567],
                    [0.2906,0.7973,0.5188],
                    [0.3406,0.8008,0.4789],
                    [0.3909,0.8029,0.4354],
                    [0.4456,0.8024,0.3909],
                    [0.5044,0.7993,0.3480],
                    [0.5616,0.7942,0.3045],
                    [0.6174,0.7876,0.2612],
                    [0.6720,0.7793,0.2227],
                    [0.7242,0.7698,0.1910],
                    [0.7738,0.7598,0.1646],
                    [0.8203,0.7498,0.1535],
                    [0.8634,0.7406,0.1596],
                    [0.9035,0.7330,0.1774],
                    [0.9393,0.7288,0.2100],
                    [0.9728,0.7298,0.2394],
                    [0.9956,0.7434,0.2371],
                    [0.9970,0.7659,0.2199],
                    [0.9952,0.7893,0.2028],
                    [0.9892,0.8136,0.1885],
                    [0.9786,0.8386,0.1766],
                    [0.9676,0.8639,0.1643],
                    [0.9610,0.8890,0.1537],
                    [0.9597,0.9135,0.1423],
                    [0.9628,0.9373,0.1265],
                    [0.9691,0.9606,0.1064],
                    [0.9769,0.9839,0.0805]]
                     #}}}
        parulaorig=np.array(parulaorig[::-1])
        rows=np.linspace(0,1,64)

        if len(args) == 0:
            parula_r=ListedColormap(parulaorig)

        elif len(args) == 1:

            N=args[0]
            rowsi=np.linspace(0,1,N)
            parula_r=np.ones((len(rowsi),3))*np.NaN

            parula_r[:,0]=np.interp(rowsi,rows,parulaorig[:,0])
            parula_r[:,1]=np.interp(rowsi,rows,parulaorig[:,1])
            parula_r[:,2]=np.interp(rowsi,rows,parulaorig[:,2])

            parula_r=ListedColormap(parula_r)

        return parula_r
    ###}}}

    ### parulaw ###{{{
    def parulaw(*args):

        # colormap #{{{
        parulaworig=[[0.9999,0.9999,0.9999],
                     [0.9460,0.9606,0.9999],
                     [0.8919,0.9212,0.9999],
                     [0.8380,0.8819,0.9998],
                     [0.7840,0.8426,0.9998],
                     [0.7301,0.8032,0.9998],
                     [0.6761,0.7639,0.9998],
                     [0.6221,0.7245,0.9998],
                     [0.5681,0.6851,0.9997],
                     [0.5142,0.6458,0.9997],
                     [0.4602,0.6065,0.9997],
                     [0.4062,0.5672,0.9997],
                     [0.3523,0.5278,0.9996],
                     [0.2983,0.4885,0.9996],
                     [0.2443,0.4491,0.9996],
                     [0.2160,0.4650,0.9958],
                     [0.1926,0.4890,0.9877],
                     [0.1824,0.5113,0.9778],
                     [0.1781,0.5326,0.9656],
                     [0.1759,0.5534,0.9490],
                     [0.1663,0.5737,0.9335],
                     [0.1522,0.5934,0.9194],
                     [0.1451,0.6121,0.9059],
                     [0.1362,0.6305,0.8959],
                     [0.1227,0.6487,0.8868],
                     [0.1092,0.6660,0.8741],
                     [0.0923,0.6821,0.8570],
                     [0.0640,0.6967,0.8363],
                     [0.0244,0.7099,0.8131],
                     [0.0021,0.7218,0.7883],
                     [0.0092,0.7326,0.7626],
                     [0.0502,0.7423,0.7360],
                     [0.1025,0.7510,0.7087],
                     [0.1452,0.7594,0.6809],
                     [0.1746,0.7680,0.6521],
                     [0.1959,0.7767,0.6216],
                     [0.2187,0.7852,0.5887],
                     [0.2508,0.7925,0.5530],
                     [0.2954,0.7977,0.5150],
                     [0.3454,0.8011,0.4749],
                     [0.3956,0.8030,0.4314],
                     [0.4507,0.8022,0.3872],
                     [0.5092,0.7989,0.3445],
                     [0.5660,0.7938,0.3010],
                     [0.6215,0.7870,0.2582],
                     [0.6758,0.7786,0.2203],
                     [0.7276,0.7692,0.1891],
                     [0.7769,0.7592,0.1633],
                     [0.8229,0.7492,0.1536],
                     [0.8657,0.7401,0.1603],
                     [0.9055,0.7327,0.1788],
                     [0.9408,0.7287,0.2117],
                     [0.9742,0.7300,0.2402],
                     [0.9960,0.7442,0.2365],
                     [0.9970,0.7667,0.2193],
                     [0.9951,0.7900,0.2022],
                     [0.9890,0.8142,0.1882],
                     [0.9784,0.8392,0.1763],
                     [0.9675,0.8644,0.1641],
                     [0.9609,0.8894,0.1535],
                     [0.9597,0.9137,0.1421],
                     [0.9628,0.9375,0.1264],
                     [0.9691,0.9607,0.1063],
                     [0.9769,0.9839,0.0805]]
                     #}}}
        parulaworig=np.array(parulaworig)
        rows=np.linspace(0,1,64)

        if len(args) == 0:
            parulaw=ListedColormap(parulaworig)

        elif len(args) == 1:

            N=args[0]
            rowsi=np.linspace(0,1,N)
            parulaw=np.ones((len(rowsi),3))*np.NaN

            parulaw[:,0]=np.interp(rowsi,rows,parulaworig[:,0])
            parulaw[:,1]=np.interp(rowsi,rows,parulaworig[:,1])
            parulaw[:,2]=np.interp(rowsi,rows,parulaworig[:,2])

            parulaw=ListedColormap(parulaw)

        return parulaw
    ###}}}

    ### parulaw_r ###{{{
    def parulaw_r(*args):

        # colormap #{{{
        parulaworig=[[0.9999,0.9999,0.9999],
                     [0.9460,0.9606,0.9999],
                     [0.8919,0.9212,0.9999],
                     [0.8380,0.8819,0.9998],
                     [0.7840,0.8426,0.9998],
                     [0.7301,0.8032,0.9998],
                     [0.6761,0.7639,0.9998],
                     [0.6221,0.7245,0.9998],
                     [0.5681,0.6851,0.9997],
                     [0.5142,0.6458,0.9997],
                     [0.4602,0.6065,0.9997],
                     [0.4062,0.5672,0.9997],
                     [0.3523,0.5278,0.9996],
                     [0.2983,0.4885,0.9996],
                     [0.2443,0.4491,0.9996],
                     [0.2160,0.4650,0.9958],
                     [0.1926,0.4890,0.9877],
                     [0.1824,0.5113,0.9778],
                     [0.1781,0.5326,0.9656],
                     [0.1759,0.5534,0.9490],
                     [0.1663,0.5737,0.9335],
                     [0.1522,0.5934,0.9194],
                     [0.1451,0.6121,0.9059],
                     [0.1362,0.6305,0.8959],
                     [0.1227,0.6487,0.8868],
                     [0.1092,0.6660,0.8741],
                     [0.0923,0.6821,0.8570],
                     [0.0640,0.6967,0.8363],
                     [0.0244,0.7099,0.8131],
                     [0.0021,0.7218,0.7883],
                     [0.0092,0.7326,0.7626],
                     [0.0502,0.7423,0.7360],
                     [0.1025,0.7510,0.7087],
                     [0.1452,0.7594,0.6809],
                     [0.1746,0.7680,0.6521],
                     [0.1959,0.7767,0.6216],
                     [0.2187,0.7852,0.5887],
                     [0.2508,0.7925,0.5530],
                     [0.2954,0.7977,0.5150],
                     [0.3454,0.8011,0.4749],
                     [0.3956,0.8030,0.4314],
                     [0.4507,0.8022,0.3872],
                     [0.5092,0.7989,0.3445],
                     [0.5660,0.7938,0.3010],
                     [0.6215,0.7870,0.2582],
                     [0.6758,0.7786,0.2203],
                     [0.7276,0.7692,0.1891],
                     [0.7769,0.7592,0.1633],
                     [0.8229,0.7492,0.1536],
                     [0.8657,0.7401,0.1603],
                     [0.9055,0.7327,0.1788],
                     [0.9408,0.7287,0.2117],
                     [0.9742,0.7300,0.2402],
                     [0.9960,0.7442,0.2365],
                     [0.9970,0.7667,0.2193],
                     [0.9951,0.7900,0.2022],
                     [0.9890,0.8142,0.1882],
                     [0.9784,0.8392,0.1763],
                     [0.9675,0.8644,0.1641],
                     [0.9609,0.8894,0.1535],
                     [0.9597,0.9137,0.1421],
                     [0.9628,0.9375,0.1264],
                     [0.9691,0.9607,0.1063],
                     [0.9769,0.9839,0.0805]]
                     #}}}
        parulaworig=np.array(parulaworig[::-1])
        rows=np.linspace(0,1,64)

        if len(args) == 0:
            parulaw_r=ListedColormap(parulaworig)

        elif len(args) == 1:

            N=args[0]
            rowsi=np.linspace(0,1,N)
            parulaw_r=np.ones((len(rowsi),3))*np.NaN

            parulaw_r[:,0]=np.interp(rowsi,rows,parulaworig[:,0])
            parulaw_r[:,1]=np.interp(rowsi,rows,parulaworig[:,1])
            parulaw_r[:,2]=np.interp(rowsi,rows,parulaworig[:,2])

            parulaw_r=ListedColormap(parulaw_r)

        return parulaw_r
    ###}}}

    ### jet ###{{{
    def jet(*args):

        # colormap #{{{
        jetorig=[[     0,     0,0.5625], 
                 [     0,     0,0.6250],
                 [     0,     0,0.6875],
                 [     0,     0,0.7500],
                 [     0,     0,0.8125],
                 [     0,     0,0.8750],
                 [     0,     0,0.9375],
                 [     0,     0,1.0000],
                 [     0,0.0625,1.0000],
                 [     0,0.1250,1.0000],
                 [     0,0.1875,1.0000],
                 [     0,0.2500,1.0000],
                 [     0,0.3125,1.0000],
                 [     0,0.3750,1.0000],
                 [     0,0.4375,1.0000],
                 [     0,0.5000,1.0000],
                 [     0,0.5625,1.0000],
                 [     0,0.6250,1.0000],
                 [     0,0.6875,1.0000],
                 [     0,0.7500,1.0000],
                 [     0,0.8125,1.0000],
                 [     0,0.8750,1.0000],
                 [     0,0.9375,1.0000],
                 [     0,1.0000,1.0000],
                 [0.0625,1.0000,0.9375],
                 [0.1250,1.0000,0.8750],
                 [0.1875,1.0000,0.8125],
                 [0.2500,1.0000,0.7500],
                 [0.3125,1.0000,0.6875],
                 [0.3750,1.0000,0.6250],
                 [0.4375,1.0000,0.5625],
                 [0.5000,1.0000,0.5000],
                 [0.5625,1.0000,0.4375],
                 [0.6250,1.0000,0.3750],
                 [0.6875,1.0000,0.3125],
                 [0.7500,1.0000,0.2500],
                 [0.8125,1.0000,0.1875],
                 [0.8750,1.0000,0.1250],
                 [0.9375,1.0000,0.0625],
                 [1.0000,1.0000,     0],
                 [1.0000,0.9375,     0],
                 [1.0000,0.8750,     0],
                 [1.0000,0.8125,     0],
                 [1.0000,0.7500,     0],
                 [1.0000,0.6875,     0],
                 [1.0000,0.6250,     0],
                 [1.0000,0.5625,     0],
                 [1.0000,0.5000,     0],
                 [1.0000,0.4375,     0],
                 [1.0000,0.3750,     0],
                 [1.0000,0.3125,     0],
                 [1.0000,0.2500,     0],
                 [1.0000,0.1875,     0],
                 [1.0000,0.1250,     0],
                 [1.0000,0.0625,     0],
                 [1.0000,     0,     0],
                 [0.9375,     0,     0],
                 [0.8750,     0,     0],
                 [0.8125,     0,     0],
                 [0.7500,     0,     0],
                 [0.6875,     0,     0],
                 [0.6250,     0,     0],
                 [0.5625,     0,     0],
                 [0.5000,     0,     0],]
                     #}}}
        jetorig=np.array(jetorig)
        rows=np.linspace(0,1,64)

        if len(args) == 0:
            jet=ListedColormap(jetorig)

        elif len(args) == 1:

            N=args[0]
            rowsi=np.linspace(0,1,N)
            jet=np.ones((len(rowsi),3))*np.NaN

            jet[:,0]=np.interp(rowsi,rows,jetorig[:,0])
            jet[:,1]=np.interp(rowsi,rows,jetorig[:,1])
            jet[:,2]=np.interp(rowsi,rows,jetorig[:,2])

            jet=ListedColormap(jet)

        return jet
    ###}}}

    ### redbluemod ###{{{
    def redbluemod(*args):

        # colormap #{{{

        # Config{{{
        bd=np.array([27 ,38 ,49 ])/255
        b =np.array([0  ,0  ,200])/255
        g =np.array([26 ,82 ,118])/255
        g =np.array([52 ,152,219])/255
        gl=np.array([171,235,198])/255
        w =np.array([255,255,255])/255
        yl=np.array([249,231,159])/255
        y =np.array([230,126,34 ])/255
        r =np.array([200,0  ,0  ])/255
        rd=np.array([100,30 ,22 ])/255

        np_bd_b=40
        np_b_g =40
        np_g_gl=35
        np_gl_w=5
        np_w_yl=5
        np_yl_y=35
        np_y_r =40
        np_r_rd=40
        np_sumi=0
        #}}}

        # bd to b {{{
        R=np.interp(np.array(range(np_bd_b+1)),np.array((0,np_bd_b)),np.array((bd[0],b[0])))
        G=np.interp(np.array(range(np_bd_b+1)),np.array((0,np_bd_b)),np.array((bd[1],b[1])))
        B=np.interp(np.array(range(np_bd_b+1)),np.array((0,np_bd_b)),np.array((bd[2],b[2])))
        np_sumi=np_sumi + np_bd_b
        np_sumf=np_sumi + np_b_g
        #}}}
        
        # b to g {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((b[0],g[0])))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((b[1],g[1])))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((b[2],g[2])))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_b_g
        np_sumf=np_sumi + np_g_gl
        #}}}

        # g to gl {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((g[0],gl[0])))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((g[1],gl[1])))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((g[2],gl[2])))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_g_gl
        np_sumf=np_sumi + np_gl_w
        #}}}

        # gl to w {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((gl[0],1)))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((gl[1],1)))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((gl[2],1)))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_gl_w
        np_sumf=np_sumi + np_w_yl
        #}}}

        # w to yl {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((1,yl[0])))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((1,yl[1])))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((1,yl[2])))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_w_yl
        np_sumf=np_sumi + np_yl_y
        #}}}

        # yl to y {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((yl[0],y[0])))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((yl[1],y[1])))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((yl[2],y[2])))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_yl_y
        np_sumf=np_sumi + np_y_r
        #}}}

        # y to r {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((y[0],r[0])))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((y[1],r[1])))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((y[2],r[2])))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_y_r
        np_sumf=np_sumi + np_r_rd
        #}}}

        # r to rd {{{
        Rtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((r[0],rd[0])))
        Gtmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((r[1],rd[1])))
        Btmp=np.interp(np.array(range(np_sumi,np_sumf+1)),np.array((np_sumi,np_sumf)),
            np.array((r[2],rd[2])))

        R=np.append(R,Rtmp[1:])
        G=np.append(G,Gtmp[1:])
        B=np.append(B,Btmp[1:])
        np_sumi=np_sumi + np_r_rd
        #}}}

        #}}}

        redbluemodorig=np.array((R,G,B)).transpose()
        rows=np.linspace(0,1,np.max(redbluemodorig.shape))

        if len(args) == 0:

            N=65
            rowsi=np.linspace(0,1,N)
            redbluemod=np.ones((len(rowsi),3))*np.NaN
            redbluemod[:,0]=np.interp(rowsi,rows,redbluemodorig[:,0])
            redbluemod[:,1]=np.interp(rowsi,rows,redbluemodorig[:,1])
            redbluemod[:,2]=np.interp(rowsi,rows,redbluemodorig[:,2])

            redbluemod=ListedColormap(redbluemod)

        elif len(args) == 1:

            N=args[0]
            rowsi=np.linspace(0,1,N)
            redbluemod=np.ones((len(rowsi),3))*np.NaN

            redbluemod[:,0]=np.interp(rowsi,rows,redbluemodorig[:,0])
            redbluemod[:,1]=np.interp(rowsi,rows,redbluemodorig[:,1])
            redbluemod[:,2]=np.interp(rowsi,rows,redbluemodorig[:,2])

            redbluemod=ListedColormap(redbluemod)

        return redbluemod
    ###}}}


#}}}

### class remtime ###{{{
class remtime:

    def remtime1(strend):

        itcont=1;
        deltat=0;
        ti=nc.date2num(datetime.now(),'days since 1-1-1')

        str=''
        for l in range(0,20):
            str=f'{str}.'

        str=f'{str} 0% - --:--:-- - {strend}'

        print(str)

        return [ti,deltat,itcont,strend]

    def remtime2(totalt,itcont,deltat,ti,strend):

        tf=nc.date2num(datetime.now(),'days since 1-1-1')

        deltat=deltat + (tf-ti)

        remt=deltat/itcont*(totalt - itcont)

        progprcnt=int(np.floor(100 - (totalt - itcont)/totalt*100))
        progprcnt2=int(np.floor(20 - (totalt - itcont)/totalt*20))

        sys.stdout.write("\033[F") #back to previous line
        sys.stdout.write("\033[K") #clear line

        str=''
        for l in range(0,progprcnt2):
            str=f'{str}#'

        for l in range(0,20-progprcnt2):
            str=f'{str}.'

        remHH=int(np.floor(remt*24))
        remMM=int(np.floor((remt*24 - remHH)*60))
        remSS=int(np.floor(((remt*24 - remHH)*60 -  remMM)*60))

        str=f'{str} {progprcnt}% - {remHH}:{remMM:02d}:{remSS:02d} - {strend}'

        print(str)

        ti=tf*1

        itcont=itcont+1

        return [ti,deltat,itcont,str]

#}}}

