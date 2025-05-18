
import numpy as np
from ase.io import read, write
from ase.build.supercells import make_supercell
import pandas as pd
from copy import deepcopy
from scipy.optimize import minimize
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import plotly.graph_objects as go
from plotly.figure_factory import create_trisurf
from colorsys import hsv_to_rgb





class pylro:
    
    def __init__(self,filename,atom_number,fileformat='vasp'):
        """
        Args:
            filename: a crystal structure file
            atom_number: the atomic number of the element chosen for miller plane analysis. Only supports unit n=1 in unit cell
            fileformat: ASE accepted file format
        
        Returns:
            pylro object.
        """
        
        
        
        self.struc = read(filename,format=fileformat)
        self.atom_number=atom_number
        original = np.array([[1,0,0],[0,1,0],[0,0,1]])
        self.struc = make_supercell(self.struc, original) 
        if atom_number==None: #lowest atomic number species by default to calculate periodicity
            unique, counts = np.unique(self.struc.numbers, return_counts=True)
            atom_number=unique[np.argmin(counts)]
        self.cell=self.struc.cell
        idx=[i for i,x in enumerate(self.struc.numbers) if x==atom_number]
        points_=self.struc.get_scaled_positions() #grab scaled, <abc> independent positions
        self.atom_locations=np.array([points_[i] for i in idx])
        

                


        

    def lattice_fit(self,n=12):
        """
        Creates a best fit lattice for the structure. Uses structure factor to determine periodicity among a,b,c directions.
        Args:
            n: the size limit of the sample supercell of lattice points
        Returns:
            self.dimensions: 1x3 array, size of cell
            self.lattice_repr: nx3 array, integer representation of atomic locations
            self.avgfac: average structure factor fit for all atoms. Perfect fit is number of atoms
            self.x,self.y,self.z: list of atomic locations in locations within supercell
            self.x_,self.y_,self.z: list of atomic locations in unit lattice representation
            
        """
        

        def structure_factor(pos, hkl):
            """ N*1 array"""
            F = 0
            h, k, l = hkl
            for xyz in pos:
                x,y,z = xyz
                F += np.exp(-2*np.pi*(1j)*(h*x + k*y+ l*z))

            return F
        
        
        def fit(dim,al,scaled=False):
            
            locs=al*dim
            locs=locs%1
            locs[locs>.5]=locs[locs>.5]-1
            locs=locs[np.argsort(locs)]
            centers=[locs[0]]
            for i,x in zip(range(1,len(locs)),locs[1:]):
                c=np.average(centers)
                periods=np.array([x+1,x,x-1])
                f=periods-c
                idx=np.argmin(np.abs(f))
                centers.append(periods[idx])

            shifted_al=al*dim-np.average(centers)
            lattice=np.round(shifted_al)

            if scaled:

                lattice+=np.average(centers)
                lattice/=dim
                return lattice
            
            lattice=np.array([int(x) for x in lattice])
            return lattice

        al=deepcopy(self.atom_locations)

        cdim=[]
        sfac=[]
        lattice_repr=[]
        basis=np.array([[1,0,0],[0,1,0],[0,0,1]])
        for j,b in enumerate(basis):
            ss=[np.abs(structure_factor(al,b*x)) for x in range(1,n)]
            ss_idx=np.argsort(-np.array(ss))
            continuous=False
            c=0
            while not continuous:
                
                dim=ss_idx[c]+1
                points=fit(dim,al[:,j])
                if (len(np.unique(points))<dim or dim==1) and (c!=n-2):
                    c+=1
                else:
                    if c==n-2:
                        dim=ss_idx[0]+1
                        c=0
                        if dim==1:
                            dim=ss_idx[1]+1
                            c=1
                        points=fit(dim,al[:,j])
                    cdim.append(dim)
                    sfac.append(ss[ss_idx[c]])
                    lattice_repr.append(points)
                    continuous=True

        self.dimensions=cdim
        self.lattice_repr=np.transpose(lattice_repr)
        self.avgfac=np.average(sfac)
    
        self.x,self.y,self.z=np.transpose(self.atom_locations)
        
        self.x_=fit(cdim[0],self.x,scaled=True)
        self.y_=fit(cdim[1],self.y,scaled=True)
        self.z_=fit(cdim[2],self.z,scaled=True)


        aa=np.unique(self.x_)[0:2]
        a_=aa[1]-aa[0]
        a_spacing=np.linalg.norm(a_*self.cell[0])
        
        bb=np.unique(self.y_)[0:2]
        b_=bb[1]-bb[0]
        b_spacing=np.linalg.norm(b_*self.cell[1])
        
        cc=np.unique(self.z_)[0:2]
        c_=cc[1]-cc[0]
        c_spacing=np.linalg.norm(c_*self.cell[2])
        
        self.d_spacings=[a_spacing,b_spacing,c_spacing]
        
        
        a_disorder=np.abs(self.x-self.x_)/a_
        b_disorder=np.abs(self.y-self.y_)/b_
        c_disorder=np.abs(self.z-self.z_)/c_
        
        self.disorders=[a_disorder,b_disorder,c_disorder]
        self.disorders=np.transpose(self.disorders)
        
        
        
        
        
        
        

                                
    def plane_order(self,plane,angstrom=False):
        """
        Calculates the order of an individual plane.
        Order is defined as the average unit lattice deviation from all atoms in any direction
        """
        
        plane=np.array(plane)/np.linalg.norm(plane) #must normalize plane so miller planes can be compared.

            
 
        s_=[np.dot(x,self.cell) for x in self.disorders]
        self.absolute_disorder=[np.abs(np.dot(x,np.array(plane))) for x in s_]
        self.relative_disorder=[np.abs(np.dot(x,np.array(plane))) for x in self.disorders]
        
        if angstrom:
            return np.average(self.absolute_disorder)
        return np.average(self.relative_disorder)
    
    def maximum_order(self,n=700):
        hkl=fibonacci_sphere(n)
        I=[self.plane_order(x) for x in hkl]
        mags=[np.linalg.norm(x) for x in I]
        return (np.min(mags), hkl[np.argmin(mags)])
    
    def minimum_order(self,n=700):
        hkl=fibonacci_sphere(n)
        I=[self.plane_order(x) for x in hkl]
        mags=[np.linalg.norm(x) for x in I]
        return (np.max(mags), hkl[np.argmax(mags)])
    
        # return np.average([1]) #The average unit lattice deviation in planar direction
    
    def miller_sphere_plot(self,n=700,c1=1.2,c2=10,cross_section=False,plot=True,angstrom=False):
        """Plotting function"""
        hkl=fibonacci_sphere(n)
        if not angstrom:
            I=[self.plane_order(x) for x in hkl]
        if angstrom:
            I=[self.plane_order(x,angstrom=True) for x in hkl]
            
        

        self.I=I
        if plot:
            O=Order_plot(hkl,I,c1,c2,angstrom=angstrom)
        if cross_section:
            O.plot_cross_sections()
        
        

    
    
class Order_plot(): 
    
    def __init__(self, hkl,I,c1=1.2,c2=10,angstrom=False):
        """
        Plots the data from LR_order
        Args:
            hkl: list of miller planes
            I: intensities from LR_order
            c1: controls relative peak intensities. Higher value exaggerates highest peaks more.
            c2: controls sphere size relative to peak heights. Higher value makes peaks smaller relative height.
        """
        
        #Load in Data
        self.hkl=hkl
        h=np.array([x[0] for x in hkl])
        k=np.array([x[1] for x in hkl])
        l=np.array([x[2] for x in hkl])
        I=np.array(I)

        #Prep Intensities for density calculation

        if angstrom:
            mmm=np.max(I)
        else:
            I=np.array([1-x for x in I])
            mmm=np.max(I)
 
        self.I=I


        

    
        sigma, n =.2 , 10000
        xyzs = fibonacci_sphere(n)
        grids = np.zeros([n, 3])
        grids[:, :2] = self.xyz2sph(xyzs)
        pts = []
        for i in range(len(h)):
            p, r = self.hkl2tp(h[i], k[i], l[i])
            pts.append([p, r, I[i]])
        pts = np.array(pts)
        vals = self.calculate_density(pts, xyzs, sigma=sigma)
        
        
        #Prep heights for sphere scaling
        valss=vals

        valss/=valss.max()
        valss*=mmm


        for i,x in enumerate(valss):
            xyzs[i]*=np.abs(x)

        
        phi=[]
        rho=[]
        for row in xyzs:
            r,p=self.hkl2tp(row[0],row[1],row[2])
            phi.append(p)
            rho.append(r)
        phi=np.array(phi)
        rho=np.array(rho)
        x=xyzs[:,0]
        y=xyzs[:,1]
        z=xyzs[:,2]
        
        self.x=x
        self.y=y
        self.z=z
        self.colorscale = [
            [0, "rgb(84,48,5)"],
            [1, "rgb(84,48,5)"],
        ]
        cmap=self.colormap_gen_(np.min(I),np.max(I))
        cscale_=[]
        for i in range(3):
            cscale_.append(tuple(x*255 for x in cmap[i]))
        cscale=[[0.,'rgb'+str(cscale_[0])],[.5,'rgb'+str(cscale_[1])],[1.,'rgb'+str(cscale_[2])]]
        if angstrom:
            cmap='Jet'   

        
        
        points2D=np.vstack([phi,rho]).T
        tri=Delaunay(points2D)
        simplices=tri.simplices
        layout = go.Layout(scene=dict(aspectmode='data',annotations=self.get_axis_names()))

        trisurf=create_trisurf(x=x,y=y,z=z,colormap=cmap, simplices=simplices,plot_edges=False,color_func=self.color_func_,show_colorbar=False)
  

        fig=go.Figure(data=trisurf, layout=layout)
        fig.add_trace(go.Scatter3d(x = [1.1,0,0], y = [0,1.1,0], z=[0,0,1.1], mode="text", text = ['a','b','c'],textfont=dict(size=29,family='Times New Roman')))
        self.add_axis_arrows(fig)
        fig.update_scenes(camera_projection_type='orthographic')
        fig.update_layout(
            scene=dict(
        xaxis=dict(
            tickvals=[-1,0,1],tickfont=dict(size=20,family='Times New Roman')  # Custom tick positions on the x-axis  # Custom tick labels
        ), 
        yaxis=dict(
            tickvals=[-1,0],tickfont=dict(size=20,family='Times New Roman')  # Custom tick positions on the y-axis
              # Custom tick labels
        ),
        zaxis=dict(
            tickvals=[1,0],tickfont=dict(size=20,family='Times New Roman')),
        xaxis_title='',
        yaxis_title='',
        zaxis_title=''))
       
        fig.update_layout(showlegend=False)
        mesh3dcbar=go.Mesh3d(x=[0,0],y=[0,0],z=[np.min(I),np.max(I)],intensity=z,showscale=True, colorscale=cscale,opacity=0,cmin=np.min(I),cmax=np.max(I),
                                colorbar=dict(len=.8,thickness=20,tickvals=np.round(np.linspace(np.min(I)*1.01,np.max(I)*.99,5),2),tickfont=dict(size=24,family='Times New Roman'),
                                title=dict(text='Order',font=dict(size=24,family='Times New Roman'),side='top'),x=.9))

        fig.add_trace(mesh3dcbar)
        fig.show(config={
            'displayModeBar': True,         
            'modeBarButtonsToRemove': ['toggleSpikelines', 'resetCameraDefault3d', 'hoverClosest3d', 'hoverClosestCartesian'],  
            'showTips': False})
        
        



    
    def colormap_gen_(self,low,high,threshold=.8):
        blue=np.array([0,0,1])
        yellow=np.array([1,1,0])
        red=np.array([1,0,0])
        mid=1-(1-threshold)/2
        if low<mid and high <=mid:
            x1=(low-threshold)/(mid-threshold)
            if x1<0:
                x1=0
            y1=1-x1
            x2=(high-threshold)/(mid-threshold)
            if x2<0:
                x2=0
            y2=1-x2
            c1=x1*yellow+y1*blue
            c2=x2*yellow+y2*blue
            c1[c1>1]=1
            c2[c2>1]=1
            
            return [tuple(c1),tuple(c2)]
            
            
        if low>mid and high>=mid:
            x1=(low-mid)/(1-mid)
            y1=1-x1
            x2=(high-mid)/(1-mid)
            y2=1-x2
            c1=x1*red+y1*yellow
            c2=x2*red+y2*yellow
            c1[c1>1]=1
            c2[c2>1]=1
            return [tuple(c1),tuple(c2)]
            
        if low<mid and high>mid:
            x1=(low-threshold)/(mid-threshold)
            if x1<0:
                x1=0
            y1=1-x1
            x2=(high-mid)/(1-mid)
            y2=1-x2
            c1=x1*yellow+y1*blue
            c2=x2*red+y2*yellow
            center=(low+high)/2
            if center>mid:
                x_=(center-mid)/(1-mid)
                y_=1-x_
                c_=x_*yellow+y_*red
                return [tuple(c1),tuple(c_),tuple(c2)]
            else:
                x_=(center-threshold)/(mid-threshold)
                y_=1-x_
                c_=x_*blue+y_*yellow
                return [tuple(c1),(1,1,0),tuple(c2)]
            
            
            
            
        
        

        
    def color_func(self,x,y,z):
        """
        Assigns color to distance
        """
        arr=np.array([x,y,z])
        arr_=[np.linalg.norm(arr-np.array(x)) for x in self.hkl]
        mag=self.I[np.argmax(arr_)]
        # mag=np.sqrt(x**2 + y**2 + z**2)
        # return np.floor(mag*255.9999)
        return mag
    def color_func_(self,x,y,z):
        """
        Assigns color to distance
        """

        mag=np.sqrt(x**2 + y**2 + z**2)
        # return np.floor(mag*255.9999)
        return mag

    def calculate_density(self,pts, xyzs, sigma=0.1):
        """
        calculate the projected order density on the unit sphere
        uses gaussain distrbution to smooth points.
        """
        vals = np.zeros(len(xyzs))
        pi = np.pi
        for pt in pts:
            t0, p0, h = pt
            x0, y0, z0 = np.sin(t0)*np.cos(p0), np.sin(t0)*np.sin(p0), np.cos(t0)
            dst = np.linalg.norm(xyzs - np.array([x0, y0, z0]), axis=1)
            vals += h*np.exp(-(dst**2/(2.0*sigma**2)))
        return vals

    def hkl2tp(self,h, k, l):
        """
        convert hkl to theta and phi
        """
        mp = [h,k,l]
        r = np.linalg.norm(mp)

        theta = np.arctan2(mp[1],mp[0])
        phi = np.arccos(mp[2]/r)

        #return theta, phi
        return phi, theta



    def xyz2sph(self,xyzs, radian=True):
        """
        convert the vectors (x, y, z) to the sphere representation (theta, phi)

        Args:
            xyzs: 3D xyz coordinates
            radian: return in radian (otherwise degree)
        """
        pts = np.zeros([len(xyzs), 2])   
        for i, r_vec in enumerate(xyzs):
            r_mag = np.linalg.norm(r_vec)
            theta0 = np.arccos(r_vec[2]/r_mag)
            if abs((r_vec[2] / r_mag) - 1.0) < 10.**(-8.):
                theta0 = 0.0
            elif abs((r_vec[2] / r_mag) + 1.0) < 10.**(-8.):
                theta0 = np.pi

            if r_vec[0] < 0.:
                phi0 = np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif 0. < r_vec[0] and r_vec[1] < 0.:
                phi0 = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif 0. < r_vec[0] and 0. <= r_vec[1]:
                phi0 = np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] == 0. and 0. < r_vec[1]:
                phi0 = 0.5 * np.pi
            elif r_vec[0] == 0. and r_vec[1] < 0.:
                phi0 = 1.5 * np.pi
            else:
                phi0 = 0.
            pts[i, :] = [theta0, phi0]
        if not radian:
            pts = np.degree(pts)

        return pts



    def get_arrow(self,axisname="x"):
        """
        Creates arrow object to plot axis lines
        """


        body = go.Scatter3d(
            marker=dict(size=1, color=self.colorscale[0][1]),
            line=dict(color=self.colorscale[0][1], width=3),
            showlegend=False,  # hide the legend
        )

        head = go.Cone(
            sizeref=0.1,
            autocolorscale=None,
            colorscale=self.colorscale,
            showscale=False,  # disable additional colorscale for arrowheads
            hovertext=axisname,
        )
        for ax, direction in zip(("x", "y", "z"), ("u", "v", "w")):
            if ax == axisname:
                body[ax] = -1,1
                head[ax] = [1]
                head[direction] = [1]
            else:
                body[ax] = 0,0
                head[ax] = [0]
                head[direction] = [0]

        return [body, head]


    def add_axis_arrows(self,fig):
        for ax in ("x", "y", "z"):
            for item in self.get_arrow(ax):
                fig.add_trace(item)

    def get_annotation_for_ax(self,ax):
        """
        plots abc axis labels
        """
        d = dict(showarrow=False, text=ax, xanchor="left", font=dict(color="#1f1f1f",size=28))

        if ax == "a":
            d["x"] = 1.1
            d["y"] = 0
            d["z"] = 0
        elif ax == "b":
            d["x"] = 0
            d["y"] = 1.1 
            d["z"] = 0
        else:
            d["x"] = 0
            d["y"] = 0
            d["z"] = 1.1 

        if ax in {"a", "b"}:
            d["xshift"] = 15

        return d


    def get_axis_names(self):
        return [self.get_annotation_for_ax(ax) for ax in ("a", "b", "c")]
    
    def plot_cross_sections(self):
        """
        Plots the axis plane cross sections of the plot 3D visualization
        Uses scattering of points under a limit as fibonacci sphere doesn't points distributed in a plane.
        """
        x=self.x
        y=self.y
        z=self.z
        
        e=2e-2
        xy_a=[x_ for i,x_ in enumerate(x) if np.abs(z[i])<e]
        xy_b=[z_ for i,z_ in enumerate(y) if np.abs(z[i])<e]
        xz_a=[x_ for i,x_ in enumerate(x) if np.abs(y[i])<e]
        xz_b=[z_ for i,z_ in enumerate(z) if np.abs(y[i])<e]
        yz_a=[x_ for i,x_ in enumerate(y) if np.abs(x[i])<e]
        yz_b=[z_ for i,z_ in enumerate(z) if np.abs(x[i])<e]

        fig,ax=plt.subplots(nrows=1,ncols=3,figsize=(15,5))
        # plt.gca().set_aspect('equal', adjustable='box')
        ax[0].scatter(xy_a,xy_b,color='r')
        ax[0].set_xlim(-1,1)
        ax[0].set_ylim(-1,1)
        ax[0].set_title('X-Y')
        ax[0].set_xlabel('X')
        ax[0].set_ylabel('Y')
        ax[1].scatter(xz_a,xz_b,color='r')
        ax[1].set_title('X-Z')
        ax[1].set_xlabel('X')
        ax[1].set_ylabel('Z')
        ax[1].set_xlim(-1,1)
        ax[1].set_ylim(-1,1)
        ax[2].scatter(yz_a,yz_b,color='r')
        ax[2].set_title('Y-Z')
        ax[2].set_xlabel('Y')
        ax[2].set_ylabel('Z')
        ax[2].set_xlim(-1,1)
        ax[2].set_ylim(-1,1)
        plt.show()
                           
def fibonacci_sphere(samples=1000):
    """
    Sampling the sphere grids

    Args:
        samples: number of pts to generate

    Returns:
        3D points array in Cartesian coordinates
    """
    points = []
    phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians
    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = np.sqrt(1 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append((x, y, z))

    return np.array(points)

