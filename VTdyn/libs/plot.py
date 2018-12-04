import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, PatchCollection
import matplotlib.patches as patches
import seaborn as sns
from shapely.ops import polygonize
from shapely.geometry import LineString, MultiPolygon, Polygon, MultiPoint, Point, LinearRing
from scipy.spatial import Voronoi, Delaunay
from descartes.patch import PolygonPatch
import os

current_palette = sns.color_palette()

#functions for plotting and animating tissue objects/histories


def plot_tri_torus(tissue,ax=None,time = None,label=False,palette=current_palette):
    """plot Delaunay triangulation of a tissue object (i.e. cell centres and neighbour connections) with torus geometry"""
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    width, height = tissue.mesh.geometry.width, tissue.mesh.geometry.height 
    centres = tissue.mesh.centres
    centres_3x3 = np.vstack([centres+[dx, dy] for dx in [-width, 0, width] for dy in [-height, 0, height]])
    tri = Delaunay(centres_3x3).simplices  
    plt.triplot(centres_3x3[:,0], centres_3x3[:,1], tri.copy(),color=palette[3])
    plt.plot(centres_3x3[:,0], centres_3x3[:,1], 'o',color = 'black')
    if label:
        for i, coords in enumerate(tissue.mesh.centres):
            plt.text(coords[0],coords[1],str(i))
    if time is not None:
        lims = plt.axis()
        plt.text(lims[0]+0.1,lims[3]+0.1,'t = %.2f hr'%time)
    plt.show()



def torus_plot(tissue,palette=np.array(current_palette),key=None,key_label=None,ax=None,show_centres=False,cell_ids=False,mesh_ids=False,areas=False,boundary=False,colours=None,animate=False):
    """plot tissue object with torus geometry
    args
    ---------------
    tissue: to plot
    palette: array of colors
    key: (str) color code cells by key. must match a tissue.properties key, e.g. type, ancestors
    key_label: (bool) plot key vals if True
    ax: matplotlib axes
    show_centres: (bool) plot cell centres if True
    cell_ids/mesh_ids: (bool) label with cell/mesh ids if True
    areas: (bool) label with cell areas if True
    boundary: (bool) if True plot square boundary over which tissue is periodic
    colours: (array) provide colours for specific key values (if key not None) OR colour of each cell 
    animate: (bool) True if plot is part of animation
    """
    width, height = tissue.mesh.geometry.width, tissue.mesh.geometry.height 
    centres = tissue.mesh.centres 
    centres_3x3 = np.vstack([centres+[dx, dy] for dx in [-width, 0, width] for dy in [-height, 0, height]])
    N = tissue.mesh.N_mesh
    mask = np.full(9*N,False,dtype=bool)
    mask[4*N:5*N]=True
    vor = Voronoi(centres_3x3)
    
    mp = MultiPolygon([Polygon(vor.vertices[region])
                for region in (np.array(vor.regions)[np.array(vor.point_region)])[mask]])
    
    if ax is None: 
        fig = plt.figure()
        ax = fig.add_subplot(111)
        minx, miny, maxx, maxy = -width/2,-height/2,width/2,height/2
        w, h = maxx - minx, maxy - miny
        ax.set_xlim(minx - 0.2 * w, maxx + 0.2 * w)
        ax.set_ylim(miny - 0.2 * h, maxy + 0.2 * h)
        ax.set_aspect(1)
        ax.set_axis_off()
    ax.set_xticks([])
    ax.set_yticks([])
    
    if key is None and colours is None: ax.add_collection(PatchCollection([PolygonPatch(p,linewidth=2.5,edgecolor='black') for p in mp],match_original=True))
    elif key is not None:
        if colours is None: colours = palette[tissue.properties[key]]
        coll = PatchCollection([PolygonPatch(p,facecolor = c,linewidth=2.5,edgecolor='black') for p,c in zip(mp,colours)],match_original=True)
        ax.add_collection(coll)
    elif colours is not None:
        coll = PatchCollection([PolygonPatch(p,facecolor = colours,linewidth=2.5,edgecolor='black') for p in mp],match_original=True)
        ax.add_collection(coll)
    
    if show_centres: 
        plt.plot(centres[:,0], centres[:,1], 'o',color='black')
    if cell_ids:
        ids = tissue.cell_ids
        for i, coords in enumerate(tissue.mesh.centres):
            plt.text(coords[0],coords[1],str(ids[i]))
    if mesh_ids:
        for i, coords in enumerate(tissue.mesh.centres):
            plt.text(coords[0],coords[1],str(i))
    if areas:
        for area, coords in zip(tissue.mesh.areas,tissue.mesh.centres):
            plt.text(coords[0],coords[1],'%.2f'%area)
    if key_label is True:
        for id, coords in zip(tissue.properties[key],tissue.mesh.centres):
            plt.text(coords[0],coords[1],str(id))
    if boundary: 
        ax.add_patch(patches.Rectangle((-width/2,-height/2),width,height,fill=False,linewidth=1.5))
    if not animate: plt.show()

def save_mpg_torus(history, name, index=None,key = None, delete_images=True):
    """save animation of tissue history with torus geometry
    args: 
    --------------- 
    history: list of tissue objects
    name: (str) output filename
    index: (int) if not None add index to output filename
    key: (str) color code cells by key. must match a tissue.properties key, e.g. type, ancestors
    delete_images: (bool) if True delete image files after creating video
    """
    outputdir="images"
    if not os.path.exists(outputdir): # if the folder doesn't exist create it
        os.makedirs(outputdir)
    width,height = history[0].mesh.geometry.width,history[0].mesh.geometry.height
    fig = plt.figure()
    ax = fig.add_subplot(111)
    minx, miny, maxx, maxy = -width/2,-height/2,width/2,height/2
    w, h = maxx - minx, maxy - miny
    ax.set_xlim(minx - 0.2 * w, maxx + 0.2 * w)
    ax.set_ylim(miny - 0.2 * h, maxy + 0.2 * h)
    ax.set_aspect(1)
    plt.ion()
    sns.set_style('white')
    ax.set_axis_off()
    ax.set_autoscale_on(False)
    frames=[]
    i = 0
    if key is not None:
        key_max = max((max(tissue.properties[key]) for tissue in history))
        if key_max>6:
            palette = np.array(sns.color_palette("husl", key_max+1))
            np.random.shuffle(palette)
        else: palette = np.array(current_palette)
        for tissue in history:
            try:
                for c in ax.collections:
                    c.remove()
            except TypeError: pass
            torus_plot(tissue,palette,key,ax=ax,animate=True)
            frame="images/image%05i.png" % i
            fig.savefig(frame,dpi=500)
            frames.append(frame)
            i+=1
    else:
        for tissue in history:
            try:
                for c in ax.collections:
                    c.remove()
            except TypeError: pass
            torus_plot(tissue,ax=ax,animate=True)
            frame="images/image%05i.png" % i
            fig.savefig(frame,dpi=500)
            frames.append(frame)
            i+=1
    if index is not None: os.system("mencoder 'mf://images/image*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "%s%0.3f.mpg" %(name,index))
    else: os.system("mencoder 'mf://images/image*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "%s.mpg" %name)
    if delete_images:
        for frame in frames: os.remove(frame)



def animate_torus(history, key = None):
    """view animation of tissue history with torus geometry
    args: 
    --------------- 
    history: list of tissue objects
    key: (str) color code cells by key. must match a tissue.properties key, e.g. type, ancestors
    """
    width,height = history[0].mesh.geometry.width,history[0].mesh.geometry.height
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    minx, miny, maxx, maxy = -width/2,-height/2,width/2,height/2
    w, h = maxx - minx, maxy - miny
    ax.set_xlim(minx - 0.2 * w, maxx + 0.2 * w)
    ax.set_ylim(miny - 0.2 * h, maxy + 0.2 * h)
    ax.set_aspect(1)
    sns.set_style('white')
    ax.set_axis_off()
    ax.set_autoscale_on(False)
    plot = []
    if key is not None:
        key_max = max((max(tissue.properties[key]) for tissue in history))
        if key_max>6:
            palette = np.array(sns.color_palette("husl", key_max+1))
            np.random.shuffle(palette)
        else: palette = np.array(current_palette)
        for tissue in history:
            try:
                for c in ax.collections:
                    c.remove()
            except TypeError: pass
            torus_plot(tissue,palette,key,ax=ax)
            plt.pause(0.001)
    else:
        for tissue in history:
            try:
                for c in ax.collections:
                    c.remove()
            except TypeError: pass
            torus_plot(tissue,ax=ax)
            plt.pause(0.001)#



#plotting functions for a PLANAR geometry tissue object IN DEVELOPMENT   
    
# def finite_plot(tissue,palette=current_palette,key=None,key_label=None,ax=None,show_centres=False,cell_ids=False,mesh_ids=False):
#     centres = tissue.mesh.centres
#     vor = tissue.mesh.voronoi()
#     center = vor.points.mean(axis=0)
#     ptp_bound = vor.points.ptp(axis=0)
#
#     regions = [Polygon(vor.vertices[region]) if -1 not in region else
#                 Polygon(get_region_for_infinite(region,vor,center,ptp_bound))
#                 for region in np.array(vor.regions)[np.array(vor.point_region)] if len(region)>=2
#                 ]
#     convex_hull = MultiPoint([Point(i) for i in centres]).convex_hull
#     mp = MultiPolygon(
#         [poly.intersection(convex_hull) for poly in regions])
#
#     if ax is None:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
#         minx, miny, maxx, maxy = mp.bounds
#         w, h = maxx - minx, maxy - miny
#         ax.set_xlim(minx - 0.2 * w, maxx + 0.2 * w)
#         ax.set_ylim(miny - 0.2 * h, maxy + 0.2 * h)
#         ax.set_aspect(1)
#
#     if key is None: ax.add_collection(PatchCollection([PolygonPatch(p) for p in mp]))
#     else:
#         colours = palette[tissue.by_mesh(key)]
#         coll = PatchCollection([PolygonPatch(p,facecolor = c) for p,c in zip(mp,colours)],match_original=True)
#         # coll.set_facecolors(palette[tissue.by_mesh(key)])
#         ax.add_collection(coll)
#
#     if show_centres:
#         plt.plot(centres[:,0], centres[:,1], 'o',color='black')
#     if cell_ids:
#         ids = tissue.by_mesh('id')
#         for i, coords in enumerate(tissue.mesh.centres):
#             plt.text(coords[0],coords[1],str(ids[i]))
#     if mesh_ids:
#         for i, coords in enumerate(tissue.mesh.centres):
#             plt.text(coords[0],coords[1],str(i))
#     if key_label is not None:
#         ids = tissue.by_mesh(key_label)
#         for i, coords in enumerate(tissue.mesh.centres):
#             plt.text(coords[0],coords[1],str(ids[i]))
#
#
# def get_vertices_for_inf_line(vor,i,pointidx,center,ptp_bound):
#     t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
#     t /= np.linalg.norm(t)
#     n = np.array([-t[1], t[0]])  # normal
#
#     midpoint = vor.points[pointidx].mean(axis=0)
#     direction = np.sign(np.dot(midpoint - center, n)) * n
#     far_point = vor.vertices[i] + direction * ptp_bound.max()
#
#     return far_point
#
# def get_region_for_infinite(region,vor,center,ptp_bound):
#     infidx= region.index(-1)
#     finite_end1 = region[infidx-1]
#     finite_end2 = region[(infidx+1)%len(region)]
#
#     try: pointidx1 = vor.ridge_points[vor.ridge_vertices.index([-1,finite_end1])]
#     except ValueError: pointidx1 = vor.ridge_points[vor.ridge_vertices.index([finite_end1,-1])]
#     far_point1 = get_vertices_for_inf_line(vor,finite_end1,pointidx1,center,ptp_bound)
#
#     try: pointidx2 = vor.ridge_points[vor.ridge_vertices.index([-1,finite_end2])]
#     except ValueError: pointidx2 = vor.ridge_points[vor.ridge_vertices.index([finite_end2,-1])]
#     far_point2 = get_vertices_for_inf_line(vor,finite_end2,pointidx2,center,ptp_bound)
#     region_vertices = []
#     for pt in region:
#         if pt != -1: region_vertices.append(vor.vertices[pt])
#         else: region_vertices.append(far_point1); region_vertices.append(far_point2)
#     return np.array(region_vertices)
#
#
# def plot_cells(tissue,current_palette=current_palette,key=None,ax=None,label=False,time = False,colors=None,centres=True):
#     ghosts = tissue.by_mesh('ghost')
#     fig = plt.Figure()
#     if ax is None:
#         ax = plt.axes()
#         plt.axis('scaled')
#         xmin,xmax = min(tissue.mesh.centres[:,0]), max(tissue.mesh.centres[:,0])
#         ymin,ymax = min(tissue.mesh.centres[:,1]), max(tissue.mesh.centres[:,1])
#         ax.set_xlim(xmin,xmax)
#         ax.set_ylim(ymin,ymax)
#         ax.xaxis.set_major_locator(plt.NullLocator())
#         ax.yaxis.set_major_locator(plt.NullLocator())
#     ax.cla()
#     vor = tissue.mesh.voronoi()
#     cells_by_vertex = np.array(vor.regions)[np.array(vor.point_region)]
#     verts = [vor.vertices[cv] for cv in cells_by_vertex[~ghosts]]
#     if colors is not None:
#         coll = PolyCollection(verts,linewidths=[2.],facecolors=colors)
#     elif key is not None:
#         colors = np.array(current_palette)[tissue.by_mesh(key)]
#         coll = PolyCollection(verts,linewidths=[2.],facecolors=colors)
#     else: coll = PolyCollection(verts,linewidths=[2.])
#     ax.add_collection(coll)
#     if label:
#         ids = tissue.by_mesh('id')
#         for i, coords in enumerate(tissue.mesh.centres):
#             if ~ghosts[i]: plt.text(coords[0],coords[1],str(ids[i]))
#     if time:
#         lims = plt.axis()
#         plt.text(lims[0]+0.1,lims[3]+0.1,'t = %.2f hr'%time)
#     if centres:
#         real_centres = tissue.mesh.centres[~ghosts]
#         plt.plot(real_centres[:,0], real_centres[:,1], 'o',color='black')
#
#
# def animate_finite(history, key = None, timestep=None):
#     xmin,ymin = np.amin([np.amin(tissue.mesh.centres,axis=0) for tissue in history],axis=0)*1.5
#     xmax,ymax = np.amax([np.amax(tissue.mesh.centres,axis=0) for tissue in history],axis=0)*1.5
#
#     plt.ion()
#     fig = plt.figure()
#     ax = plt.axes()
#     plt.axis('scaled')
#     ax.set_xlim(xmin,xmax)
#     ax.set_ylim(ymin,ymax)
#     fig.set_size_inches(6, 6)
#     ax.set_autoscale_on(False)
#     plot = []
#     if key is not None:
#         key_max = max((max(tissue.by_mesh(key)) for tissue in history))
#         palette = np.array(sns.color_palette("husl", key_max+1))
#         np.random.shuffle(palette)
#         for tissue in history:
#             ax.cla()
#             finite_plot(tissue,palette,key,ax=ax)
#             plt.pause(0.001)
#     else:
#         for tissue in history:
#             ax.cla()
#             finite_plot(tissue,ax=ax)
#             plt.pause(0.001)
#
#
# def animate_video_mpg(history,name,index,facecolours='Default'):
#     v_max = np.max((np.max(history[0].mesh.centres), np.max(history[-1].mesh.centres)))
#     if key: key_max = np.max(history[0].properties[key])
#     size = 2.0*v_max
#     outputdir="images"
#     if not os.path.exists(outputdir): # if the folder doesn't exist create it
#         os.makedirs(outputdir)
#     fig = plt.figure()
#     ax = plt.axes()
#     plt.axis('scaled')
#     lim = [-0.55*size, 0.55*size]
#     ax.set_xlim(lim)
#     ax.set_ylim(lim)
#     frames=[]
#     i = 0
#     for cells in history:
#         plot_cells(cells,key,ax)
#         i=i+1
#         frame="images/image%04i.png" % i
#         fig.savefig(frame,dpi=500)
#         frames.append(frame)
#     os.system("mencoder 'mf://images/image*.png' -mf type=png:fps=20 -ovc lavc -lavcopts vcodec=wmv2 -oac copy  -o " + "%s%0.3f.mpg" %(name,index))
#     for frame in frames: os.remove(frame)
#