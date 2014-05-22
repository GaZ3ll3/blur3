

#----------------------------------------------------------
def read_params(outputdir,params):

    import string

    Fparams = "".join((outputdir,"/mesh_params.dat"     ))
    Rparams = open(Fparams,'r')

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[0] = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[1] = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[2] = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[3] = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[4] = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[5] = int(linelist[0])

    linestring = Rparams.readline()
    linelist = string.split(linestring)
    params[6] = int(linelist[0])

    Rparams.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_tnode(outputdir,params,tnode_full,tnode):

    import string
    
    Ftnode  = "".join((outputdir,"/mesh_tnode.dat"      ))
    Rtnode  = open(Ftnode, 'r')

    NumElems      = params[0]
    NumPhysElems  = params[1]
    
    for i in range(0,NumElems):
        linestring = Rtnode.readline()
        linelist = string.split(linestring)
        tnode_full[i,0] = int(linelist[0])-1
        tnode_full[i,1] = int(linelist[1])-1
        tnode_full[i,2] = int(linelist[2])-1

    for i in range(0,NumPhysElems):
        tnode[i,0] = tnode_full[i,0]
        tnode[i,1] = tnode_full[i,1]
        tnode[i,2] = tnode_full[i,2]

    Rtnode.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_node(outputdir,params,x_full,y_full,x,y):    

    import string

    Fnode   = "".join((outputdir,"/mesh_node.dat"       ))
    Rnode   = open(Fnode,  'r')

    NumNodes      = params[3]
    NumPhysNodes  = params[4]
    
    for i in range(0,NumNodes):
        linestring = Rnode.readline()
        linelist = string.split(linestring)
        x_full[i] = float(linelist[0])
        y_full[i] = float(linelist[1])
        
    for i in range(0,NumPhysNodes):
        x[i] = x_full[i]
        y[i] = y_full[i]
        
    Rnode.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_bnd(outputdir,params,bnd):

    import string

    Fbnode  = "".join((outputdir,"/mesh_bnd_node.dat"   ))
    Rbnode  = open(Fbnode, 'r')

    NumBndNodes   = params[5]
    
    for i in range(0,NumBndNodes):
        linestring = Rbnode.readline()
        linelist = string.split(linestring)
        bnd[i] = int(linelist[0])-1
    
    Rbnode.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_ghost(outputdir,params):

    import string

    Fghost  = "".join((outputdir,"/mesh_ghost_link.dat" ))
    Rghost  = open(Fghost, 'r')
    
    Rghost.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_area(outputdir,params):

    import string

    Farea   = "".join((outputdir,"/mesh_area_prim.dat"  ))
    Rarea   = open(Farea,  'r')
    
    Rarea.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_dual(outputdir,params):

    import string

    Fdual   = "".join((outputdir,"/mesh_area_dual.dat"  ))
    Rdual   = open(Fdual,  'r')
    
    Rdual.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_edge(outputdir,params):

    import string

    Fedge   = "".join((outputdir,"/mesh_edge.dat"       ))
    Redge   = open(Fedge,  'r')
    
    Redge.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_tedge(outputdir,params):

    import string

    Ftedge  = "".join((outputdir,"/mesh_tedge.dat"      ))
    Rtedge  = open(Ftedge, 'r')    

    Rtedge.close()
#----------------------------------------------------------


#----------------------------------------------------------
def read_eelem(outputdir,params):

    import string
    
    Feelem  = "".join((outputdir,"/mesh_eelem.dat"      ))
    Reelem  = open(Feelem, 'r')

    Reelem.close()
#----------------------------------------------------------


#----------------------------------------------------------
def node_minmax(NumNodes,x,y,xminmax):

    xminmax[0] = x[0]   #xmin
    xminmax[1] = x[0]   #xmax
    xminmax[2] = y[0]   #ymin
    xminmax[3] = y[0]   #ymax
    
    for i in range(1,NumNodes): 
        xtmp = x[i]        
        if xtmp < xminmax[0]:
            xminmax[0] = xtmp
        elif xtmp > xminmax[1]:
            xminmax[1] = xtmp
        
        ytmp = y[i]
        if ytmp < xminmax[2]:
            xminmax[2] = ytmp
        elif ytmp > xminmax[3]:
            xminmax[3] = ytmp

#----------------------------------------------------------
