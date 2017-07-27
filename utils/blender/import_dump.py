#!/usr/bin/env python3

#tool for rendering 1cpn in blender

import bpy
import math
import copy
from mathutils import Vector
from mathutils import Quaternion
import pdb

path_zewdie = "/media/int_3TB/Work/dna_chromatin/1cpn/1cpn-model/utils/blender/zewdie.stl"
#path_traj = "/media/int_3TB/Work/dna_chromatin/1cpn/viz/trunk/blender/dump/in.dump"
path_traj = "/tmp/frame.dump"

#path_zewdie = "/home/lequieu/Work/depablo/1cpn/viz/trunk/blender/utils/zewdie.stl"
#path_traj = "/home/lequieu/Work/depablo/1cpn/viz/trunk/blender/simple.dump"

SCALEFACTOR=0.05
# if true, draw bonds within nucl, and create different materials
DRAWNDNA = True
DYADBONDS = False
NDNABONDS = False
DNABONDS = True

#                  R,G,B,alpha 
COLORS={'fvec'  :(1,0.05,0.05, 1), #red
        'vvec'  :(0.05,1,0.05, 1), # green
        'uvec'  :(0.05,0.05,1, 1), # dark blue
        'bead'  :(0.111,0.264,1, 1), 
        'nucldna'  :(0.111,0.264,1, 1.0),
#       'dyad' :(1.0,0.386,0.145, 1),
        'dyad' :(1.0,0.921,0.0, 1),
        'plane' :(1, 1, 1, 1),
        'worldbg' :(1, 1, 1, 1),
        'zewdie':    (1,0.194,0.194, 0.2),
        'zewdiecore':(1,0.194,0.194, 1.0),
        'bond-12':(0.892, 0.282, 0.051, 1.0) , #orange
        'bond-22-nucl':(0.046, 0.500, 0.045, 1.0), # green
#        'bond-22-nucl':(0.216, 0.013, 1.000, 1.0),
        'bond-22':(0.111,0.264,1, .7), #blue 
        'bond-22-main':(0.111,0.264,1, 1.0), #blue 
        'bond-22-secondary':(.134,.5, .692, 1.0), #light blue
        'bond-13':(0.049, 0.466, 0.805, 1.0) , #light blue
#        'bond-23':(0.892, 0.282, 0.051, 1.0) , 
        'bond-23':(0.216, 0.013, 1.000, 1.0)} #purple
#        'zewdie':    (0.575,0.057,1, 0.2),
#        'zewdie':(0.05,0.50,1, 0.2)}

try:
    import bpy
    bpy_avail = True
except ImportError:
    bpy_avail = False

def quat_conj(q): # return conj(q)
    return [q[0], -q[1], -q[2], -q[3]]

def quat_vec_rot(src,q):
  # quat_vec_rot(vector &dest, vector &src, quaternion &q) dest = q*src*conj(q)
  #dest = [0]*3
  dest = Vector([0,0,0])
  aa=[q[0]*q[0], q[1]*q[1], q[2]*q[2], q[3]*q[3]]
  ab=[q[0]*q[1], q[0]*q[2], q[0]*q[3], q[1]*q[2], q[1]*q[3], q[2]*q[3]]
  dest[0] = (aa[0]+aa[1]-aa[2]-aa[3])*src[0]+\
            ((ab[3]-ab[2])*src[1]+(ab[1]+ab[4])*src[2])*2.0;
  dest[1] = (aa[0]-aa[1]+aa[2]-aa[3])*src[1]+\
            ((ab[2]+ab[3])*src[0]+(ab[5]-ab[0])*src[2])*2.0;
  dest[2] = (aa[0]-aa[1]-aa[2]+aa[3])*src[2]+\
            ((ab[4]-ab[1])*src[0]+(ab[0]+ab[5])*src[1])*2.0;
  return dest;


def quat_multiply(p,q): # out = p*q
    #out = [0]*4
    out = Vector([0,0,0,0])
    out[0] =  p[0]*q[0]-p[1]*q[1]-p[2]*q[2]-p[3]*q[3];
    out[1] =  p[0]*q[1]+p[1]*q[0]+p[2]*q[3]-p[3]*q[2];
    out[2] =  p[0]*q[2]-p[1]*q[3]+p[2]*q[0]+p[3]*q[1];
    out[3] =  p[0]*q[3]+p[1]*q[2]-p[2]*q[1]+p[3]*q[0];
    return out

#def tu2rotquat(theta,u): # theta and u to rotational quat
#    u[0] = u[0]/vector_norm(u)
#    u[1] = u[1]/vector_norm(u)
#    u[2] = u[2]/vector_norm(u)
#    return [math.cos(theta/2.0),\
#               u[0]*math.sin(theta/2.0),\
#               u[1]*math.sin(theta/2.0),\
#               u[2]*math.sin(theta/2.0)]
#
#def rotquat2tu(q): # rotational quat to theta and u
#    theta = 2.0*math.acos(q[0])
#    if theta < 1e-3: theta = 1e-3
#    u = [0]*3
#    u[0] = q[1] / math.sin(theta/2.0)
#    u[1] = q[2] / math.sin(theta/2.0)
#    u[2] = q[3] / math.sin(theta/2.0)
#    return theta,u

class Bond:
    __slots__ = ('id','x1','x2','q1','q2','typenameA','typenameB','drawmode')
    def __init__(self,x1,x2,typenameA,q1=None,q2=None,typenameB=None):
        self.x1 = x1
        self.x2 = x2
        self.typenameA = typenameA
        if q1 != None and q2 != None and typenameB != None:
          self.q1 = q1
          self.q2 = q2
          self.typenameB = typenameB
          self.drawmode = "draw_dna_bond"
        elif q1 == None and q2 == None and typenameB == None:
          self.drawmode = "draw_bond"
        else:
          print('ERROR: Invalid bond constructor')
          exit(1)

class DumpReader:
    def __init__(self,fnme_):
        self.fnme = fnme_
        self.is_first_frame = True;
        self.f = open(self.fnme,"r")
        self.nframes = 0
        while (self.read_next_frame() == 0):
          self.nframes += 1
        self.f.close() 
        self.f = open(self.fnme,"r")
        #self.x = []
        #self.type = []
        #self.q = []
    def close(self):
        self.f.close()

    def compute_bonds(self):
      self.bonds = []
      for i in range(0,self.natoms-1):
          t1 = t2 = t3 = t4 = t5 = -1
          x1 = x2 = x3 = x4 = x5 = (0,0,0)

          x1 = self.x[i]
          t1 = self.type[i]
          q1 = self.q[i] 
          x2 = self.x[i+1]
          t2 = self.type[i+1]
          q2 = self.q[i+1] 
          if i+2 < self.natoms:
            x3 = self.x[i+2]
            t3 = self.type[i+2]
          if i+3 < self.natoms:
            x4 = self.x[i+3]
            t4 = self.type[i+3]

          if DYADBONDS:
            if (t1 == 1 and t2 == 3): self.bonds.append(Bond(x1, x2, 'bond-13'))
            if (t1 == 3 and t2 == 2): self.bonds.append(Bond(x1, x2, 'bond-23'))
            if (t1 == 2 and t3 == 3): self.bonds.append(Bond(x1, x3, 'bond-23'))
          if NDNABONDS:
            if (t1 == 2 and t2 == 1): self.bonds.append(Bond(x1, x2, 'bond-12'))
            if (t1 == 1 and t3 == 2): self.bonds.append(Bond(x1, x3, 'bond-12'))
            if (t1 == 2 and t2 == 1 and t4 == 2): self.bonds.append(Bond(x1, x4, 'bond-22-nucl'))
          if DNABONDS:
            #if (t1 == 2 and t2 == 2): draw_bond(x1, x2, 'bond-22', rcyl=5)
            if (t1 == 2 and t2 == 2): 
              self.bonds.append(Bond(x1,x2,'bond-22-main',q1=self.q[i],q2=self.q[i+1],typenameB='bond-22-secondary'))

      self.nbonds = len(self.bonds)
      #print("nbonds: %d" % self.nbonds)
        
    def read_next_frame(self):
        f = self.f
        while True:
          line=f.readline()
          if "ITEM: TIMESTEP" in line:
            line=f.readline()
            self.timestep = int(line)
          elif "ITEM: NUMBER OF ATOMS" in line:
            line=f.readline()
            n = int(line)
            if (not self.is_first_frame) and (self.natoms != n):
                print("Warning! natoms changes between frames. Make sure this is correct!")
            self.natoms = int(line)
          elif "ITEM: BOX BOUNDS" in line:
            l=f.readline().split()
            self.xboxmin = float(l[0])
            self.xboxmax = float(l[1])
            l=f.readline().split()
            self.yboxmin = float(l[0])
            self.yboxmax = float(l[1])
            l=f.readline().split()
            self.zboxmin = float(l[0])
            self.zboxmax = float(l[1])
          elif "ITEM: ATOMS" in line:
            if self.is_first_frame:
                self.x = [Vector([0,0,0]) for i in range(self.natoms)]
                self.q = [Quaternion([0,0,0,0]) for i in range(self.natoms)]
                self.type = ['']*self.natoms

            for i in range(self.natoms):
                l = f.readline().split()
                index = int(l[0]) - 1
                self.type[index] = int(l[1])
                self.x[index] = Vector([float(l[2]), float(l[3]), float(l[4])])
                self.q[index][0] = float(l[5])
                self.q[index][1] = float(l[6])
                self.q[index][2] = float(l[7])
                self.q[index][3] = float(l[8])
            
            #compute bonds
            self.compute_bonds()
            
            # calc avg x
            self.xavg = Vector([0,0,0])
            for i in range(self.natoms):
                self.xavg[0] += self.x[i][0]
                self.xavg[1] += self.x[i][1]
                self.xavg[2] += self.x[i][2]
            self.xavg /= self.natoms 

            if self.is_first_frame:
                self.x0 = copy.copy(self.x)
                self.q0 = copy.copy(self.q)
                self.is_first_frame = False

            break
          elif '' == line:
            print("Dump File Ended")
            return 1
        #print("Finished reading frame")
        return 0
    
# create_material - Creates a material with a color if it doesn't exist already
def create_material(name, color):
    #group = bpy.data.node_groups.get(group_name)
    #if group is None:
    #    group = create_node_group(group_name, color)

    mat = bpy.data.materials.get(name)
    if mat is not None:
      return mat

    mat = bpy.data.materials.new(name)
    mat.use_nodes = True

    for node in mat.node_tree.nodes:
        mat.node_tree.nodes.remove(node)

    # Create diffuse shader.
    diffuse = mat.node_tree.nodes.new('ShaderNodeBsdfDiffuse')
    diffuse.location = (0, 600)
    diffuse.inputs[0].default_value = (color[0], color[1], color[2], 1)
    diffuse.name = 'diffuse'

    # Create transparent shader.
    transparent = mat.node_tree.nodes.new('ShaderNodeBsdfTransparent')
    transparent.location = (0, 400)
    transparent.name = 'transparent'

    # Create shader mixer.
    mixer = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer.inputs[0].default_value = color[3]
    mixer.location = (200, 600)
    mixer.name = 'mixer'

    # Create input and output sockets.
    #mat.node_tree.inputs.new('NodeSocketFloat', 'Transparency')
    #mat.node_tree.outputs.new('NodeSocketShader', 'Shader')
    
    #bpy.data.materials['uvec'].node_tree.nodes['Diffuse BSDF'].inputs[0].default_value = (1,1,1,1)

    outmat = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
    outmat.location = (400, 600)
    outmat.name = 'outmat'

    # Link everything together.
    mat.node_tree.links.new(transparent.outputs[0], mixer.inputs[1])
    mat.node_tree.links.new(diffuse.outputs[0], mixer.inputs[2])
    mat.node_tree.links.new(mixer.outputs[0], outmat.inputs[0])

    return mat

def create_material_emission_glossy(name,color):

    mat = bpy.data.materials.get(name)
    if mat is not None:
      return mat

    mat = bpy.data.materials.new(name)
    mat.use_nodes = True

    for node in mat.node_tree.nodes:
        mat.node_tree.nodes.remove(node)

    #create rgb
    rgb = mat.node_tree.nodes.new('ShaderNodeRGB')
    rgb.location = (0,600)
    rgb.outputs[0].default_value = (color[0], color[1], color[2], 1)
    
    # Create emission shader.
    emission = mat.node_tree.nodes.new('ShaderNodeEmission')
    emission.location = (0, 400)
    emission.inputs[1].default_value = 5.0
    emission.name = 'emission'

    # Create glossy shader.
    glossy = mat.node_tree.nodes.new('ShaderNodeBsdfGlossy')
    glossy.location = (200, 600)
    glossy.name = 'glossy'

    # Create shader mixer.
    mixer = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer.inputs[0].default_value = 0.5
    mixer.location = (400, 600)
    mixer.name = 'mixer'

    #how to change value 
    #bpy.data.materials['uvec'].node_tree.nodes['Diffuse BSDF'].inputs[0].default_value = (1,1,1,1)

    outmat = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
    outmat.location = (400, 0)
    outmat.name = 'outmat'

    # Link everything together.
    mat.node_tree.links.new(rgb.outputs[0], emission.inputs[0])
    mat.node_tree.links.new(rgb.outputs[0], glossy.inputs[0])
    mat.node_tree.links.new(emission.outputs[0], mixer.inputs[1])
    mat.node_tree.links.new(glossy.outputs[0], mixer.inputs[2])
    mat.node_tree.links.new(mixer.outputs[0], outmat.inputs[0])

    return mat


def create_material_holograph():
    '''This is a complex material, see the hologram .blend file for an example to play with different materials '''

    name = 'hologram'
    #red
    val = {'layerweight': 0.05, 'transparent1': (1,.262,.282,1), 'emission1': ((.8,.3,.341,1),2), 'wireframe': 0.3, 'volabsorp': ((.8,.414,.414,1),1), 'emission2': ((.8,.51,.459,1),.5),'mixer4': 0.5}
    #original blue
    #val = {'layerweight': 0.05, 'transparent1': (.566,.708,1,1), 'emission1': ((.298,.676,.8,1),2), 'wireframe': 0.3, 'volabsorp': ((.255,.308,.8,1),1), 'emission2': ((.298,.676,.8,1),0.5), 'mixer4': 0.5}

    mat = bpy.data.materials.get(name)
    if mat is not None:
        return mat #if mat already exists, just return it

    mat = bpy.data.materials.new(name)
    mat.use_nodes = True
    for node in mat.node_tree.nodes:
        mat.node_tree.nodes.remove(node)

    #mat.node_tree = bpy.data.node_mat.node_trees.new(name,'ShaderNodeTree')


    layerweight = mat.node_tree.nodes.new('ShaderNodeLayerWeight')
    layerweight.location = (0,800)
    layerweight.inputs[0].default_value = val['layerweight'] #0.150 #blend
    layerweight.name = 'layerweight'

    transparent1 = mat.node_tree.nodes.new('ShaderNodeBsdfTransparent')
    transparent1.location = (0,600)
    transparent1.inputs[0].default_value = val['transparent1'] #color
    transparent1.name = 'transparent1'

    emission1 = mat.node_tree.nodes.new('ShaderNodeEmission')
    emission1.location = (0,400)
    emission1.inputs[0].default_value = val['emission1'][0] #color
    emission1.inputs[1].default_value = val['emission1'][1]  #strength
    emission1.name = 'emission1'

    add = mat.node_tree.nodes.new('ShaderNodeMath')
    add.location = (200,800)
    add.operation = 'ADD'
    add.use_clamp = True
    add.name = 'add'

    #link layerweight to add
    mat.node_tree.links.new(layerweight.outputs[0], add.inputs[0])
    mat.node_tree.links.new(layerweight.outputs[1], add.inputs[1])

    wireframe = mat.node_tree.nodes.new('ShaderNodeWireframe')
    wireframe.location = (400,800)
    wireframe.use_pixel_size = True
    wireframe.inputs[0].default_value = val['wireframe'] #0.1 #size
    wireframe.name = 'wireframe'

    mixer1 = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer1.location = (400,600)
    mixer1.name = 'mixer1'
    
    #link mixer1 to add, transparent1 and emission1
    mat.node_tree.links.new(add.outputs[0], mixer1.inputs[0])
    mat.node_tree.links.new(transparent1.outputs[0], mixer1.inputs[1])
    mat.node_tree.links.new(emission1.outputs[0], mixer1.inputs[2])

    lightpath = mat.node_tree.nodes.new('ShaderNodeLightPath')
    lightpath.location = (600,1000)
    lightpath.name = 'lightpath'
    
    mixer2 = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer2.location = (600,600)
    mixer2.name = 'mixer2'

    #link mixer2 to wireframe, mixer1, emission1
    mat.node_tree.links.new(wireframe.outputs[0], mixer2.inputs[0])
    mat.node_tree.links.new(mixer1.outputs[0], mixer2.inputs[1])
    mat.node_tree.links.new(emission1.outputs[0], mixer2.inputs[2])

    transparent2 = mat.node_tree.nodes.new('ShaderNodeBsdfTransparent')
    transparent2.location = (600,400)
    transparent2.inputs[0].default_value = (1,1,1,1) #color
    transparent2.name = 'transparent2'

    mixer3 = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer3.location =  (800,800)
    mixer3.name = 'mixer3'
    
    volabsorp = mat.node_tree.nodes.new('ShaderNodeVolumeAbsorption')
    volabsorp.location = (800,600)
    volabsorp.inputs[0].default_value = val['volabsorp'][0] #color
    volabsorp.inputs[1].default_value = val['volabsorp'][1] #density
    volabsorp.name = 'volabsorp'

    emission2 = mat.node_tree.nodes.new('ShaderNodeEmission')
    emission2.location = (800,400)
    emission2.inputs[0].default_value = val['emission2'][0]
    emission2.inputs[1].default_value = val['emission2'][1]
    emission2.name = 'emission2'
    
    #link mixer 3 to lightpath, mixer2 and transparent2
    mat.node_tree.links.new(lightpath.outputs[9], mixer3.inputs[0]) #transparent depth
    mat.node_tree.links.new(mixer2.outputs[0], mixer3.inputs[1])
    mat.node_tree.links.new(transparent2.outputs[0], mixer3.inputs[2])

    mixer4 = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer4.location = (1000,400)
    mixer4.inputs[0].default_value = val['mixer4']
    mixer4.name = 'mixer4'
    
    #link mixer4 to volabsorp and emission2
    mat.node_tree.links.new(volabsorp.outputs[0], mixer4.inputs[1])
    mat.node_tree.links.new(emission2.outputs[0], mixer4.inputs[2])
   
    # ------- net transparency begin -------
    #mixer 5 controls transparency of object surface for appear/disappear annimation)
    mixer5 = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer5.location = (1000,600)
    mixer5.name = 'mixer5'

    # for appear/disappear annimation
    transparent3 = mat.node_tree.nodes.new('ShaderNodeBsdfTransparent')
    transparent3.location = (1000,800)
    transparent3.inputs[0].default_value = (1,1,1,1) #color
    transparent3.name = 'transparent3'

    #link mixer 5 to mixer 3 and transparent3 
    mat.node_tree.links.new(transparent3.outputs[0], mixer5.inputs[1])
    mat.node_tree.links.new(mixer3.outputs[0], mixer5.inputs[2])

    #mixer 6 controls transparency of object surface for appear/disappear annimation)
    mixer6 = mat.node_tree.nodes.new('ShaderNodeMixShader')
    mixer6.location = (1200,400)
    mixer6.name = 'mixer6'
    
    #link mixer4 to mixer6
    mat.node_tree.links.new(mixer4.outputs[0], mixer6.inputs[2]) # note mixer6.imputs[1] is empty

    value = mat.node_tree.nodes.new('ShaderNodeValue')
    mixer6.name = 'value'
    value.location = (1200,800)
    value.outputs[0].default_value = 1.0

    #link value to mixer 5 and 6
    mat.node_tree.links.new(value.outputs[0], mixer6.inputs[0])
    mat.node_tree.links.new(value.outputs[0], mixer5.inputs[0])


    # ------- net transparency end -------


    outmat = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
    outmat.location = (1400, 400)
    outmat.name = 'outmat'
    
    #link outmat to mixer5 and mixer6
    mat.node_tree.links.new(mixer5.outputs[0], outmat.inputs[0])
    mat.node_tree.links.new(mixer6.outputs[0], outmat.inputs[1])
    
    return mat



def create_lamp(name, strength):
    mat = bpy.data.materials.get(name)
    if mat is None:
      mat = bpy.data.materials.new(name)
      mat.use_nodes = True

      for node in mat.node_tree.nodes:
          mat.node_tree.nodes.remove(node)

      emission = mat.node_tree.nodes.new('ShaderNodeEmission')
      emission.location= (0,0)
      emission.inputs[0].default_value = (1,1,1,1)   # color light
      emission.inputs[1].default_value = strength  # strength

      outmat = mat.node_tree.nodes.new('ShaderNodeOutputMaterial')
      outmat.location = (200, 0)

      mat.node_tree.links.new(emission.outputs[0], outmat.inputs[0])

    return mat

def draw_rounded_cyl(loc,quat,height=1,radius=1):
    obs = []
    
    #create rounded cylinder height=2, rad=1
    bpy.ops.mesh.primitive_uv_sphere_add(size=1, location=(loc[0]+0.5,loc[1],loc[2]))
    obs.append(bpy.context.object)
    bpy.ops.mesh.primitive_uv_sphere_add(size=1, location=(loc[0]-0.5,loc[1],loc[2]))
    obs.append(bpy.context.object)
    bpy.ops.mesh.primitive_cylinder_add(radius=1,depth=1,location=loc,rotation=(0,math.pi/2,0))
    obs.append(bpy.context.object)

    # join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()

    # scale
    bpy.ops.transform.resize(value=(2./3/2*height,radius,radius))

    ob = bpy.context.object
    return ob

def draw_cyl(start=(0,0,0), end=(5,0,0), rcyl=1.0,end_fill_type='NOTHING'):
    
    r = [ end[i] - start[i] for i in range(3)]
    norm = math.sqrt(sum([r[i]**2 for i in range(3)]))

    dcyl = norm # depth cylinder

    cylshift = dcyl/2

    #obs = []
    vec0 = [1,0,0]
    # axis is centered on (0,0,0). No rotation initially
    #bpy.ops.mesh.primitive_uv_sphere_add(size=rcyl,location=(0,0,0))            #sphere at corner
    #objs.append(bpy.context.object)
    bpy.ops.mesh.primitive_cylinder_add(radius=rcyl,depth=dcyl,location=(cylshift,0,0),rotation=(0,math.pi/2.0,0),end_fill_type=end_fill_type)
    ob = bpy.context.object
    #obs.append(bpy.context.object)
    
    # join all objects
    #for ob in obs:
    #    ob.select=True
    #bpy.ops.object.join()

    # set rotation origin of axis
    bpy.context.scene.cursor_location = (0.0, 0.0, 0.0)  # set cursor to origin
    bpy.ops.object.mode_set(mode = 'OBJECT')             # make sure in object mode
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')      # set rotation origin to cursor location
    
    # get theta and u for rotation
    dot = sum([vec0[i]*r[i] for i in range(3)])/norm
    if dot > 1: dot = 1
    elif dot <-1: dot = -1

    if (dot == 1): #no rotation
        theta,u= 0,[1,0,0]
    elif (dot == -1): #pi rotation
        theta,u = math.pi,[0,1,0]
    else:
      u = [0]*3
      u[0] = vec0[1]*r[2] - vec0[2]*r[1]
      u[1] = vec0[2]*r[0] - vec0[0]*r[2]
      u[2] = vec0[0]*r[1] - vec0[1]*r[0]
      theta = math.acos(dot)
    
    # rotate
    ob = bpy.context.object
    #DO NOT SET ob.rotation_mode = 'QUATERNION'
    bpy.ops.transform.rotate(value=theta,axis=u)

    #translate
    bpy.ops.transform.translate(value=(start[0],start[1],start[2]))

    return ob


def draw_nucl_dna(loc,quat,nucl_bp_unwrap = 0,npoints=100,radius=46):

    obs = []
 
    #q = Quaternion((0,1,0,0))*Quaternion((quat[0],quat[1],quat[2],quat[3])) # rotate pi around x axis
    #quat = q
    #quat = q

    lengths_1kx5= {
          0: {'c': 55.8,      'd': 33.1,      'e': 67.1}, 
          9:  {'c': 48.178119, 'd': 36.932414, 'e': 90.105419} ,
          10: {'c': 48.017333, 'd': 35.357808, 'e': 91.871604} ,
          11: {'c': 47.175048, 'd': 32.370058, 'e': 91.850785} }

    c = lengths_1kx5[nucl_bp_unwrap]['c']  # distance of stem site from nucl center
    d = lengths_1kx5[nucl_bp_unwrap]['d']  # vertical distance between entering and exiting DNA sites
    e = lengths_1kx5[nucl_bp_unwrap]['e']  #real space distance between entering and exiting DNA sites
    f = math.sqrt(e*e - d*d); #horizontal distance between entering and exiting DNA
    g = math.sqrt(c*c - d*d/2.0/2.0)
    h = math.sqrt(g*g - f*f/2.0/2.0)
    gamma = 2*math.asin(0.5*f/g) # then rotate by asin around f to point u at dyad
    
    radtotal = 2*2*math.pi - gamma
    dtheta = radtotal/(npoints-1)
    theta0 = 0.5*gamma + 0.5*math.pi
    z0 = 0.5*d * radius / 50.
    dz = d/npoints * radius / 50.
    magnitude = 4 # how much I 'bulge' the DNA radius so that the viz looks nice
    
    if DNABONDS:
      cylrad =7. 
    else:
      cylrad = radius / 15. #pre 7/12/17
    
    #draw DNA
    for i in range(npoints):
      theta = theta0 + i*dtheta
      x = z0 - i*dz
      factor = magnitude*(math.cos(x/(2.*z0)*math.pi)) 
      z = math.sin(theta)*(radius + factor)
      y = math.cos(theta)*(radius + factor)

      bpy.ops.mesh.primitive_uv_sphere_add(size=cylrad, location=(x,y,z))
      obs.append(bpy.context.object)
      
      if i != 0:
        ob = draw_cyl((x,y,z),(xprev,yprev,zprev),rcyl = cylrad)
        obs.append(ob)

      xprev = x
      yprev = y
      zprev = z
    
    #draw ellipsoid
    #draw_rounded_cyl((0,0,0),(1,0,0,0),radius=radius,height=radius)
    #obs.append(bpy.context.object)

    #draw single sphere at center use use for rotations (wanted to use an 'empty' but that broke the join() command)
    bpy.ops.mesh.primitive_uv_sphere_add(size=0.01, location=(0,0,0))
    #bpy.ops.object.empty_add(type='PLAIN_AXES',location=(0,0,0))
    obs.append(bpy.context.object)

    # join all objects
    for ob in obs:
        ob.select=True
    #bpy.ops.object.modifier_apply(apply_as='DATA',modifier="Remesh")
    bpy.ops.object.join()

    #bpy.ops.object.parent_set(type="OBJECT",keep_transform=False)

    ob = bpy.context.object

    bpy.context.scene.cursor_location = (0.0, 0.0, 0.0)  # set cursor to origin
    bpy.ops.object.mode_set(mode = 'OBJECT')             # make sure in object mode
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')      # set rotation origin to cursor location

    #translate
    ob.location = loc[0],loc[1],loc[2]
    # rotate
    ob.rotation_mode = 'QUATERNION'
    ob.rotation_quaternion = Quaternion((quat[0], quat[1],quat[2],quat[3]))

    #theta,u = rotquat2tu(quat)
    #bpy.ops.transform.rotate(value=theta,axis=u)

    #translate
    #ob.select=True
    #bpy.ops.transform.translate(value=(loc[0],loc[1],loc[2]))
    #ob.select=False

    return ob


def draw_arrow(start=(0,0,0), end=(5,0,0), dcone = 2, cone_rfrac=1.5, rcyl=0.8):
    
    r = [ end[i] - start[i] for i in range(3)]
    norm = math.sqrt(sum([r[i]**2 for i in range(3)]))

    rcone = rcyl*cone_rfrac    # radius cone
    dcone = dcone     # depth cone
    dcyl = norm-dcone # depth cylinder

    cylshift = dcyl/2
    coneshift = dcyl+dcone/2

    obs = []
    vec0 = [1,0,0]
    # axis is centered on (0,0,0). No rotation initially
    #bpy.ops.mesh.primitive_uv_sphere_add(size=rcyl,location=(0,0,0))            #sphere at corner
    #objs.append(bpy.context.object)
    bpy.ops.mesh.primitive_cylinder_add(radius=rcyl,depth=dcyl,location=(cylshift,0,0),rotation=(0,math.pi/2.0,0))            #x-axis cyl
    obs.append(bpy.context.object)
    bpy.ops.mesh.primitive_cone_add(radius1=rcone,radius2=0,depth=dcone,location=(coneshift,0,0),rotation=(0,math.pi/2.0,0))  #x-axis cone
    obs.append(bpy.context.object)
    
    # join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()

    # set rotation origin of axis
    bpy.context.scene.cursor_location = (0.0, 0.0, 0.0)  # set cursor to origin
    bpy.ops.object.mode_set(mode = 'OBJECT')             # make sure in object mode
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')      # set rotation origin to cursor location
    
    # get theta and u for rotation
    dot = sum([vec0[i]*r[i] for i in range(3)])/norm
    if dot > 1: dot = 1
    elif dot <-1: dot = -1

    if (dot == 1): #no rotation
        theta,u= 0,[1,0,0]
    elif (dot == -1): #pi rotation
        theta,u = math.pi,[0,1,0]
    else:
      u = [0]*3
      u[0] = vec0[1]*r[2] - vec0[2]*r[1]
      u[1] = vec0[2]*r[0] - vec0[0]*r[2]
      u[2] = vec0[0]*r[1] - vec0[1]*r[0]
      theta = math.acos(dot)
    
    # rotate
    ob = bpy.context.object
    # DO NOT SET rotation_mode = 'QUATERNION'
    bpy.ops.transform.rotate(value=theta,axis=u)

    #translate
    bpy.ops.transform.translate(value=(start[0],start[1],start[2]))

    return ob


def draw_coordframe(loc,quat,length=(8,8,8),drawflag=(True,True,True),dcone=2,cone_rfrac=1.5,rcyl=0.8):
    obs = []
    for i,name in enumerate(['fvec','vvec','uvec']):
      if drawflag[i] != True: continue

      if name == 'fvec':   vec0 = Vector([1,0,0])
      elif name == 'vvec': vec0 = Vector([0,1,0])
      elif name == 'uvec': vec0 = Vector([0,0,1])
      #vec = quat_vec_rot(vec0,quat)
      #myend = [loc[j] + vec[j]*length[i] for j in range(3)]
      #ob = draw_arrow(start=loc, end=myend)
      myend = vec0*length[i]
      ob = draw_arrow(start=(0,0,0), end=myend,dcone=dcone,cone_rfrac=cone_rfrac,rcyl=rcyl)

      mat = create_material(name,COLORS[name])
      ob.data.materials.append(mat)
      obs.append(ob)

    # now join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()

    ob = bpy.context.object
    ob.rotation_mode = 'QUATERNION'
    ob.rotation_quaternion = Quaternion((quat[0], quat[1],quat[2],quat[3]))
    ob.location = loc

    ob.name='frame'

    return ob

def draw_dna_bead(loc,quat):
    obs = []
    
    #sphere
    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, size=3, location=loc)
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = 'bead'
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)    
    
    #coord frame
    ob = draw_coordframe(loc,quat,drawflag=(True,False,True),length=(6,6,6))
    obs.append(ob)    

    # now join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()
    
    ob.name = 'dna'

    return ob

def draw_dyad_bead(loc,quat):
    obs = []
    
    #sphere
    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, size=3, location=loc)
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = 'dyad'
    if NDNABONDS or DYADBONDS:
      mat = create_material(name,COLORS[name])
    else:
      mat = create_material_emission_glossy(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)    
    
    #coord frame
    ob = draw_coordframe(loc,quat,drawflag=(True,True,True),length=(8,8,8))
    obs.append(ob)    

    # now join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()
    
    ob.name = 'dyad'
    return ob


def draw_nucl(loc,quat):
    obs = []

    
    #sphere
    bpy.ops.mesh.primitive_uv_sphere_add(segments=32, size=3, location=loc)
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = 'zewdiecore'
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)    
 
    ob = draw_coordframe(loc,quat)
    #ob = draw_coordframe(loc,quat, length=(50,75,75),dcone=15,rcyl=5)
    obs.append(ob)    

    #nucl dna
    if DRAWNDNA:
      ob = draw_nucl_dna(loc,quat,nucl_bp_unwrap = 11)
      name = 'nucldna'
      mat = create_material(name,COLORS[name])
      ob.data.materials.append(mat)
      obs.append(ob)
    
    # zewdie
    if NDNABONDS or DYADBONDS:
      pass
    else:
      bpy.ops.import_mesh.stl(filepath=path_zewdie)
      bpy.ops.object.shade_smooth()
      ob = bpy.context.object
      name = 'zewdie'
      #mat = create_material(name,COLORS[name])
      mat = create_material_holograph()
      ob.data.materials.append(mat)
      ob.name = name
      ob.location = loc
      ob.rotation_mode = 'QUATERNION'
      ob.rotation_quaternion = Quaternion((quat[0], quat[1],quat[2],quat[3]))
      obs.append(ob)    


    # now join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()

    #bpy.ops.object.parent_set(type="OBJECT",keep_transform=False)

    ob = bpy.context.object
    ob.name = "nucl"


    return ob

def set_world_bg(color):
  # setup world background
  w = bpy.context.scene.world
  #w = bpy.data.worlds['World']
  w.use_nodes = True

  for node in w.node_tree.nodes:
    w.node_tree.nodes.remove(node)

  background = w.node_tree.nodes.new('ShaderNodeBackground')
  background.inputs[0].default_value = color

  outmat = w.node_tree.nodes.new('ShaderNodeOutputWorld')
  outmat.location = (200, 0)
  # ling together
  w.node_tree.links.new(background.outputs[0], outmat.inputs[0])



def setup_scene():
  
  # bottom plane
  #bpy.ops.mesh.primitive_plane_add(radius=100,location=(0,0,-10))
  #ob = bpy.context.object
  #name = 'plane'
  #mat = create_material(name,COLORS[name])
  #ob.data.materials.append(mat)

  obs = []
  
  # key light
  ###bpy.ops.object.lamp_add(type='SUN', radius=10, location=(0,-10,10))
  bpy.ops.mesh.primitive_plane_add(radius=3,location=(0,-10,5),rotation=(math.pi/3,0,0))
  ob = bpy.context.object
  mat = create_lamp('keylight',10)
  ob.data.materials.append(mat)
  obs.append(ob)

  # fill light
  bpy.ops.mesh.primitive_plane_add(radius=3,location=(10,0,5),rotation=(0,math.pi/3,0))
  ob = bpy.context.object
  mat = create_lamp('filllight',5)
  ob.data.materials.append(mat)
  obs.append(ob)
  ob.select = False

  # camera
  bpy.ops.object.camera_add(location=(10,-10,8),rotation=(1.06732,0,0.716928))
  ob = bpy.context.object
  obs.append(ob)
  
  for ob in obs:
    ob.select = True
  bpy.ops.object.parent_set(type='OBJECT',keep_transform=False)

def draw_bond(x1,x2,myname,rcyl=1*SCALEFACTOR):

    draw_cyl(start = x1, end=x2,rcyl=rcyl,end_fill_type='NOTHING')
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = myname
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)

    # set rotation origin of axis
    bpy.context.scene.cursor_location = x1               # set cursor to x1
    bpy.ops.object.mode_set(mode = 'OBJECT')             # make sure in object mode
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')      # set rotation origin to cursor location

    return ob 

# dna bonds are different, because they use the quaternion to show the twist of the DNA
# nameA - color name of primary cylinder connecting sites
# nameB - color name of secondary c connecting sites
# flength - how long along the fvec I go
def draw_dna_bond(x1,x2,q1,q2,nameA,nameB,rcylA=7*SCALEFACTOR,rcylB=3*SCALEFACTOR,flength=9*SCALEFACTOR):
    obs = []

    f0 = [1,0,0]
    f1 = quat_vec_rot (f0,q1)
    f2 = quat_vec_rot (f0,q2)
    x1B = x1 + f1*flength
    x2B = x2 + f2*flength

    bpy.ops.mesh.primitive_uv_sphere_add(size=rcylB, location=x1B)
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = nameB
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)

    bpy.ops.mesh.primitive_uv_sphere_add(size=rcylB, location=x2B)
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = nameB
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)

    draw_cyl(start = x1B, end=x2B,rcyl=rcylB,end_fill_type='NOTHING')
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = nameB
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)

    bpy.ops.mesh.primitive_uv_sphere_add(size=rcylA, location=x1)
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = nameA
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)


    draw_cyl(start = x1, end=x2,rcyl=rcylA,end_fill_type='NOTHING')
    bpy.ops.object.shade_smooth()
    ob = bpy.context.object
    name = nameA
    mat = create_material(name,COLORS[name])
    ob.data.materials.append(mat)
    obs.append(ob)




    # now join all objects
    for ob in obs:
        ob.select=True
    bpy.ops.object.join()

    #bpy.ops.object.parent_set(type="OBJECT",keep_transform=False)

    ob = bpy.context.object
    #ob.name = "bond"

    # set rotation origin of axis
    bpy.context.scene.cursor_location = x1               # set cursor to x1
    bpy.ops.object.mode_set(mode = 'OBJECT')             # make sure in object mode
    bpy.ops.object.origin_set(type='ORIGIN_CURSOR')      # set rotation origin to cursor location


    #ob.rotation_mode = 'QUATERNION'
    #ob.rotation_quaternion = q1

    return ob 

def import_traj(dump,frames):
    fps = 24

    # move all objects into layer 10 (should just be lamp, camera and block)
    obs = bpy.context.scene.objects
    for ob in obs:
        ob.layers[9] = True
        ob.layers[0] = False

    isfirstframe = True

    obs = []
    
    for ifrme in range(dump.nframes):

      dump.read_next_frame()
      if ifrme not in frames:
        continue

      #else load frame
      print("Loading frame ",ifrme)
      bpy.context.scene.frame_set(ifrme*fps)

      # create objects
      if isfirstframe: 
        for i in range(dump.natoms):
          x = (dump.x[i] - dump.xavg)*SCALEFACTOR
          q = dump.q[i]
          if (dump.type[i] == 1): # nucl
            ob = draw_nucl(x,q)
          elif (dump.type[i] == 2): # dna
            ob = draw_dna_bead(x,q)
          elif (dump.type[i] == 3): #dyad 
            ob = draw_dyad_bead(x,q)
          else:
              print("Error! Undefined atom type!")

          ob.scale = Vector((SCALEFACTOR,SCALEFACTOR,SCALEFACTOR))
          ob.name += str("-%05d" % i) #update name
          obs.append(ob) # append

        #draw bonds
        for i in range(dump.nbonds):
            mybond = dump.bonds[i]
            x1 = (mybond.x1 - dump.xavg)*SCALEFACTOR
            x2 = (mybond.x2 - dump.xavg)*SCALEFACTOR
            if mybond.drawmode == 'draw_bond':
              ob = draw_bond(x1, x2, mybond.typenameA)
            elif mybond.drawmode == 'draw_dna_bond':
              ob = draw_dna_bond(x1, x2, mybond.q1, mybond.q2, mybond.typenameA, mybond.typenameB)
            # dont scale bonds, otherwise they're too short, only scale cylinder radius
            #ob.scale = Vector((SCALEFACTOR,SCALEFACTOR,SCALEFACTOR))
            ob.name = str("bond-%05d" % i) #update name
            obs.append(ob) # append


      
      for ob in obs:
        ob.select = False

      # set keyframes
      for i in range(dump.natoms+dump.nbonds):
        name = obs[i].name   # get ob by name, in case order changes
        ob = bpy.context.scene.objects[name]
        

        if i < dump.natoms:
          ob.location = (dump.x[i] - dump.xavg) * SCALEFACTOR
          ob.keyframe_insert(data_path="location")
          ob.rotation_quaternion = Quaternion((dump.q[i][0], dump.q[i][1],dump.q[i][2],dump.q[i][3]))
          ob.keyframe_insert(data_path="rotation_quaternion")
        else: #is bond
          #pass
          #NEED TO FIX ME, NEED TO DEFINE CENTER OF EACH BOND AS X1+X2/2 AND USE THIS AS THE NEW BOND LOCATION
          mybond = dump.bonds[i-dump.natoms]
          x1 = (mybond.x1 - dump.xavg)*SCALEFACTOR
          #x2 = (mybond.x2 - dump.xavg)*SCALEFACTOR
          ob.location = x1
          ob.keyframe_insert(data_path="location")

          #This rotation doesn't do anything currently because rotation_mode='QUATERNION' isnt set
          # however if I set it, then the DNA spiraling location gets messed up
          ob.rotation_quaternion = mybond.q1 
          ob.keyframe_insert(data_path="rotation_quaternion")
                     


      isfirstframe = False
    

    # clean up
    bpy.context.scene.frame_end = dump.nframes*fps

    # switch interpolation type
    #old_type = bpy.context.area.type
    #bpy.context.area.type = 'GRAPH_EDITOR'
    #bpy.ops.graph.interpolation_type(type='LINEAR')
    ##bpy.ops.graph.interpolation_type(type='BEZIER')
    #bpy.context.area.type = old_type

    # deselect objects
    obs = bpy.context.scene.objects
    for ob in obs:
        ob.select=False
 

##################################
# START THE BUSINESS
##################################
#bpy.context.space_data.viewport_shade = 'MATERIAL'


bpy.context.scene.render.engine = 'CYCLES'

# this prevented many stacked transparent objects from turning black
bpy.context.scene.cycles.transparent_max_bounces= 24

#clear materials
items=bpy.data.materials.items()
for i in range(len(items)):
    bpy.data.materials.remove(items[i][1],do_unlink=True)

dump = DumpReader(path_traj)
frames = range(10)
import_traj(dump,frames)

set_world_bg(COLORS['worldbg'])

setup_scene()


#draw_nucl((0,0,0),(0,1,0,0))
#draw_nucl_dna((0,0,0),(0,1,0,0),nucl_bp_unwrap=11)
#draw_dna_bead((0,0,0),(0,1,0,0))
#draw_coordframe((0,0,10),(1,0,0,0))
#draw_arrow()
