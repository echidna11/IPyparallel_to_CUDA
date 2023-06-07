'''
use IMP to perform spherical packing
'''

import sys, copy, math
import numpy as N
import IMP as I
import IMP.core as IC
import IMP.container as ICO
import IMP.algebra as IAL
import IMP.atom as IAT
import tomominer.statistics.regression.linear as SRL

'''
main function of perform packing using IMP
'''
def do_packing(conf, op, pymol_file_name=None):       # parameters: conf: information of individual sphere, op: option

    m = I.Model()   # Create model
    box_size = [op['box']['x'], op['box']['y'], op['box']['z']]
    corner0 = IAL.Vector3D(0, 0, 0)
    corner1 = IAL.Vector3D(box_size[0], box_size[1], box_size[2])
    box = IAL.BoundingBox3D(corner0, corner1)   # Create Box

    ps = [None] * len(conf) # Create list of particle initialized to None
    # Initialize each particle at random position in box and with radius (bounding spheres)
    for i, b in enumerate(conf):
        p = I.Particle(m)   # Add particle
        v = IAL.get_random_vector_in(box)   # Get random position in box
        s = IAL.Sphere3D(v, b['radius'])    # Get radius of that particle (depending on what type of complex it is)
        IC.XYZR.setup_particle(p, s)   # Setup particle with coordinates and radius
        IAT.Mass.setup_particle(p, b['mass']) # Add Mass (required for molecular dynamics)
        IAT.LinearVelocity.setup_particle(p, IAL.Vector3D(0, 0, 0))
        ps[i] = p
    
    # Use ps[i].show() to see all the attributes of particle
    # create object for restraints
    r = I.RestraintSet(m)
    
    # Add first restriant to avoid spheres from overlapping each other
    r.add_restraint(IC.ExcludedVolumeRestraint(ps))
    
    '''
    #restricting the maximum pairwise distance between balls
    diameter = 5
    m.add_restraint(IC.DiameterRestraint(IC.HarmonicUpperBound(0,1), I.container.ListSingletonContainer(ps), diameter))

    '''
    # This code here was written by Min and bounds center of each sphere to be within the bounding box
    # But we need each sphere as whole within bounding box
    '''
    sscell = IC.BoundingBox3DSingletonScore(IC.HarmonicUpperBound(0, 1), box)
    rcell = ICO.SingletonsRestraint(sscell, ps)
    r.add_restraint(rcell)
    '''
    inds = {}
    ps_new = {}
    #for i in range(len(op['molecule']['types'])):
    for i in op['molecule']['types'].keys():
        inds[i] = [_['id'] for _ in conf if _['type']==i]
        ps_new[i] = [ps[_] for _ in inds[i]]
        boundary_dist = math.ceil(conf[inds[i][0]]['radius']) + 1
        sscell = IC.BoundingBox3DSingletonScore(IC.HarmonicUpperBound(boundary_dist, 1), box)
        rcell = ICO.SingletonsRestraint(sscell, ps_new[i])
        r.add_restraint(rcell)

    #----------------------------------------------------------
    # perform simulation

    temprature = op['temprature']
    sf = IC.RestraintsScoringFunction(r, "scoring function")
    o = IAT.MolecularDynamics(m)
    o.set_particles(ps)
    o.assign_velocities(temprature)
    o.set_scoring_function(sf)     # this is important to initialize the particles' movements, otherwise they will not move
    s = N.inf
    #cg=IC.ConjugateGradients(m)
    #s=cg.optimize(op['step'])
    recent_scores = []
    while (s > op['min score']) and (temprature > 0):
        #o.assign_velocities(temprature)
        md = IAT.VelocityScalingOptimizerState(m, ps, temprature)
        o.add_optimizer_state(md)
        s = o.optimize(op['step'])
        # count the number of particles inside the box.    IMPORTANT: watch this number at the steady state. If at steady state the number of balls inside the box is much smaller than total ball number, then the balls inside are supposed to be compact
        cx = N.array(get_center_coordinates(ps))
        is_in = N.array([True]*len(ps))
        for dim_i in range(3):
            is_in[cx[:, dim_i].flatten() < 0] = False
            is_in[cx[:, dim_i].flatten() >= box_size[dim_i]] = False
        

        print '\rscore:', s, '\ttemprature: ', temprature, '\tcenter inside box particles: ', is_in.sum(), '\t',
        if True:
            # use linear regression of recent scores to determine if we need to decrease the tempreature
            recent_scores.append(s)
            while len(recent_scores) > op['recent_scores number']:      recent_scores.pop(0)
            if len(recent_scores) == op['recent_scores number']:
                recent_scores_slope = SRL.one_var(N.array(range(len(recent_scores))), N.array(recent_scores) / recent_scores[0])[0]
                print 'slope', recent_scores_slope, '                  ',
                if N.abs(recent_scores_slope) < op['recent_scores slope min']:
                    temprature -= op['temprature decrease']
                    recent_scores = []
        sys.stdout.flush()

    #o.remove_optimizer_state(md)
    cx = get_center_coordinates(ps)

    if pymol_file_name is not None:        pymol_export(filename=pymol_file_name, model=m, particles=ps)
    # somehow I cannot process ps when this function returns. So let's export the balls with this function
    # add generated coordinate inforamtion back to conf
    conf = copy.deepcopy(conf)
    for i, b in enumerate(conf):    b['x'] = cx[i]
    return {'conf':conf, 'score':s, 'temprature':temprature, 'inside box num':is_in.sum()}


# export model in into a pym script that can be executed using pymol
def pymol_export(filename, model, particles):
    '''get the pym file'''
    pym2 = I.display.PymolWriter(filename)
    g3 = IC.XYZRsGeometry(ICO.ListSingletonContainer(model, particles))
    g3.set_name("beads")
    g3.set_color(I.display.Color(1,1,1))
    pym2.add_geometry(g3)



def get_center_coordinates(ps):
    x = []
    for p in ps:
        p = IC.XYZ(p)
        x.append([p.get_x(), p.get_y(), p.get_z()])
    
    return x

'''
# test code

import numpy as N
n = 1000
c = []
for i in range(n):   c.append(  {'r': N.random.random()}      )
len = 20
import tomominer.geometry.pack.imp_sphere as GPI
cx = GPI.do_packing(conf=c, op={'box':{'x':len, 'y':len, 'z':len}, 'temprature':1.0, 'temprature decrease':0.1, 'step':1, 'recent_scores number':10, 'recent_scores slope min':0.01, 'min score':0.1}, pymol_file_name='/tmp/t.pym')
'''

