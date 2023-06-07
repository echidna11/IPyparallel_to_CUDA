

'''

source ~/proj/structure/modelling/util/imp/2.2.0/imp.bashrc

'''

import time
import numpy as N

import IMP as I
import IMP.core as IC
import IMP.container as ICO
import IMP.atom as IA



use_LennardJones = False

#I.set_log_level(I.SILENT)
m = I.kernel.Model()


cube_half = 100
corner1=I.algebra.Vector3D(-cube_half, -cube_half, -cube_half)
corner2=I.algebra.Vector3D(cube_half, cube_half, cube_half)
box=I.algebra.BoundingBox3D(corner1,corner2)


n_ball = 10000

'''
xyzr = IC.create_xyzr_particles(m, n_ball, 1.0)
ball_set = ICO.ListSingletonContainer(xyzr)
'''

r = list(range(n_ball));        r = [N.random.random()*10 + 1 for _ in r]

#ps = ICO.ListSingletonContainer(m)
ps = []
for i in range(n_ball):
    p = I.Particle(m)
    v = I.algebra.get_random_vector_in(box)
    s = I.algebra.Sphere3D(v, r[i])
    pt = IC.XYZR.setup_particle(p, s)
    if use_LennardJones:        IA.LennardJones.setup_particle(p, cube_half*2)      # LJ potential used by IA.LennardJonesPairScore()
    #IA.Diffusion.setup_particle(p, 1)                   # needed by Brownian dynamics simulation
    pt = IA.Mass.setup_particle(p,1)
    p.add_attribute(I.kernel.FloatKey('vx'), 0.0)       # initial velocity attribute is needed for molecular dynamics simulation
    p.add_attribute(I.kernel.FloatKey('vy'), 0.0)       # initial velocity attribute is needed for molecular dynamics simulation
    p.add_attribute(I.kernel.FloatKey('vz'), 0.0)       # initial velocity attribute is needed for molecular dynamics simulation
    ps.append(p)
    #ps.add_particle(p)



m.add_restraint(I.core.ExcludedVolumeRestraint(ps))         # avoid the balls from overlaping each other



'''
# restricting the maximum pairwise distance between balls

diameter = 5
m.add_restraint(IC.DiameterRestraint(IC.HarmonicUpperBound(0,1), I.container.ListSingletonContainer(ps), diameter))
'''

if True:
    # restricting the location of balls inside a sphere with radius space_radius
    # Note: such restraint cannot work with Brownian dynamics
    space_radius = cube_half 
    center = I.algebra.Vector3D(0,0,0)
    ubcell = IC.HarmonicUpperBound(space_radius, 1.0)
    sscell = IC.DistanceToSingletonScore(ubcell, center)
    rcell = I.container.SingletonsRestraint(sscell, ps)
    m.add_restraint(rcell) #2

if False:
    harminic_bound = IC.Harmonic(0, cube_half)
    sscell = IC.DistanceToSingletonScore(harminic_bound, I.algebra.Vector3D(0,0,0))
    for p in ps:    m.add_restraint(IC.SingletonRestraint(sscell, p))



if use_LennardJones:
    # adding LennardJones attraction, so that the nearby particles will attract each other to form clusters of particles
    # NOTE: this is slow when there are large number of particles
    nbl = ICO.ClosePairContainer(ICO.ListSingletonContainer(ps), cube_half*2)
    sf = IA.ForceSwitch(cube_half, cube_half*2)
    pair_score = IA.LennardJonesPairScore(sf)
    m.add_restraint(ICO.PairsRestraint(pair_score, nbl))




'''
# display coordinates
for p in ps:    print p.show()
'''


def dsp_xyz(ps):
    for p in ps:
        pt = IC.XYZ(p)
        print pt.get_x(), pt.get_y(), pt.get_z()


dsp_xyz(ps)


temprature = 1000
o = IA.MolecularDynamics(m)
o.set_particles(ps)
o.assign_velocities(temprature)     # this is important to initialize the particles' movements, otherwise they will not move
md = IA.VelocityScalingOptimizerState(m, ps, temprature)  # replace 300 K with 500 K
o.add_optimizer_state(md)

step = 1000
cur_time = time.time()
s=o.optimize(step)
#dsp_xyz(ps)
print 'time used', time.time() - cur_time, 'sec'

#o.remove_optimizer_state(md)


# export model in into a pym script that can be executed using pymol
def pymol_export(filename, particles):
    '''get the pym file'''
    pym2 = I.display.PymolWriter(filename)
    g3 = I.core.XYZRsGeometry(ICO.ListSingletonContainer(particles))
    g3.set_name("beads")
    g3.set_color(I.display.Color(1,1,1))
    pym2.add_geometry(g3)


pymol_export('/home/rcf-47/mxu/tmp-staging/tmp.pym', ps)



'''
# to display, run following

scp hpc-cmb.usc.edu:/home/rcf-47/mxu/tmp-staging/tmp.pym /tmp/tmp.pym ;  pymol /tmp/tmp.pym

'''


'''
# BrownianDynamics simulation

bd = IA.BrownianDynamics(m)
bd.set_particles(ps)
bd.set_time_step(10000)
step = 1000
bd.optimize(step)

'''

'''
# ConjugateGradients simulation
step = 1000
o = IC.ConjugateGradients(m)
o.set_particles(ps)
f=o.optimize(step)
'''



'''

ball_set = ICO.ListSingletonContainer(ps)

for i in range(n_ball):
    p = ball_set.get_particle(i)
    v = I.algebra.get_random_vector_in(box)
    s = I.algebra.Sphere3D(v, r[i])
    IC.XYZR.setup_particle(p, s)
    IA.LennardJones.setup_particle(p, 0.5)
    IA.Mass.setup_particle(p,1)

m.add_restraint(I.core.ExcludedVolumeRestraint(ball_set))    # to avoid balls to overlap

'''


'''
for i in range(n_ball):
 
    #v = I.algebra.Vector3D(1.0, 2.0, 3.0)
    #p.set_coordinates(coor)
    #p.set_radius(r)
    v = I.algebra.get_random_vector_in(box)
    p.set_coordinates(v)
    p.set_radius(r[i])
    #s = I.algebra.Sphere3D(v, r[i])
    #IC.XYZR.setup_particle(p, s)
    IA.LennardJones.setup_particle(p, 0.5)
    IA.Mass.setup_particle(p,1)

    #p0 = ball_set.get_particle(i)
    #p = I.core.XYZR(p0)
    #p = I.kernel.Particle(m)
    #d = I.core.XYZR.setup_particle(p)
    #d.set_coordinates(I.algebra.Vector3D(10.0, 10.0, 10.0))
'''








'''
adding restraint and energy terms, see following examples
https://www.google.com/search?client=ubuntu&channel=fs&q=imp+add_restraint&ie=utf-8&oe=utf-8

http://www.salilab.org/~drussel/unstructured/imp_doc/core_examples.html
http://www.salilab.org/~drussel/unstructured/imp_doc/atom_examples.html

http://svn.salilab.org/imp/branches/stable/modules/atom/test/test_lennard_jones_decorator.py

'''


