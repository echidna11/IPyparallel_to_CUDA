

# a wrapper to call Zach's sphere packing method


import StringIO



# generate a conf file according to given options,
# see /home/rcf-47/mxu/ln/electron/util/geometry/packing/sphere/potential-brownian/pack_c/conf.README

def generate_conf_file(op):

    out = StringIO.StringIO()

    print >>out, op['file']['ftype'], op['file']['version']          # file type and version (ftype, version)

    if 'fixed' not in op['molecule']:   op['molecule']['fixed'] = []
    print >>out, len( op['molecule']['fixed'] )              # number of molecules with fixed positions. (n_mol)

    for r in op['molecule']['fixed']:       print >>out, r['id'], r['x'], r['y'], r['z'], r['fixed']              # id x y z fixed  (ids have to start at 1, fixed if it should not move during equilibration)

    if 'random' not in op['molecule']:   op['molecule']['random'] = []
    print >>out, len( op['molecule']['random'] )            # number of random records to read

    for r in op['molecule']['random']:       print >>out, r['type'], r['number'], r['ignore_overlap']              # type, number of random particles, if we ignore overlaps while placing particles.

    re = out.getvalue()

    out.close()

    return re



def generate_param_file(op):

    out = StringIO.StringIO()

    print >>out, op['file']['ftype'], op['file']['version']                                     # [file type] [file version] (ftype, version) 
    print >>out, op['n_max']                                                    # max number of particles we will ever have in the system. (n_max)
    print >>out, op['dt']                                                # time step (ns)    (dt)
    print >>out, op['temp']                                              # temperature (K)   (temp)
    print >>out, op['sim_time']                                          # total simulated time (ns) (sim_time)
    print >>out, op['force_eq_time']                                    # force equilibration time. (ns) (force_eq_time)
    print >>out, op['eq_time']                                          # equilibration time. (ns)      (eq_time)
    print >>out, op['log_freq']                                                 # frequency of log dumps. (ns)  (log_freq)
    print >>out, op['vol_high']['x'], op['vol_high']['y'], op['vol_high']['z']                                                  # range of volume (nm nm nm)    (vol.high)
    print >>out, op['bc']['x'], op['bc']['y'], op['bc']['z']                                                 # boundary conditions (x,y,z) (0 = periodic, 1 = reflecting, 2 = absorbing ) (bc)
    print >>out, op['RCUT']                                                 # rcut, maximum force ranges.  determines box size. (nm)    (RCUT)
    print >>out, op['molecule']['version']                                                 # version of molecule types section. (version)
    print >>out, len(op['molecule']['types'])                                                 # number of molecule types coming (number of lines to parse as molecule descriptions.) (n_types)

    for r in op['molecule']['types']:      print >>out, r['id'], r['mass'], r['radius'], r['diffusion_rate']                                                 # [type id] [mass (kDa)] [radius (nm)] [diffusion rate (nm^2/ns)] (type id must be sequential from 1.)


    re = out.getvalue()
    
    out.close()
    
    return re
    

def parse_result(f):

    running_info = []

    coordinates = []

    while True:
        line = f.readline()

        if line is None: break

        if line == '':      break

        line = line.strip()

        if line.startswith('INFO:'):
            running_info.append(line)
            continue

        ss = line.split(' ')
        if len(ss) != 5:
            print 'failed parsing:', line
            continue

        coordinates.append( {'id':int(ss[0]), 'type':int(ss[1]), 'x':[float(ss[2]), float(ss[3]), float(ss[4])]} )


    return {'running_info':running_info, 'coordinates':coordinates}



import os
import tempfile
import subprocess

# prepare input files, call Zach's program, and return packing results
def do_packing(op):

    tmp_f, conf_file = tempfile.mkstemp(prefix='tmp-conf-')
    os.close(tmp_f)
    with open(conf_file, 'wb') as f:        f.write(generate_conf_file(op['conf']))

    tmp_f, param_file = tempfile.mkstemp(prefix='tmp-param-')
    os.close(tmp_f)
    with open(param_file, 'wb') as f:        f.write(generate_param_file(op['param']))


    proc = subprocess.Popen([op['program'], param_file, conf_file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    re = parse_result(proc.stdout)

    #out, err = proc.communicate()
    #assert      proc.returncode == 0


    return re



        




