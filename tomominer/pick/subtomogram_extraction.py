#!/usr/bin/env python


# functions for extracting subtomograms
import os, uuid
import tomominer.io.path_util as IP
import numpy as N
import tomominer.image.vol.util as CV
import tomominer.image.vol.wedge.util as W
import tomominer.model.util as MU
import tomominer.io.file as IOF



# given a large tomogram v, extract subtomograms at corresponding locations x, and save extracted subtomograms inside a dir
# return a vector indicating whether individual subtomograms can be successfully extracted
# parameters: v: tomogram,    x: peak locations,     siz: subtomogram size to be cut,      out_dir: output dir,    r: radius of spherical mask in real space,     ids: id of subtomograms
def subtomogram_extraction(v, x, siz, out_dir, r=None, ids=None, max_file_num_per_dir=100):

    print 'subtomogram_extraction()', 'extracting', len(x), 'subtomograms'

    if ids is None:     ids = list(range(len(x)))

    d = CV.grid_distance_to_center( CV.grid_displacement_to_center(siz) )
    print "d calculated"
    file_paths = {}
    for i in range(len(x)):

        f_dir = os.path.join(out_dir, IP.id_digit_seperation(i=ids[i], max_item_per_dir=max_file_num_per_dir))
        out_file = os.path.join( f_dir, '%06d.mrc'%(ids[i]) )

        if os.path.isfile(out_file):
            file_paths[i] = out_file
            continue

        print '\r', i, '\t', float(i) / len(x), '              ',
        vc = CV.cut_from_whole_map(whole_map=v, c=x[i,:].flatten(), siz=siz)
        if vc is None:      continue

        if (r is not None) and (r[i] > 0):
            #raise Exception('currently use a hard mask, which may introduce boundary effect (as can be seen by visual inspection) that bias the alignment. Need to replace with soft margin mask. Also, the isolation of individual complex regions will result in a simulation that is similiar to the simulation of individual complexes, right?')
            vc[d > r[i]] = vc[d <= r[i]].mean()           # put a spherical mask according to radius given in r


        vc = N.array(vc, dtype=float, order='F')

        if not os.path.isdir(f_dir):
            try:
                os.makedirs(f_dir)
            except:
                pass

        file_paths[i] = out_file
        IOF.put_mrc(vc, out_file)


    return file_paths





def main():

    import json
    with open('subtomogram_extraction__config.json') as f:      op = json.load(f)
    op['out_dir'] = os.path.abspath(op['out_dir'])

    if 'selected_tomogram_ids' in op:
        print 'selected_tomogram_ids', op['selected_tomogram_ids']
        selected_tomogram_ids = set(op['selected_tomogram_ids'])
    else:
        selected_tomogram_ids = None


    with open('particle_picking_dog__config.json') as f:     pp_op = json.load(f)
    with open(pp_op['tomogram_info_file']) as f:        tom_info = json.load(f)
    with open(op['peak_file']) as f:     pp = json.load(f)

    # My Code added to extract only top few peaks
    find_maxima = pp_op['find_maxima']
    if find_maxima:
        pp = sorted(pp, key=lambda _:(-_['val']))      # sort pp according to increase of peak value because our peaks are minimas
    else:
        pp = sorted(pp, key=lambda _:(_['val']))
    
    if 'top_num' in op:     pp = pp[:op['top_num']]
    assert      len(pp) > 0
    # My modification ends

    print 'initially loaded', len(pp), 'peaks'

    import csv
    def load_num_list(fname):
        l = []
        with open(fname) as f:
            reader = csv.reader(f)
            for row in reader:

                if row is None:     continue
                if len(row) == 0:   continue

                row = row[0]

                if row is None:     continue
                
                row = row.strip()
                if len(row) == 0:   continue

                try:
                    l.append(float(row))
                except:
                    print 'error parsing', row
        return l

    def make_wedge_file__from_angle(tilt_angle_file, out_file):
        print '********* make sure the generated wedge mask has correct tilt angle!!! ***************', 
        tilt_angs = load_num_list(tilt_angle_file)
        tilt_ang = N.max(N.abs(tilt_angs))
        mask = W.wedge_mask(size=size, ang1=90-tilt_ang, verbose=True)            # important!! the missing wedge angle is 90 - tilt angle range!!!!
        IOF.put_mrc(mask, out_file)

    def make_wedge_file__from_info(wedge, out_file):
        print '********* make sure the generated wedge mask has correct tilt angle!!! ***************',
        mask = W.wedge_mask(size=size, ang1=wedge['ang1'], ang2=wedge['ang2'], direction=wedge['direction'], verbose=True)
        IOF.put_mrc(mask, out_file)


    def make_wedge_file__from_average(subtomogram_paths, out_file):
        # calculate the average, then use the fft magnitude to generate mask
        vs = None
        for i in subtomogram_paths:
            if not os.path.isfile(subtomogram_paths[i]):        continue
            v = IOF.read_mrc(subtomogram_paths[i])['value']
            if vs is None:      vs = N.zeros(v.shape, dtype=N.float)
            vs += v

        vs /= len(subtomogram_paths)
        vs = N.fft.fftshift(    N.fft.fftn(vs)  )
        vs = N.abs(vs)

        mask = N.zeros(vs.shape, dtype=N.float)
        mask[vs >= N.median(vs)] = 1.0
        IOF.put_mrc(mask, out_file)

        IOF.put_mrc(vs, os.path.join('/tmp', str(uuid.uuid1())+'.mrc'))


        

    out_json = {}
    for tom in tom_info:
        if (selected_tomogram_ids is not None) and (tom['id'] not in selected_tomogram_ids): continue

        pp_t = [_ for _ in pp if (_['tomogram_id']==tom['id']) ]
        if len(pp_t) == 0:            continue

        out_rec = {}

        print 'loading', tom['vol_file']
        v = IOF.read_mrc(tom['vol_file'], show_progress=True)['value']

        if 'tomogram_normalize' in op:
            # randomly sample a number of voxels, then use these sampled values to normalize the tomogram. The reason to random sample is to avoid effect of golden particles, which are outliers but just occupy very small region, and thus should be hardly included into the sample
            print 'normalizing...',
            v_s = N.random.choice(v.flatten(), size=op['tomogram_normalize']['sample_size'])
            v_s__mean = float(v_s.mean())           ;       print 'mean:', v_s__mean,
            v_s__std = float(v_s.std())             ;       print 'std:', v_s__std

            v = v.astype(N.float32)
            v = (v - v_s__mean) / v_s__std

            out_rec['sample'] = {}
            out_rec['sample']['mean'] = v_s__mean
            out_rec['sample']['std'] = v_s__std

            

        locs = N.array(   [_['x'] for _ in pp_t]   )
        vals = [_['val'] for _ in pp_t]
        ids = [_['id'] for _ in pp_t]

        for size_ratio in op['size_ratios']:

            out_dir = os.path.join(op['out_dir'], str(pp_op['sigma1']), str(size_ratio), str(tom['id']))
            if not os.path.isdir(out_dir):      os.makedirs(out_dir)

            sigma1 = pp_op['sigma1'] / tom['voxel_spacing'] 
            size = N.array(         [  int(N.ceil(2.0 * sigma1 * size_ratio))  ] * 3         )
            print 'generating subtomograms of size', size

            subtomogram_paths = subtomogram_extraction(v=v, x=locs, ids=ids, siz=size, out_dir=out_dir)

            wedge_mask_file = os.path.join(out_dir, 'wedge-mask.mrc')
            if 'tilt_angle_file' in tom:
                assert      'wedge' not in tom
                make_wedge_file__from_angle(tilt_angle_file=tom['tilt_angle_file'], out_file=wedge_mask_file)
            elif 'wedge' in tom:
                assert  'tilt_angle_file' not in tom
                make_wedge_file__from_info(wedge=tom['wedge'], out_file=wedge_mask_file)
            else:
                make_wedge_file__from_average(subtomogram_paths=subtomogram_paths, out_file=wedge_mask_file)
                raise Exception('not wedge info provided')
            
            data_json = []
            for sp_i in subtomogram_paths:
                r = {}
                r['subtomogram'] = subtomogram_paths[sp_i]
                r['mask'] = wedge_mask_file
                r['uuid'] = 'st--' + str(uuid.uuid4())
                r['peak'] = {}
                r['peak']['id'] = pp_t[sp_i]['id']          # this id corresponds to peak id
                r['peak']['loc'] = pp_t[sp_i]['x']
                r['peak']['val'] = pp_t[sp_i]['val']
                r['tomogram_id'] = tom['id']


                data_json.append(r)

            info_file = os.path.join(out_dir, 'info.json')
            with open(info_file, 'w') as f:     json.dump(data_json, f, indent=2)

            if pp_op['sigma1'] not in out_rec:      out_rec[pp_op['sigma1']] = {}
            out_rec[pp_op['sigma1']][size_ratio] = {}
            out_rec[pp_op['sigma1']][size_ratio]['info_file'] = info_file

        out_json[tom['id']] = out_rec

    with open('subtomogram_extraction__out.json', 'w') as f:        json.dump(out_json, f, indent=2)



if __name__ == "__main__":
    main()

