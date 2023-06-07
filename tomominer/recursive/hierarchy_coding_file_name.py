
# walk through a hierarchy of recursive classification output averaging files, and rename the files by putting the hierarchy code in front of the file name


def extract_hierarchy_code(d):
    c_ss = d.split('clus_')

    hc_s = ''
    hc_i = []
    for c_s in c_ss:
        c_s_sp = c_s.split('/')
        try:
            clus_i = int(c_s_sp[0])
        except ValueError:
            continue   
        hc_i.append(clus_i)
        hc_s += ( repr(clus_i) + '-')

    return (hc_s, hc_i) 


import os
import sys
import shutil

if __name__ == '__main__':

    src_dir = sys.argv[1]
    des_dir = sys.argv[2]


    for root, subFolders, files in os.walk(src_dir):
        (hc_s, hc_i) = extract_hierarchy_code(root)
        if len(hc_i) == 0:        continue
        #print hc + '\t' + root    


        for f in files:
            clus_code = f.split('_')
            if len(clus_code) > 2:
                clus_code = repr(int(clus_code[2]))
            else:
                clus_code = ''

            f_ext = f.split('.');       f_ext = f_ext[len(f_ext)-1]

            shutil.copyfile(os.path.join(root, f), os.path.join(des_dir, hc_s  + clus_code + '.' + f_ext))




