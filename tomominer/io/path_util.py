
# functions for manipulating file paths

import os


# walk through every subdir, find and record last passes of a iterative process
def last_pass_dirs(scan_dir=os.getcwd(), selected_pass_i=-1):

    selected_folders = []
    for root, sub_folders, files in os.walk(scan_dir):

        max_t = {}
        max_t['pass_i'] = -1
        for sub_folder in sub_folders:
            sub_folder_path = os.path.join(root, sub_folder)

            pred_member_file = os.path.join(sub_folder_path, 'data_config.json')
            if not os.path.exists(pred_member_file):    continue
    
            cluster_average_select_file = os.path.join(sub_folder_path, 'cluster_average_select.pickle')
            if not os.path.exists(cluster_average_select_file):    continue
            
            cluster_info_list_file = os.path.join(sub_folder_path, 'cluster_info_list.json')
            if not os.path.exists(cluster_info_list_file):                cluster_info_list_file = None

            cluster_info_stat_file = os.path.join(sub_folder_path, 'cluster_info_stat.pickle')
            if not os.path.exists(cluster_info_stat_file):                cluster_info_stat_file = None

            cluster_info_file = os.path.join(sub_folder_path, 'cluster_info.pickle')
            if not os.path.exists(cluster_info_file):                cluster_info_file = None


            pass_i = int(sub_folder.split('_')[1])

            if (selected_pass_i >= 0) and (pass_i != selected_pass_i):  continue

            if pass_i > max_t['pass_i']:
                max_t['pass_i'] = pass_i
                max_t['subfolder_path'] = sub_folder_path
                max_t['root_path'] = root
                max_t['pred_member_file'] = pred_member_file
                max_t['cluster_average_select_file'] = cluster_average_select_file

                if cluster_info_list_file is not None:      max_t['cluster_info_list_file'] = cluster_info_list_file
                if cluster_info_stat_file is not None:      max_t['cluster_info_stat_file'] = cluster_info_stat_file
                if cluster_info_file is not None:      max_t['cluster_info_file'] = cluster_info_file


        if max_t['pass_i'] >= 0:        selected_folders.append(max_t)


    return selected_folders




# given an id, decompose it to digits to make subfolders, so that each subfolder contain only a limited number of files. This is used to control the maximum file per dir
# paramters:        i: id
def id_digit_seperation(i, max_item_per_dir):

    dig = []
    it = int(i / max_item_per_dir)
    while it > 0:
        dig.append(     str(it % 10)     )
        it = int(it / 10)
    dig.reverse()

    f_dir = ''
    for d in dig:
        f_dir = os.path.join(f_dir, d)

    if len(f_dir) == 0:     f_dir = '0'


    return f_dir



