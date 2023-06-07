
import json
import numpy as np

class config_options:
    """
    This class holds all of the configuration parameters for a classification
    simulation.  The components are all set to default values where possible in
    __init__.  The function parse_config() is used to parse a json formatted
    file which overwrites these defaults.
    """

    def __init__(self):
        # set some defaults.
        self.root = "./"
        """The location of the simulation root."""

        self.pass_num = 3
        """Number of iterations to run"""

        self.log_level = 0
        """Currently unused.  In the future will influence how much information
        is presented to the user"""



        self.dim_reduction_dims             = 100
        """Number of dimensions we will reduce to in the dimension reduction
        code."""
        self.dim_reduction_use_fft_avg      = True
        """Determines if we use the FFT-space averaging.  If false a simple
        average of the volumes ignoring the masks is used."""
        self.dim_reduction_avg_chunk_size   = 50
        self.dim_reduction_mat_chunk_size   = 150

        self.cluster_avg_chunk_size     = 50
        self.cluster_use_fft_avg        = True
        self.cluster_min_size_cutoff    = 0
        self.cluster_max_size_cutoff    = 0

        self.cluster_object_weighting   = 0
        self.cluster_kmeans_k           = 10
        self.cluster_kmeans_recursive_outlier_removal = 0

        self.given_templates = False
        self.default_templates_dir = ""
        self.template_align_mode = 2
        self.template_mask_ratio_cutoff = 0
        self.large_cluster_template_real_mask = 0

        self.L = 36

    def parse_config(self, opt_file):

        f = open(opt_file)
        conf = json.load(f)
        f.close()

        if "options" in conf:
            if 'pass_num' in conf['options']: self.pass_num = int(conf['options']['pass_num'])
            if 'log_level' in conf['options']: self.log_level = int(conf['options']['log_level'])
            if 'template_align_mode' in conf['options']: self.template_align_mode = int(conf['options']['template_align_mode'])

        if "dim_reduce" in conf:
            if 'dim_reduction_dims' in conf['dim_reduce']: self.dim_reduction_dims = int(conf['dim_reduce']['dim_reduction_dims'])
            if 'dim_reduction_avg_chunk_size' in conf['dim_reduce']: self.dim_reduction_avg_chunk_size = int(conf['dim_reduce']['dim_reduction_avg_chunk_size'])
            if 'dim_reduction_mat_chunk_size' in conf['dim_reduce']: self.dim_reduction_mat_chunk_size = int(conf['dim_reduce']['dim_reduction_mat_chunk_size'])
            if 'dim_reduction_use_fft_avg' in conf['dim_reduce']: self.dim_reduction_use_fft_avg = int(conf['dim_reduce']['dim_reduction_use_fft_avg'])

        if "cluster" in conf:
            if 'cluster_kmeans_k' in conf['cluster']: self.cluster_kmeans_k = int(conf['cluster']['cluster_kmeans_k'])
            if 'cluster_avg_chunk_size' in conf['cluster']: self.cluster_avg_chunk_size = int(conf['cluster']['cluster_avg_chunk_size'])
            if 'cluster_use_fft_avg' in conf['cluster']: self.cluster_use_fft_avg = int(conf['cluster']['cluster_use_fft_avg'])

        if "alignment" in conf:
            if 'L' in conf['alignment']: self.L = int(conf['alignment']['L'])

def parse_data(data_file):

    f = open(data_file)
    conf = json.load(f)
    f.close()

    subs = {}
    for record in conf:
        sub = {}

        if "subtomogram" not in record: raise
        assert record['subtomogram'] not in subs    # mxu: make sure that subtomogram id is unique

        sub['subtomogram' = record['subtomogram']

        if "mask" not in record: raise
        sub['mask'] = record['mask']

        if "angle" in record:
            if len(record["angle"]) != 3:
                raise
            sub['ang'] = np.array([float(_) for _ in record["angle"]])
        else:
            sub['ang'] = np.zeros(3, dtype=np.float)

        if "loc" in record:
            if len(record["loc"]) != 3:
                raise
            sub['loc'] = np.array([float(_) for _ in record["loc"]])
        else:
            sub['loc'] = np.zeros(3, dtype=np.float)
        
        subs[sub['subtomogram']] = sub

    return subs

