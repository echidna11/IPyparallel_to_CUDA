# use memory cached plus file io cache to cache files

import os
import stat
import shutil
import time
import uuid
import pickle

import numpy as N

import tomominer.io.file as IF

class Cache:
    def __init__(self, cache_dir=None,  tmp_dir=None, logger=None):

        self.logger = logger

        self.cache_dir = cache_dir
        #if self.cache_dir is not None:        assert  os.path.isdir(self.cache_dir)

        self.tmp_dir = tmp_dir
        #if self.tmp_dir is not None:        assert  os.path.isdir(self.tmp_dir)

    # randomly generate a tempory file name inside the temporary directory
    def get_temp_file_path(self, prefix=None, fn_id=None, suffix=None, ext=None):
        if fn_id is not None:
            fn = fn_id
        else:
            fn = str(uuid.uuid4())

        if prefix is not None:       fn = prefix + fn
        if suffix is not None:       fn = fn + suffix
        if ext is not None:       fn += ext

        assert  self.tmp_dir is not None
        
        return os.path.join(self.tmp_dir, fn)

    # generate temp file path (can according to task id given in fn_id), then check if the temp file exists. If not, save the data
    def save_tmp_data(self, d, fn_id=None):
        key = self.get_temp_file_path(prefix='tmp-', fn_id=fn_id, ext='.pickle')

        assert not os.path.isfile(key)      # this assertion is VERY IMPORTANT, it prevents a duplicated parallel job from returning values and supressing the successful return of the original job

        with open(key, 'wb') as f:   pickle.dump(d, f, protocol=-1)
    
        # change mod so that other users in the same group can read and delete
        os.chmod(key, 0666)

        return key

    
       


    # the file has to be small enough to fit into cache
    def get_mrc(self, path):

        path = str(path)

        return self.get_mrc_cache_fs(path)


    def get_mrc_cache_fs(self, path):
        return self.load_file_cache_fs(IF.read_mrc_vol, path)


    # WARNING: there can be cases that one program is copying a file while another starts to read it, you have to be careful to avoid this situation!!!
    def load_file_cache_fs(self, load_func, path):
        #assert(cache_dir is not None)

        path = os.path.abspath(path)
        if self.cache_dir is None:
            v = load_func(path)
            return v

        if not os.path.isdir(self.cache_dir):
            try:
                # if cache dir does not exist, try to create it
                os.makedirs(self.cache_dir)
            except:
                pass

            if not os.path.isdir(self.cache_dir):
                # if still a failure, quit
                raise OSError('cache_dir   ' + self.cache_dir)

        cache_path = self.cache_dir + path       # os.path.join does not concat two paths

        while self.load_file_cache_fs__is_miss(path=path, cache_path=cache_path):

            #print 'miss ' + cache_path

            # if the file is not in the cache, copy it from original location to the cache
            cache_path__dir = os.path.dirname(cache_path)

            if not os.path.isdir(cache_path__dir):
                try:
                    os.makedirs(cache_path__dir)
                except:
                    pass

            shutil.copyfile(path, cache_path)
            time.sleep(1)       # some file systems are slow, just wait for a little bit

        v = load_func(cache_path)

        return v


    def load_file_cache_fs__is_miss(self, path, cache_path):

        miss = False
        if not os.path.isfile(cache_path):
            miss = True
        else:
            # set miss if the cached file is older than original file
            path__st = os.stat(path)
            cache_path__st = os.stat(cache_path)

            if path__st.st_mtime > cache_path__st.st_mtime:
                miss = True

            if path__st.st_size != cache_path__st.st_size:
                miss = True

        return miss

