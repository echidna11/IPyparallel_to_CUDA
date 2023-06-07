

from work_queue.queue_worker import QueueWorker

import tomominer.average.util as avgu
import tomominer.align.util as aliu

class IterativeRefineWorker(QueueWorker):


    def __init__(self, host, port):
        QueueWorker.__init__(self, host, port) #, level=logging.WARN)

    
    def global_average_parallel(self, data, vol_shape, n_chunk, pass_dir, runner, use_fft_avg):
        return avgu.global_average_parallel(data, vol_shape, n_chunk, pass_dir, runner, use_fft_avg)    

    def vol_avg_fft_parallel_local(self, vol_shape, vmal_in, v_out_key, m_out_key):
        return avgu.vol_avg_fft_parallel_local(vol_shape, vmal_in, v_out_key, m_out_key)

    def vol_avg_fft_parallel_global(self, vol_shape, n_vol, vm_in, vol_avg_key, mask_avg_key):
        return avgu.vol_avg_fft_parallel_global(vol_shape, n_vol, vm_in, vol_avg_key, mask_avg_key)


    def segment_and_align_vols_against_one_template(self, vmal=None, inds=None, vm_tem=None, gauss_sigma=None, L=None):
        return aliu.segment_and_align_vols_against_one_template(vmal=vmal, inds=inds, vm_tem=vm_tem, gauss_sigma=gauss_sigma, L=L)
        

if __name__ == '__main__':

    import sys
    host = sys.argv[1]
    port = 5011

    worker = PursuitWorker(host, port)
    worker.run()
 

