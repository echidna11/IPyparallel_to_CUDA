import os

from work_queue.queue_worker import QueueWorker

import util


class PursuitWorker(QueueWorker):


    def __init__(self, host, port):
        self.cache_dir = os.getenv('CACHE_DIR');     assert(self.cache_dir is not None);    assert(os.path.isdir(self.cache_dir));
        QueueWorker.__init__(self, host, port) #, level=logging.WARN)
    

    def segment(self, vol_key, t, tm, t_seg, op):
        return util.segment(vol_key=vol_key, t=t, tm=tm, t_seg=t_seg, op=op, cache_dir=self.cache_dir)



if __name__ == '__main__':

    import sys
    host = sys.argv[1]
    port = 5011

    worker = PursuitWorker(host, port)
    worker.run()
 



