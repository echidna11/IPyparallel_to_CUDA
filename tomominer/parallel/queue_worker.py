


'''
worker program
dynamically import modules and run its functions

~/ln/tomominer/tomominer/parallel/queue_worker.py
'''


import traceback
import logging
import os
import sys
import time
import warnings
import uuid

from multiprocessing.pool import Pool

import importlib

from tomominer.parallel.RPCClient import RPCClient
from tomominer.parallel.RPCLoggingHandler import RPCLoggingHandler
from tomominer.parallel.Task import Task

from tomominer.io.cache import Cache

class QueueWorker:
    """
    Connect to a QueueServer running on a RPCServer using the RPCClient.  Run
    tasks retrieved from the QueueServer, returning results to the server.

    This instance exposes a set of functions from the tomograhy package tomo wrapped by swig.
    """


    def __init__(self, host=None, port=None, instance=None, pool=None, tmp_dir=None):

        self.worker_id = str(uuid.uuid4())
        self.work_queue = RPCClient(host,port)


        self.handler = RPCLoggingHandler(self.work_queue)
        self.logger = logging.getLogger()
        self.logger.addHandler(self.handler)
        self.logger.setLevel(logging.DEBUG)
        self.logger  = logging.LoggerAdapter(logging.getLogger(), {'host': os.environ.get('HOSTNAME', 'unknown'), 'job_id' : os.environ.get('PBS_JOBID', 'N/A').split('.')[0], 'source_type' : 'queue_worker'})

        if tmp_dir is None:     tmp_dir = os.getenv('TOMOMINER_TMP_DIR')
        assert      tmp_dir is not None

        self.cache = Cache(tmp_dir=tmp_dir,  logger=self.logger)
        self.cache_none = Cache(logger=self.logger)

        self.pool = pool



    def run(self, interval=5):

        while True:

            task = self.work_queue.get_broadcast_task(worker_id=self.worker_id)             # we first check if there is any broadcast task for this particular worker to process

            if not task:    task = self.work_queue.get_task(worker_id=self.worker_id)       # if no broadcast task, the worker will try to retrive a normal task

            if not task:
                time.sleep(interval)   # sleep a bit to prevent system to be too busy, block until a certain timeout
                continue

            #self.logger.debug("got task: %s", task)
            self.task = task        # we send the task id to the caller function, so that the result can be stored according to task id

            err, err_msg, result = self._dispatch()
            if err:
                self.logger.warning("task failed: %s, %s", repr(task), err_msg)        # mxu: this is broken, the error massage does not go to the queue server, so we ignore the logging
                continue

            
            while True:
                if self.work_queue.done_tasks_contains(task.task_id):    break       # stop if the task queue already masked the task as completion, this is used to reduce communication cost on sending back redundant results
                if self.work_queue.put_result(worker_id=self.worker_id, task_id=task.task_id, error=err, error_msg=err_msg, result=result):  break       # stop after succesfully sent to result queue
                time.sleep(10)       # keep trying to put result to queue server until success

    def _dispatch(self):
        if self.task.method.startswith('_'):
            return (True, "method starts with '_'", None)
        else:
            #func = getattr(__import__(self.task.module), self.task.method)
            assert      self.task.module is not None
            try:
                modu = importlib.import_module(self.task.module)
            except Exception:
                ex_type, ex, tb = sys.exc_info()
                return (True, 'loading module error: %s  in sys path  %s     ;    exception  %s'%(self.task.module, repr(sys.path), repr(traceback.format_tb(tb))), None)

            try:
                func = getattr(modu, self.task.method)
            except:
                return (True, 'method not found: %s '%(self.task.method), None)

            if not callable(func):
                return (True, "method not callable", None)

            try:
                assert 'self' not in self.task.kwargs
                # For MPP (next two lines)
                #self.task.kwargs['self'] = self
                #result = func(*self.task.args, **self.task.kwargs)
                # For my writen functions (next one line)
                result = func(**self.task.kwargs)
                
                return (False, None, result)

            except Exception as ex:
                return (True, traceback.format_exc(), None)

                sys.stderr.write(traceback.format_exc())
                self.logger.error("Exception: %s", traceback.format_exc())      # mxu: this cannot be reached?

if __name__ == '__main__':

    warnings.filterwarnings('error')
 
    host = sys.argv[1]
    port = 5011

    if True:
        worker = QueueWorker(host=host, port=port)          # pool must be initialized at beginning of script to prevent memory blow up
    else:
        worker = QueueWorker(host=host, port=port, pool=Pool(processes=2))          # pool must be initialized at beginning of script to prevent memory blow up

    worker.run()
 
