import time
import os
import sys
import uuid
import pickle
import thread

import numpy as N

from tomominer.parallel.RPCClient import RPCClient
from tomominer.parallel.RPCLoggingHandler import RPCLoggingHandler
from tomominer.parallel.Task import Task

import logging


class QueueMaster:
    """
    A handler for submitting a number of jobs to a Work Queue implemented on
    top of RPCServer, and blocking until all complete.
    """

    def __init__(self, host, port):
        """

        Connect to the server.  Setup a logging handler so all logging events
        are sent to the server.  A new project is created for the duration of
        this objects lifecycle.

        :param host
        :param port
        """
        self.work_queue = RPCClient(host,port)

        self.handler = RPCLoggingHandler(self.work_queue)
        self.logger = logging.getLogger()
        self.logger.addHandler(self.handler)
        self.logger.setLevel(logging.DEBUG)
        self.logger  = logging.LoggerAdapter(logging.getLogger(), {'host': os.environ.get('HOSTNAME', 'unknown'), 'job_id' : os.environ.get('PBS_JOBID', 'N/A').split('.')[0], 'source_type' : 'queue_master'})


        self.proj_id = str(uuid.uuid4())
        self.work_queue.new_project(self.proj_id)       # IMPORTANT: the keep alive signal is sent using an independent thread, therefore must create an independent RPCClient for an independent channel, in order to avoid conflication with the main channel

        thread.start_new_thread(QueueMaster.keep_alive, (self, RPCClient(host,port)))


    # TODO: this has been deprecated, now that we understand pickle better.  Delete from future version.
    def task(self, priority=1000, module=None, method=None, args=[], kwargs={}):
        return Task(priority=priority, proj_id=self.proj_id, module=module, method=method, args=args, kwargs=kwargs)

    def __del__(self):
        self.work_queue.del_project(self.proj_id)      # ignore this step, the project will be automatically deleted when at queue server side according to time out
        return


    def keep_alive(self, work_queue_alive, interval=10):

        while True:
            work_queue_alive.project_keep_alive(self.proj_id)
            time.sleep(interval)



    def run(self, tasks, max_time=None, max_retry=1000, one_at_a_time=False):       # mxu: change max_retry=1 to max_retry=5
        """
        Run a number of tasks by submitting to the QueueServer for processing by QueueWorker processes.

        :param tasks The list of Task() objects to submit.
        :param max_time Currently unused.
        :param max_retry Number of times to attempt the function call in the event of an error being returned, or a timeout occuring.
        """

        task_dict = {}
        """ A mapping from task_id to Task objects. """

        state     = {}
        """ """

        results   = {}
        """ """

        for t in tasks:
            t.max_time = max_time
            task_dict[t.task_id] = t
            state[t.task_id] = max_retry
            #self.logger.debug("sending task %s to queue", (t.task_id,t.method))

        if one_at_a_time:
            for t in tasks:     self.work_queue.put_task(t)
        else:
            self.work_queue.put_tasks(tasks)        # mxu: feed the whole set of tasks to queue at once, otherwise the queue server is busy dealing with all kinds of requests and only a very small number of tasks are reveived by workers

        while len(state):
            #self.logger.debug("%s jobs out to workers", len(state))

            results = self.work_queue.get_results(self.proj_id)

            for res in results:

                if res.task_id not in state:
                    # from a previous run that got cancelled?
                    # or a timed out job?
                    self.logger.warning("recieved result from an unknown task: %s", (res,))
                    continue

                if res.error:
                    self.logger.debug("result %s, contains error flag.  task raised exception remotely: %s", res.task_id, res)
                    state[res.task_id] -= 1
                    if state[res.task_id] > 0:
                        self.logger.warning("resubmitting crashed task: %s", res.task_id)
                        self.work_queue.put_task(task_dict[res.task_id])
                        continue
                    else:
                        self.logger.warning("task failed too many times: %s", res.task_id)

                #self.logger.debug("recieved result: %s", res)
                del state[res.task_id]
                yield res

            # mxu: when task queue (pending tasks) is empty, re-assign into unfinished tasks to idle workers
            if (self.work_queue.task_queue_size() == 0) and (len(state) > 0):

                undone_tasks = []
                for task_id in state.keys():
                    if state[task_id] <= 0: continue
                    undone_tasks.append(task_id)

                if len(undone_tasks) > 0:
                    re_submitted_tasks = self.work_queue.resubmit_undone_tasks(undone_tasks)
                    if len(re_submitted_tasks) > 0:
                        self.logger.warning("resubmitted in progress %d  tasks", len(re_submitted_tasks))

                        for task_id in re_submitted_tasks:                    state[task_id] -= 1
                        #with open('/tmp/undone_tasks-'+ self.proj_id, 'w') as f:       pickle.dump([ task_dict[_] for _ in re_submitted_tasks   ], f)          # occationally some tasks stuck at worker side. we record resubmitted tasks for diagonisis purpose


#            # Ask queue for jobs that have been submitted, but have not returned in the cutoff time?
#            tasks = self.work_queue.get_timed_out_tasks(self.proj_id)
#
#            self.logger.debug("Queue timeout.  Assuming %s jobs lost.", len(tasks))
#                # TODO: resubmit all failed jobs.
#            for task in tasks:
#                self.logger.debug("resending task %s to queue", (task))
#                self.put_task.put(task)


    # mxu: run tasks, and if there is an error, print it and raise exception
    def run__except(self, tasks, one_at_a_time=False):
        task_num = float(len(tasks))
        count = 0
        for res in self.run(tasks, one_at_a_time=one_at_a_time):
            count += 1
            sys.stdout.write('%d   %0.3f    \r'%(count, (count / task_num)));     sys.stdout.flush()

            if res.error:
                print "Computation Failed!"
                print res
                print "task_id      :", res.task_id
                print "method       :", res.method
                print "args         :", res.args
                print "kwargs       :", res.kwargs
                print "error_msg    :", res.error_msg
                raise Exception
            yield res


    # given a problem size n, see how to distribute evenly among all workers avaliable
    def estimate_chunk_size(self, n, worker_number_multiply_factor=1.5):
        n_chunk = float(n) / (self.work_queue.get_worker_number() * float(worker_number_multiply_factor))
        n_chunk = N.max(    (n_chunk, 1)    )
        n_chunk = int(N.ceil(n_chunk))

        return n_chunk
     
    #----------------------------------------------------------
    # simple broadcasting task, mainly used to clean up persistence data in worker side.
    def broadcast(self, task):


        task_dict = {}
        state     = {}
        results   = {}

        state = self.work_queue.put_broadcast_task(task)

        return state

    def broadcast_collect(self, state, timeout=120):

        cur_time = time.time()
        while len(state):
            #self.logger.debug("%s jobs out to workers", len(state))

            results = self.work_queue.get_results(self.proj_id)

            for res in results:

                if res.task_id not in state:
                    # from a previous run that got cancelled?
                    # or a timed out job?
                    self.logger.warning("recieved result from an unknown task: %s", (res,))
                    continue

                if res.error:
                    self.logger.debug("result %s, contains error flag.  task raised exception remotely: %s", res.task_id, res)

                #self.logger.debug("recieved result: %s", res)
                del state[res.task_id]
                yield res

            if (time.time() - cur_time) > timeout:      break           # we do not guarantee to get all results from all workers returned

       
    def run_broadcast__except(self, task):
        state = self.broadcast(task)

        task_num = float(len(state))
        count = 0

        for res in self.broadcast_collect(state):
            count += 1
            sys.stdout.write('%d   %0.3f    \r'%(count, (count / task_num)));     sys.stdout.flush()

            if res.error:
                print "Computation Failed!"
                print res
                print "task_id      :", res.task_id
                print "method       :", res.method
                print "args         :", res.args
                print "kwargs       :", res.kwargs
                print "error_msg    :", res.error_msg
                raise Exception
            yield res


 

if __name__ == '__main__':

    master = QueueMaster("127.0.0.1", 5011)

    tasks = [master.Task(master.proj_id, "echo", i) for i in range(10)]
    wait = [_ for _ in master.run(tasks)]

    for t in wait:
        print "Timings: task_queued: %4.2f   calculation: %4.2f,  result_queued: %4.2f" % (t.todo_queue_total, t.calc_total, t.done_queue_total)


    del master
