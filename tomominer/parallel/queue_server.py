#!/usr/bin/env python

import os
import sys
import time
import copy
import json
import uuid
import thread
import threading
import Queue
import logging
import logging.handlers
import psutil

from tomominer.parallel.Task import Task
from tomominer.parallel.RPCServer import RPCServer

class QueueServer:
    """
    This object is exposed by the RPCServer() as an instance.
    It is a task-queue for running work across a cluster of machines.
    It acts as middle-ware.  Any process with work to be done can connect
    with a project_id, and submit work to be completed.
    When the work is done it can be retrieved by the project.
    """

    def __init__(self):
        """
        Setup the QueueServer for task processing.
        """
        self.todo_queue = Queue.PriorityQueue()
        self.done_queues = {}
        self.done_tasks_time = dict()
        self.done_tasks_time_max = 60.0*20

        self.broadcast_todo_queue = {}

        self.start_calc  = {}
        self.out_tasks   = {}

        self.proj_alive_time = {}
        self.proj_alive_time_max = 60*5

        self.worker_alive_time = {}
        self.worker_alive_time_max = 60*30


        # configure logging.
        self.pub_logger = logging.getLogger()
        self.pub_logger.setLevel(logging.WARN)
        if False:
            h = logging.handlers.RotatingFileHandler(filename='./server.log', mode='a', maxBytes=1000000, backupCount=10)
        else:
            h = logging.handlers.RotatingFileHandler(filename='./server.log', mode='a')     # logging without rollover

        f = logging.Formatter('%(asctime)s %(host)-16s %(job_id)-16s %(source_type)-8s %(levelname)-8s %(message)s')
        h.setFormatter(f)
        self.pub_logger.addHandler(h)

        self.logger  = logging.LoggerAdapter(logging.getLogger(), {'host': os.environ.get('HOSTNAME', 'unknown'), 'job_id' : os.environ.get('PBS_JOBID', 'N/A').split('.')[0], 'source_type' : 'queue_server'})

        self.process = psutil.Process(os.getpid())

        thread.start_new_thread(QueueServer.remove_dead_projects_daemon, (self,))


    def log_queue_stats(self):
        self.logger.debug("TODO_QUEUE: %6d    IN_PROGRESS: %6d    DONE_QUEUE: %6d", self.todo_queue.qsize(), len(self.out_tasks), sum(_.qsize() for _ in self.done_queues.values()))
    
    # mxu: added for mornitoring
    def queue_stats_string(self):
        return "CPU %3.0f  PROJ %2d  WORKER %4d    TODO_QUEUE: %6d  IN_PROGRESS: %4d  DONE_QUEUE: %3d "%(self.process.cpu_percent(interval=1), len(self.done_queues), self.get_worker_number(), self.todo_queue.qsize(), len(self.out_tasks), sum(_.qsize() for _ in self.done_queues.values()))



    def new_project(self, proj_id, max_queue_size=0):
        lock = threading.Lock()
        with lock:  
            if proj_id not in self.done_queues:
                self.done_queues[proj_id] = Queue.Queue(maxsize=max_queue_size)

        self.project_keep_alive(proj_id)
        return True

    def del_project(self, proj_id=None):
        """
        Mark a project as complete. This will clean up the project data
        structures and help mark any outstanding work as done.
        If the project is not known to the server, this will return False.
        :param proj_id The project id to delete.
        """

        if proj_id in self.proj_alive_time:         del self.proj_alive_time[proj_id]

        # clean up any performing tasks with such project id
        out_tasks_to_delete = []
        for task_id in self.out_tasks:
            if self.out_tasks[task_id].proj_id in self.proj_alive_time:  continue
            out_tasks_to_delete.append(task_id)

        for task_id in out_tasks_to_delete:
            try: del self.out_tasks[task_id]
            except: pass

        # clean up all results of such project id
        #self.log_queue_stats()
        #self.logger.debug("del_project %s", proj_id)
        done_queues_to_delete = []
        for proj_id in self.done_queues:
            if proj_id in self.proj_alive_time:  continue
            done_queues_to_delete.append(proj_id)

        # remove expired tasks
        cur_time = time.time()
        for task_id in copy.deepcopy(self.done_tasks_time.keys()):
            if (task_id in self.done_tasks_time) and ( (cur_time - self.done_tasks_time[task_id]) > self.done_tasks_time_max ):
                try:
                    del self.done_tasks_time[task_id]
                except:
                    pass

        try: 
            lock = threading.Lock()
            with lock:
                for proj_id in done_queues_to_delete:                del self.done_queues[proj_id]
            return True

        except:
            return False
    

    def project_keep_alive(self, proj_id):
        self.proj_alive_time[proj_id] = time.time()
        return True

    def get_project_alive_time(self):
        return self.proj_alive_time


    # keep watching and remove dead projects using an independent thread
    def remove_dead_projects_daemon(self, interval=60):
        while True:
            time.sleep(interval)
            self.remove_dead_projects()


    def remove_dead_projects(self):
        dead_projects = []
        for proj_id_t in self.proj_alive_time:
            time_diff = time.time() - self.proj_alive_time[proj_id_t]
            if time_diff > self.proj_alive_time_max:
                dead_projects.append(proj_id_t)

        for proj_id_t in dead_projects:
            del self.proj_alive_time[proj_id_t]

        self.del_project()
        return dead_projects


    def get_worker_number(self):
        cur_time = time.time()
        worker_ids = copy.copy(self.worker_alive_time.keys())     # sometimes the self.worker_alive_time changes during the iteration, so we make a copy of keys instead of using   for _ in self.worker_alive_time
        c = 0
        for _ in worker_ids:
            if _ not in self.worker_alive_time: continue
            if (cur_time - self.worker_alive_time[_]) > self.worker_alive_time_max:     continue
            c += 1
        return c


    # mxu: added
    def task_queue_size(self):
        return self.todo_queue.qsize()
    
    # mxu: added
    def in_progress_task_number(self):
        return len(self.out_tasks)
        

    # mxu: used to feed the whole set of tasks to queue at once, otherwise the queue server is busy dealing with all kinds of requests and only a very small number of tasks are reveived by workers
    def put_tasks(self, tasks):
       print '\r' + self.queue_stats_string(),       # mxu: added for monitoring
       for task in tasks:
            self.todo_queue.put( (task.priority, task) )
            self.logger.debug("put_task %s", task)

    def put_task(self, task):
        print '\r' + self.queue_stats_string(),       # mxu: added for monitoring

        self.todo_queue.put(   (task.priority, task)   )
        self.logger.debug("put_task %s", task)
        #self.log_queue_stats()

    def cancel_task(self, task):
        raise NotImplementedError


    def get_task(self, worker_id=None, interval=5, timeout=10):
        print '\r' + self.queue_stats_string(),       # mxu: added for monitoring
        self.worker_alive_time[worker_id] = time.time()
        start_time = time.time()
        #self.log_queue_stats()
        while (time.time() - start_time) < timeout:

            try:
                priority, task = self.todo_queue.get_nowait()

                # mxu: there are cases that the task is already completed, in this case ignore it
                if task.task_id in self.done_tasks_time:
                    continue

                # record when the task was sent out for processing.
                self.start_calc[task.task_id] = time.time()

                # add it to the internal list of running tasks.
                self.out_tasks[task.task_id] = task

                # actually send the task out to the worker.
                return task

            except Queue.Empty:
                time.sleep(interval)   # sleep a bit to prevent system to be too busy, block until a certain timeout
        
        return None


#    def out_tasks_contains(self, task_id):
#        return (task_id in self.out_tasks)

    def done_tasks_contains(self, task_id):
        return (task_id in self.done_tasks_time)

    def resubmit_undone_tasks(self, task_ids):
        task_ids_t = []
        for task_id in task_ids:
            if task_id in self.out_tasks:
                self.put_task(self.out_tasks[task_id])
                task_ids_t.append(task_id)

        return task_ids_t


    def put_result(self, worker_id, task_id, error, error_msg, result):
        print '\r' + self.queue_stats_string(),       # mxu: added for monitoring

        self.worker_alive_time[worker_id] = time.time()
 
        if task_id in self.done_tasks_time:    return True     # this means the task is already processed
        self.done_tasks_time[task_id] = time.time()

        if not task_id in self.out_tasks:
            #self.logger.warning("Task is not listed as out for completion: %s", task_id)        # mxu: bug correction
            return True

        task = self.out_tasks[task_id]
        del self.out_tasks[task_id]

        task.error     = error
        task.error_msg = error_msg
        task.result    = result

        task.calc_total = time.time() - self.start_calc[task.task_id]
        del self.start_calc[task.task_id]

        if task.proj_id not in self.done_queues:    self.new_project(task.proj_id)
        self.done_queues[task.proj_id].put(task)
        self.logger.debug("put_result: %s", task)

        return True



    def get_results(self, proj_id):

        if proj_id not in self.done_queues:    self.new_project(proj_id)
        self.log_queue_stats()

        print '\r' + self.queue_stats_string(),       # mxu: added for monitoring

        results = []
        while True:
            try:
                task = self.done_queues[proj_id].get_nowait()
                results.append(task)
                
            except Queue.Empty:

                if len(results) > 0:
                    self.logger.debug("get_results: (%s)", len(results))
                else:
                    time.sleep(1)      # mxu: if no results returned, wait a while

                break
                

        # Look for AWOL tasks.
        now = time.time()

        for task_id, start_time in self.start_calc.items():
            if task_id in self.out_tasks:
                task = self.out_tasks[task_id]
                if task.proj_id not in self.done_queues:
                    del self.out_tasks[task_id]
                    continue

                if task.max_time > 0 and now - start_time > task.max_time:
                    self.logger.error("Task %s has been running for %d! Time to resubmit", task, now - start_time)


#
#        awol_tasks = [task_id for task_id,start_time in self.start_calc.items() if self.out_tasks[task_id] and now - start_time > self.out_tasks[task_id].max_time]
#        awol_tasks = [self.out_tasks[task_id] for task_id in awol_tasks]
#
#        for task in awol_tasks:
#            self.logger.error("Calculated AWOL Task: %s", task)
#
        return results


    #--------------------------------------------------------------------------------------------------
    # for broadcast task, every worker only run such task once. If the worker breaks, we dont care and we do repeat the task and send to another worker. The result of such task is still returned to the . Such broadcast task is mainly used to handle data presistance, like loading presistance data, or deleting presistance data.

    def put_broadcast_task(self, task):

        print '\r' + self.queue_stats_string(),       # mxu: added for monitoring
        
        task_worker_id = {}
        for worker_id in self.broadcast_todo_queue:
            if (time.time() - self.worker_alive_time[worker_id]) > self.worker_alive_time_max:        continue            # we do not bother a worker if it is not responding for a long time

            t = copy.deepcopy(task)                 # duplicate the task
            t.task_id_original = t.task_id
            t.task_id = str(uuid.uuid4())           # we re-assign a task id

            self.broadcast_todo_queue[worker_id].put(   (t.priority, t)   )
            self.logger.debug("put_task %s", t)
            #self.log_queue_stats()

            task_worker_id[t.task_id] = worker_id

        return task_worker_id



    def get_broadcast_task(self, worker_id):

        self.worker_alive_time[worker_id] = time.time()
 
        if worker_id not in self.broadcast_todo_queue:
            self.broadcast_todo_queue[worker_id] = Queue.PriorityQueue()
            return None

        try:
            priority, task = self.broadcast_todo_queue[worker_id].get_nowait()

            self.out_tasks[task.task_id] = task         # add it to the internal list of running tasks.
            self.start_calc[task.task_id] = time.time() # record start calculating time

            # actually send the task out to the worker.
            return task

        except Queue.Empty:
            pass

        return None



    def log(self, record):
        self.pub_logger.handle(record)

if __name__ == '__main__':

    host = "0.0.0.0"
    port = 5011

    server = RPCServer((host, port))
    server.register_instance(QueueServer())

    try:
        sys.stdout.write('starting server....');    sys.stdout.flush()
        server.serve_forever()

    except KeyboardInterrupt:
        print "Caught ctrl-c, exiting."
        server.server_close()

