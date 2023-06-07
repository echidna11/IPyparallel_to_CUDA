#!/usr/bin/python
# -*- coding: utf-8 -*-

import SocketServer
import pickle
import time
import sys
import select
import os
import psutil


class RPCHandler(SocketServer.StreamRequestHandler):
    """
    Connection object that communicates with a RPCClient object.
    """

    def handle(self):
        """
        Connect to the socket and send pickled request/response across until
        the connection is closed.

        Within the connection, we enter a loop loading pickled data, processing
        it and returning results.
        """

        self.server.active_connections+=1
        while True:

            try:
                data = pickle.load(self.rfile)
            except EOFError:
                # EOF means we're done with this request.
                break

            try:
                result = self.server._dispatch(data)
            except Exception, e:
                pickle.dump(('ERR', e), self.wfile, protocol=2)

                if True:
                    import traceback
                    traceback.print_exc()
            else:
                pickle.dump(('OK', result), self.wfile, protocol=2)

            self.wfile.flush()

        self.server.active_connections-=1

class RPCServer(SocketServer.ThreadingTCPServer):
    """
    The actual server itself.  This is the TCP server which creates a new
    thread for every connection.

    MixIn style class for the server.  This adds the capability for registering
    an object that supplies the methods that are looked up.  It also supplies
    the dispatch function for calling the mapped method.
    """
    daemon_threads      = True
    allow_reuse_address = True
    active_connections  = 0


    def __init__(self, addr, requestHandler=RPCHandler, bind_and_activate=True):
        self.instance = None
        SocketServer.ThreadingTCPServer.__init__(self, addr, requestHandler, bind_and_activate)


        self.previous_method = None
        self.same_method_call_count = 0

        self.process = psutil.Process(os.getpid())

    def register_instance(self, obj):
        """
        Register an object which will provides the functions to be called remotely.

        :param obj  the object that supplies all of the methods that will be
                    called remotely.
        """
        self.instance = obj

    def _dispatch(self, data, cpu_usage_threshold=200):
        """
        Run a requested function on provided arguments.  Each request contains
        a function name, arguments, and keyword arguments.  The function name
        must match a method of the registered instance class.  It must also not
        start with underscores, as an attempt to protect private methods.

        :param data The request to process.
        """

        try:
            method = data['method']
            args = data['args']
            kwargs = data['kwargs']
        except:
            raise

        if self.instance is None:
            raise Exception("No instance installed on the server.")


        if False:
            # display the method name and number of repeated calls
            if method != self.previous_method:
                self.previous_method = method
                self.same_method_call_count = 0
                sys.stdout.write('\n')

            self.same_method_call_count += 1
            sys.stdout.write('\r'+method + ' ' + repr(self.same_method_call_count) + '\t')

        while True:
            # monitor CPU usage of the server process, if the process is too busy, wait for a while
            cpu_usage = self.process.cpu_percent(interval=1)
            if cpu_usage < cpu_usage_threshold:     break

            #sys.stdout.write(' ' + repr(cpu_usage) + ' ')
            time.sleep(1)


        if method.startswith("_"):
            raise AttributeError("Cannot call method (%s) with leading '_'" % (method))
        if hasattr(self.instance, method):
            func = getattr(self.instance, method)
            if not callable(func):
                raise AttributeError("Requested function (%s) is not callable" % (method))
            return func(*args, **kwargs)
        else:
            raise AttributeError("Requested function (%s) not found in instance", (method))



    '''
    def serve_forever(self, poll_interval=0.5):
        """Handle one request at a time until shutdown.

        Polls for shutdown every poll_interval seconds. Ignores
        self.timeout. If you need to do periodic tasks, do them in
        another thread.
        """
        self._BaseServer__is_shut_down.clear()
        try:
            while not self._BaseServer__shutdown_request:
                # XXX: Consider using another file descriptor or
                # connecting to the socket to wake this up instead of
                # polling. Polling reduces our responsiveness to a
                # shutdown request and wastes cpu at all other times.
                r, w, e = select.select([self], [], [], poll_interval)
                if self in r:
                    self._handle_request_noblock()
                    #sys.stdout.write(':')
                else:
                    #sys.stdout.write('.')
                    pass

        finally:
            self.__shutdown_request = False
            self._BaseServer__is_shut_down.set()
    '''

    # this overrides BaseServer.handle_error to reduce massy output
    def handle_error(self, request, client_address):
        return

        print '-'*40
        print 'Exception happened during processing of request from',
        print client_address
        import traceback
        traceback.print_exc() # XXX But this goes to stderr!
        print '-'*40


if __name__ == '__main__':

    host = "localhost"
    port = 8000

    class Test(object):
        def echo(self, data):
            '''Method for returning data got from client'''
            return data

        def div(self, dividend, divisor):
            '''Method to divide 2 numbers'''
            return dividend/divisor

        def is_computer_on(self):
            return True

        def raising_function(self, arg):
            if not arg:
                raise Exception

        def _private_fn(self, x):
            return x + 3

    server = RPCServer((host, port), RPCHandler)
    server.register_instance(Test())

    try:
        if False:
            server.serve_forever()
        else :
            t = threading.Thread(target=server.serve_forever)
            t.setDaemon(True) # don't hang, otherwise the threads do not respond to ctrl-c
            t.start()

            while True: time.sleep(100)

    except KeyboardInterrupt:
        print 'Exiting...'
        server.server_close()
