#!/usr/bin/env python

'''
simple script to call the queue server to remove dead projects
'''

if __name__ == '__main__':

    from tomominer.parallel import RPCClient
    q = RPCClient('127.0.0.1', 5011)

    import time
    at = q.get_project_alive_time()
    at = [(time.time() - at[_], _) for _ in at]
    at = sorted(at, key=lambda _:_[0], reverse=True)

    print 'project alive time'
    for _ in at:    print '%0.1f'%(_[0],), '\t\t', _[1]

    import sys
    if len(sys.argv) > 1:
        for i in range(1, len(sys.argv)):
            proj_id = sys.argv[i]
            print 'trying to remove project', proj_id
            q.del_project(proj_id)
    else:
        print 'trying to remove following projects', q.remove_dead_projects()

