#!/usr/bin/env python

# delete a particular field

if __name__ == '__main__':

    import sys
    field = sys.argv[1]

    import json
    dj = json.load(sys.stdin)

    for r in dj:
        if field in r:
            del r[field]

    json.dump(dj, sys.stdout, indent=2)


