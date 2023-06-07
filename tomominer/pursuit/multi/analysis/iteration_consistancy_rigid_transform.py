#!/usr/bin/env python


# given classification iterations, inspect the rigid transform differences between consective iterations

if __name__ == '__main__':

    import os
    import json
    import classify.util as CU

    pass_i = 0
    while True:
        dj0f  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i), 'data_config.json')
        dj1f  = os.path.join(os.getcwd(), 'pass_%03d' % (pass_i + 1), 'data_config.json')

        if not os.path.exists(dj0f) :       break
        if not os.path.exists(dj1f) :       break

        with open(dj0f) as f:   dj0 = json.load(f)
        with open(dj1f) as f:   dj1 = json.load(f)

        print pass_i, CU.rigid_transform_difference(dj0, dj1)

        pass_i += 1

