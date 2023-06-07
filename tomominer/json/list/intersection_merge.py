

# given a field name that defines keys, and multiple json files both containint list of records, find the interection set of records, then combine all fields
# WARNING: the overlap fields will be overwritten/updated without warning, so the order of file_names is important!!

if __name__ == '__main__':

    import sys
    key_field = sys.argv[1]
    file_names = sys.argv[2:]

    import json
    rec_json_s = [None]*len(file_names)
    for fi, fn in enumerate(file_names):
        with open(fn) as f:     rec_json_s[fi] = json.load(f)

    rec_keys = [None]*len(rec_json_s)
    for i, rec_json in enumerate(rec_json_s):
        rec_keys[i] = [rec[key_field] for rec in rec_json]


    common_key_set = set(rec_keys[0])
    for i in range(1, len(rec_keys)):
        common_key_set.intersection_update(set(rec_keys[i]))


    from collections import defaultdict
    rec_merge = defaultdict(dict)
    for rec_json in rec_json_s:
        for rec in rec_json:
            key_t = rec[key_field]
            if key_t not in common_key_set:     continue

            for field_t in rec:
                rec_merge[key_t][field_t] = rec[field_t]


    rec_merge_json = [rec_merge[k] for k in rec_merge]
    json.dump(rec_merge_json, sys.stdout, indent=2)


