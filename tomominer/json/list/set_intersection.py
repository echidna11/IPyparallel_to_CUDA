

# given a field name, collect all items from one JSON file
# then use this as a filter to filter another JSON file (fed in through stdin)


if __name__ == '__main__':
    
    import sys
    field = sys.argv[1]
    set_file = sys.argv[2]

    import json
    with open(set_file) as f:   set_json = json.load(f)

    set_items = set( [rec[field] for rec in set_json] )


    filter_json = json.load(sys.stdin)


    filter_json_f = [rec for rec in filter_json if (rec[field] in set_items) ]
    
    json.dump(filter_json_f, sys.stdout, indent=2)


