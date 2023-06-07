

'''
get value of a particular field (recursively) from a record
'''
def get_field_val(rec, fields):
    rec_t = rec
    for field_t in fields:        rec_t = rec_t[field_t]
    return rec_t

