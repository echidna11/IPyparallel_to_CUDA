import tomominer.io.file as IF
import copy
import tomominer.pursuit.multi.util as CU

def align_to_template(self, rec, tem_key, align_op, return_key=True):
    v =  IF.get_mrc(rec['subtomogram'])
    vm =  IF.get_mrc(rec['mask'])

    if align_op['with_missing_wedge']:
        t = IF.get_mrc(tem_key['subtomogram'])
        tm = IF.get_mrc(tem_key['mask'])
        re = CU.align_vols_with_wedge(v1=t, m1=tm, v2=v, m2=vm, op=align_op)
    else:
        t = IV.get_mrc(t_key['subtomogram'])
        re = CU.align_vols_no_wedge(v1=t, v2=v, op=align_op)

    # re have angle, loc, score, err
    re["subtomogram"] = rec["subtomogram"]
    re["mask"] = rec["mask"]
    re["id"] = rec["id"]
    if not isinstance(re["loc"], list):
        re["loc"] = re["loc"].tolist()
    if not isinstance(re["angle"], list):
        re["angle"] = re["angle"].tolist()

    if return_key:
        re_key = self.cache.save_tmp_data(re, fn_id=self.task.task_id)
        assert re_key is not None
        return {'key':re_key}
    else:
        return re