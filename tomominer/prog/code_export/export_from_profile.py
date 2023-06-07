#!/usr/bin/env python



'''

first, run profiler to collect all functions and libraries that have been accessed, there can be multiple profiler results
then use inspect module to collect the corresponding source code, and export them to another directory, 


to collect all profile files, use
find /tmp -name "profile-*.prof" > profile-file-list.txt


TODO: add check of validaty of a profile by checking the AST node type of the corresponding line



~/ln/tomominer/tomominer/prog/code_export/export_from_profile.py

'''

import os, json, pstats, inspect, copy, Queue, astor, _ast


def collect_functions(profile_files, dirs, verbose=False):
    fs = []
    for pf in profile_files:
        ps = pstats.Stats(pf)
        fst = [{'src':_[0], 'line_no':_[1], 'func':_[2]} for _ in ps.stats]        # collect all functions and their residing files
        fs.extend(fst)

    del ps

    # keep only python scripts
    fs = [_ for _ in fs if (_['src'].endswith('.py'))]

    # re-arrange according to src
    src_f = {}          # grouping of functions by souce files, the functions are indexed by their line numbers
    for f in fs:
        s = f['src']
        if s not in src_f:     src_f[s] = {'func':{}}

        f_line_no = f['line_no']
        f_func = f['func']

        if (f_line_no, f_func) in src_f[s]['func']:     continue

        src_f[s]['func'][(f_line_no, f_func)] = ({'line_no':f_line_no, 'func':f['func']})


    # print all collected functions
    if verbose:
        for s in src_f:
            print s
            for fi in src_f[s]['func']:
                f = src_f[s]['func'][fi]
                print '\t', f['func'], '\t', f['line_no']

    # form dest path
    dest_src = {}
    for s in src_f:
        hit = None

        for r in dirs:
            if s.startswith(r['src']):
                hit = r
                break
        if hit is None:     continue

        dest_t = os.path.join(hit['dest'], s[len(hit['src']):])

        # we require one to one mapping between src and dest file
        assert dest_t not in dest_src

        src_f[s]['dest'] = dest_t
        dest_src[dest_t] = s
        

    assert  len(dest_src) > 0
    src_f = {_:src_f[_] for _ in src_f if ('dest' in src_f[_])}

    # print all collected functions, with destinations
    if verbose:
        for s in src_f:
            print s, src_f[s]['dest']
            for fi in src_f[s]['func']:
                f = src_f[s]['func'][fi]
                print '\t', f['func'], '\t', f['line_no']

    return {'src_f':src_f, 'dest_src':dest_src}




def collect_ast_trees(fs):


    src_f = fs['src_f']

    for s in src_f:
        with open(s) as f:      src = f.readlines()
        src_f[s]['src'] = src

        ast_root = astor.parsefile(s)

        src_f[s]['ast'] = ast_root





'''
regarding to the module imports, for every xxx.py file executed or imported, all imported modules inside xxx.py should already appear in fs_module because they are executed, and they are supposed to be copied to the destination
so we don't need to clean up imports but just put them into the destnation xxx.py
'''

def export_functions(fs, header='', verbose=False):

    src_f = fs['src_f']

    for s in src_f:
        dest = src_f[s]['dest']
        assert      not os.path.isfile(dest)
        if not os.path.isdir(os.path.dirname(dest)):        os.makedirs(os.path.dirname(dest))

        if verbose:     print s, '\t', dest
        ast_t = extract_ast(ast_root=src_f[s]['ast'], func=src_f[s]['func'])
        
        if ast_t is None:
            open(dest, 'a').close()     # just create an empty file
            continue

        export_code(src=astor.to_source(ast_t), file_out=dest, header=header)


    # copy rest accessed modules and files to destination dir, so far, it only removes comments starts with '#'. In future, need to remove comments from quotes
    for s in src_f:
        dest = src_f[s]['dest']

        if os.path.isfile(dest):       continue
        if not os.path.isdir(os.path.dirname(dest)):        os.makedirs(os.path.dirname(dest))

        print s, '\t', dest
        export_code(src=astor.to_source(src_f[s]['ast']), file_out=dest, header=header)




def ast_walk(r):
    q = Queue.Queue()
    q.put(r)

    while not q.empty():
        n = q.get()

        yield n

        if 'body' in set(dir(n)):
            for nt in n.body:   q.put(nt)


def ast_add_parent_field(r):
    r.node_parent = None

    q = Queue.Queue()
    q.put(r)

    while not q.empty():
        n = q.get()

        if 'body' in set(dir(n)):
            for nt in n.body:
                nt.node_parent = n
                q.put(nt)




def extract_ast(ast_root, func):

    line_nos__func_map = {func[_]['line_no']:func[_]['func'] for _ in func}

    line_nos = [func[_]['line_no'] for _ in func]
    line_nos = set(line_nos)

    if len(line_nos) == 0:      return None


    ast_root = copy.deepcopy(ast_root)
    ast_add_parent_field(ast_root)

    for n in ast_walk(ast_root):        n.node_keep = False         # this field is used to indicate whether the node should be kept


    # modify ast tree to keep functions
    q = Queue.Queue()
    q.put(ast_root)

    while not q.empty():
        n = q.get()
        n_fields = set(dir(n))

        if 'body' in n_fields:
            for nt in n.body:   q.put(nt)


        if 'lineno' not in n_fields:    continue
        if n.lineno not in line_nos:    continue
        if (type(n) is not _ast.FunctionDef) and (type(n) is not _ast.ClassDef):         # we only extract functions and classes
            if not line_nos__func_map[n.lineno].startswith('<'):   print 'WARNING: omitting', n.lineno, line_nos__func_map[n.lineno], type(n)
            continue           # we only consider functions, in such case, the static variables in class def may be omitted, be aware!!

        nt = n
        while nt is not None:
            nt.node_keep = True
            nt = nt.node_parent

        for nt in ast_walk(n):      nt.node_keep = True


    if not ast_root.node_keep:  return None         # this means no functions are actually called

    # include import nodes
    for n in ast_root.body:
        if type(n) is _ast.Import:        n.node_keep = True
        if type(n) is _ast.ImportFrom:        n.node_keep = True

   

    # trim tree
    q = Queue.Queue()
    q.put(ast_root)

    while not q.empty():
        n = q.get()
        n_fields = set(dir(n))

        if 'body' in n_fields:
            body_t = []
            for nt in n.body:
                if not nt.node_keep:    continue
                q.put(nt)
                body_t.append(nt)
            n.body = body_t

    return ast_root
   

def export_code(src, file_out, header='', executable=False):

    if os.path.isfile(file_out):
        raise Exception('Destination file exists, need to remove first ' + file_out)

    with open(file_out, 'w') as f:
        if executable:
            print >>f, '#!/usr/bin/env python'          # write this header if we need automatic starting, because it is removed in parsing

        print >>f, ''
        print >>f, header
        print >>f, ''

        f.write(quote_comment_removal(src))


'''
after getting source code using astor.to_source(), the comment of triple-quoted strings is compressed into a single line, which can be easily identified and omitted
'''
def quote_comment_removal(txt):
    import StringIO
    buf = StringIO.StringIO(txt)

    ss = buf.readlines()
    
    txt_new = ''
    for s in ss:
        if s.lstrip().startswith('\''):     continue
        if s.lstrip().startswith('\"'):     continue
        txt_new += s

    return txt_new
 

def main():
    with open('export_from_profile__op.json') as f:     op = json.load(f)

    print 'WARNING: function calls from multiprocessing will not be recorded by cProfile, so they cannot be automatically exported'

    # make sure that the destination dirs do not exist, this is used to prevent overwriting of exising files and dirs
    for r in op['dirs']:
        assert  not os.path.isdir(r['dest'])

    for d in op['dirs']:        d['dest'] = os.path.abspath(d['dest'])

    with open(op['profile files list']) as f:       pf = f.readlines()
    pf_new = []
    for pf_t in pf:
        pf_t = pf_t.strip()
        if not os.path.isfile(pf_t):        continue
        pf_new.append(pf_t)

    print 'collecting functions'
    fs = collect_functions(profile_files=pf_new, dirs=op['dirs'], verbose=op['verbose'])

    print 'loading and parsing source code'
    collect_ast_trees(fs)

    with open(op['header file']) as f:        header_t = f.read()

    print 'exporting code'
    export_functions(fs=fs, header=header_t)


if __name__ == '__main__':
    main()


