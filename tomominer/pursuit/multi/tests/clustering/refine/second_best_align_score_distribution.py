
'''

given one pattern, find the set of subtomograms that best aligns to it.
then draw histogram / the emperical distribution of best and second best alignment scores
then use simple linear SVM to find the best discrimination score


see idea 6 of
https://docs.google.com/document/d/1NmYhEczJ8A--KI5PStf8hyuPM5QPfRLz7RftOfkn8_k/edit



'''


# commands



import os, sys
import cPickle as pickle

import numpy as N


data_dir = '/home/rcf-47/mxu/ln/frequent_structure/data/out/method/imputation-classification/multi/classification/30/0.005/0130/classify/pass_029'
with open(os.path.join(data_dir, 'align_template.pickle'), 'rb') as f:     al = pickle.load(f)
al = [_.result['align'] for _ in al]



import tomominer.pursuit.multi.recursive.filter.second_largest_cut as PMRFS

template_id=0
alb = PMRFS.collect_items(template_id, al)
print len(alb)



x = []
y = []
for albt in alb:
    s0, s1 = PMRFS.get_best_two_scores(albt)
    x.append([s0])            ;       y.append(1)
    x.append([s1])            ;       y.append(-1)

x = N.array(x)
y = N.array(y)

from sklearn import svm
c = svm.LinearSVC()
c.fit(x, y)

c.predict(x)
c.predict(x) - y

c.intercept_
c.coef_

cutoff = (- c.intercept_) / c.coef_
cutoff = float(cutoff)

N.sign(x.flatten() - cutoff) - c.predict(x)


