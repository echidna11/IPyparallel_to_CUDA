#!/usr/bin/env python



'''

export scores for plotting box plots, need to order the columns according to decreasing of median

~/ln/tomominer/tomominer/template/library/score_export.py

'''

import os, json, csv
import numpy as N


def main():

    with open('score_export__op.json') as f:        op = json.load(f)

    if not os.path.isdir(op['out dir']):        os.makedirs(op['out dir'])
    
    alop_dir = os.path.dirname(os.path.abspath(op['reference multi op']))
    with open(op['reference multi op']) as f:      alop = json.load(f)

    with open(os.path.join(alop_dir, alop['out_file'])) as f:      al = json.load(f)

    with open(os.path.join(alop_dir, alop['reference_json_file'])) as f:    ref = json.load(f)

    if 'selected data' in op:
        with open(op['selected data']) as f:      dj_sel = json.load(f)
        with open(os.path.join(alop_dir, alop['data_json_file'])) as f:      dj = json.load(f)
        dj_i = {dj[_]['subtomogram']:_ for _ in range(len(dj))}
        ind_sel = [dj_i[_['subtomogram']] for _ in dj_sel]

    else:
        ind_sel = None


    s = score_matrix(al)

    pids = [_['pid'] for _ in ref]

    export_tab(s, ids=pids, out_file=os.path.join(op['out dir'], op['score all out']))
    if ind_sel is not None:
        export_tab(s[ind_sel,:], ids=pids, out_file=os.path.join(op['out dir'], op['score seleted out']))
        export_tab(N.delete(s, ind_sel, axis=0), ids=pids, out_file=os.path.join(op['out dir'], op['score rest out']))

        # export rank of selected subtomograms
        export_tab(rank(s)[ind_sel,:] / float(s.shape[0]), ids=pids, out_file=os.path.join(op['out dir'], op['score seleted rank out']))

        # export enrichment
        export_enrichment(s=s, ind_sel=ind_sel, ids=pids, out_dir=os.path.join(op['out dir'], 'enrichment'))


def score_matrix(al):
    # collect a matrix of scores
    s = []
    for i, alt in enumerate(al):
        ss = []
        for ri in range(len(alt)):
            ss.append(alt[str(ri)]['score'])

        s.append(ss)

    s = N.array(s)

    return s


def export_tab(s, ids, out_file):
    print 'export_tab()', s.shape, out_file

    # order scores and references in descending of median
    sm = N.median(s, axis=0)
    sm_i = N.argsort(-sm)

    s_s = s[:,sm_i]
    ids_s = [ids[_] for _ in sm_i]


    with open(out_file, 'w') as f:
        cw = csv.writer(f)
        cw.writerow([_ for _ in ids_s])
        cw.writerows(s_s.tolist())



def rank(s):
    ss = N.argsort(s, axis=0)

    r = N.empty(s.shape, dtype=N.int)
    for j in range(s.shape[1]):
        r[ss[:,j] ,j] = N.arange(s.shape[0])

    return r


# calculate and export enrichments, in terms of the number of iterms in ind_sel divided with the total number of iterms, given certain score cutoff
def export_enrichment(s, ind_sel, ids, out_dir):
    if not os.path.isdir(out_dir):      os.makedirs(out_dir)

    ind_sel = set(ind_sel)

    for col in range(s.shape[1]):
        
        st = s[:, col].flatten()
        st_i = N.argsort(-st)

        score = []
        proportion_against_all = []         # among total selected items above cutoff, how many are inside ind_sel
        proportion_against_sel = []         # among total iterms in ind_sel, how many are above cutoff

        sel_num = 0
        rest_num = 0

        for i in st_i:
            if i in ind_sel:
                sel_num += 1
            else:
                rest_num += 1

            score.append(st[i])
            proportion_against_all.append(float(sel_num) / (sel_num + rest_num))
            proportion_against_sel.append(float(sel_num) / len(ind_sel))
   
        out_file = os.path.join(out_dir, ids[col] + '.csv')
        with open(out_file, 'w') as f:
            cw = csv.writer(f)
            cw.writerow(['score', 'all', 'sel'])
            cw.writerows(N.array([score, proportion_against_all, proportion_against_sel]).T.tolist())


if __name__ == '__main__':
    main()



'''
# R code for plotting boxplot


setwd('out')


#-------------------------------------------------------------------------------------------
# plot distribution of selected scores
ss = read.csv('score_export__score_selected_out.csv', row.names=NULL, check.names=FALSE)
setEPS()
postscript("score_export__score_selected_out__boxplot.eps")
boxplot(ss, xlab='PDB ID', ylab='Score', las=2, main='Alignment against templates of multiple complexes')
dev.off()


# calculate wilcoxon test p-values pairwisely between templates
mat = c()
for(pid0 in colnames(ss)) {
    ws = c()
    for(pid1 in colnames(ss)) {
        w = wilcox.test(x=ss[,pid0], y=ss[,pid1], alternative='greater')
        ws = c(ws, w$p.value)
    }
    mat = rbind(mat, ws)
}

max(mat[1, 2:dim(mat)[1]])      # the highest p-value between the top hit and the rest hits
max(mat[1:2, 3:dim(mat)[1]])    # the highest p-value between the first two hits that rest hits

#-------------------------------------------------------------------------------------------
# plot distribution of rank of selected scores
r = read.csv('score_export__score_selected_rank_out.csv', row.names=NULL, check.names=FALSE)
setEPS()
postscript("score_export__score_selected_rank_out__boxplot.eps")
boxplot(r, ylab='Score rank ratio', las=2, main='Alignment against templates of multiple complexes')
dev.off()



#-------------------------------------------------------------------------------------------
# plot for a single complex

s = read.csv('score_export__score_all_out.csv', row.names=NULL, check.names=FALSE)
sr = read.csv('score_export__score_rest_out.csv', row.names=NULL, check.names=FALSE)


for(pid in colnames(s)) {

#pid = '1KP8'
#pid = '2J00-2J01'


setEPS()
postscript(paste('enrichment/', pid, '--qq-against-all', '.eps', sep=''))
qqplot(s[,pid], ss[,pid], xlab='All', ylab='Selected')           # plot QQ plot of distribution of selected scores, against distribution of all scores
dev.off()

setEPS()
postscript(paste('enrichment/', pid, '--qq-against-rest', '.eps', sep=''))
qqplot(sr[,pid], ss[,pid], xlab='Rest', ylab='Selected')          # plot QQ plot of distribution of selected scores, against distribution of rest scores
dev.off()



setEPS()
postscript(paste('enrichment/', pid, '--against-rest', '.eps', sep=''))
boxplot(list(Selected=ss[,pid], Rest=sr[,pid]), ylab='Score', las=2, main=paste('Alignment against template of', pid))
dev.off()

# calculate the corresponding p-value of above box plot, which is the alignment scores of the pattern's subtomograms vs all rest subtomograms
w = wilcox.test(x=ss[,pid], y=sr[,pid], alternative='greater')
w$p.value



# plot enrichment

rin = read.csv(paste('enrichment/', pid, '.csv', sep=''), row.names=NULL, check.names=FALSE)

setEPS()
postscript(paste('enrichment/', pid, '--enrich', '.eps', sep=''))
plot(0,0, xlim=c(min(rin[,'score']), max(rin[,'score'])), ylim=c(0,1), type='n', xlab='Score cutoff', ylab='Proportion of selected subtomograms above score cutoff', main=paste('Enrichment of template search against', pid))
lines(rin[,'score'], rin[,'all'], type='l', col='blue')
#lines(rin[,'score'], rin[,'sel'], type='l', col='green')
dev.off()


}





# plot two densities, see  http://stackoverflow.com/questions/3541713/how-to-plot-two-histograms-together-in-r
library('ggplot2')
score_sel <- data.frame('Alignment score' = ss[,pid])
score_rest <- data.frame('Alignment score' = sr[,pid])

score_sel$set = 'Selected'
score_rest$set = 'Rest'

score_combined = rbind(score_sel, score_rest)
ggplot(score_combined, aes('Alignment score', fill = set)) + geom_density(alpha = 0.2)



'''

