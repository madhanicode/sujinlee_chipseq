import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from dani_analysis import which_hist, get_strain
from matplotlib.patches import Patch

#usage notes:
#activate dani environment (conda activate ~/miniconda3/envs/dani_ipynb_dec_29/)
#check for lines with comment "change me to". fix directory
#ensure that text "cen" or "tel" does not appear in path to files!
#see dani_analysis for coverage over pseudo-centromere (Figure 4 coverage tracks)

def execute():
    kind ='SLHM25_26'

    #change me!
    RPKM_input_dir = '/home/manning/scripts/sujin_testing/RPKM_calculations_input.txt'
    rpkm_input = rpkm_calcs(RPKM_input_dir)

    #change me too!
    BASE_DIR = '/home/manning/scripts/sujin_testing/'
    PLOT_DIR = BASE_DIR + 'temp_bargraphs_2/'
    isExist = os.path.exists(PLOT_DIR)
    if not isExist:
        os.makedirs(PLOT_DIR)

    r = rpkm_input[rpkm_input['hist']!='WCE']
    header = ["strain","hist","rpkm_cen","rpkm_tel","rpkm_not","tel_norm","cen_norm"]
    r.to_csv('out.csv', index=False, columns=header)
    #return

    tit5 = 'RPKM for Telomeres vs. Centromeres\nInput Normalized + Background Subtracted\n '

    settings = {'SLHM23':   {'s1':{'ylim':[-200,300], 'nrow':1, 'ncol':4, 'legend_i':4},
                            's2':{'ylim':[-4,4],     'nrow':1, 'ncol':4, 'legend_i':4}},
                'SLHM20_21':{'s1':{'ylim':[-800,500], 'nrow':2, 'ncol':4, 'legend_i':8},
                            's2':{'ylim':[-4,4],     'nrow':2, 'ncol':4, 'legend_i':8}},
                'SLHM17_18':{'s1':{'ylim':[-800,500], 'nrow':2, 'ncol':4, 'legend_i':8},
                            's2':{'ylim':[-4,4],     'nrow':2, 'ncol':4, 'legend_i':8}},
                
                'SLHM19':   {'s1':{'ylim':[-200,300], 'nrow':1, 'ncol':4, 'legend_i':4},
                            's2':{'ylim':[-2,3],     'nrow':1, 'ncol':4, 'legend_i':4}},
                
                'SLHM25_26':{'s1':{'ylim':[-200,420], 'nrow':1, 'ncol':4, 'legend_i':4},
                            's2':{'ylim':[-4,4],     'nrow':1, 'ncol':4, 'legend_i':4}}}

    fig5 = plt.figure(figsize=[20,5])
    fig5 = up_down_plot_all(rpkm_input, fig5, PLOT_DIR, title=tit5, y1lab='tel_norm', y2lab='cen_norm', 
                            ylim=settings[kind]['s1']['ylim'], nrow=settings[kind]['s1']['nrow'], 
                            ncol=settings[kind]['s1']['ncol'], legend_i=settings[kind]['s1']['legend_i'])


def rpkm_calcs(rpkm_file):
    rpkm = pd.read_csv(rpkm_file, sep='\s')
    rpkm['hist'] = rpkm['sample'].apply(lambda x: which_hist(x))
    rpkm['strain'] = rpkm['sample'].apply(lambda x: get_strain(x))
    rpkm['rpkm_cen'] = rpkm.R_cen*1000000000 / (rpkm.L_cen * rpkm.S )
    rpkm['rpkm_tel'] = rpkm.R_tel*1000000000 / (rpkm.L_tel * rpkm.S )
    rpkm['rpkm_not'] = rpkm.R_not*1000000000 / (rpkm.L_not * rpkm.S )
    rpkm['log_cen'] = np.log( rpkm.rpkm_cen / rpkm.rpkm_not )
    rpkm['log_tel'] = np.log( rpkm.rpkm_tel / rpkm.rpkm_not )
    rpkm['tel_norm'] = rpkm['rpkm_tel'] - rpkm['rpkm_not']
    rpkm['cen_norm'] = rpkm['rpkm_cen'] - rpkm['rpkm_not']
    rpkm['tel_norm_fc'] = np.log( rpkm['rpkm_tel'] / rpkm['rpkm_not'] )
    rpkm['cen_norm_fc'] = np.log( rpkm['rpkm_cen'] / rpkm['rpkm_not'] )
    return rpkm

# "Up-down" plots

def autolabel(rects, texts,above_below='above'):
    heights = [rect.get_height() for rect in rects]
    height  = max(heights) if above_below=='above' else min(heights)
    for i,rect in enumerate(rects):
        x = rect.get_x() + rect.get_width()/2.
        y = 1.05*height
        plt.text(x, y, texts[i], ha='center', va='bottom')
        
def get_labels(h, r, y1lab, y2lab):
    mapping = {1.0:'',-1.0:'(-)'}
    labels1 = [mapping[np.sign(y)] for y in r[y1lab]]
    labels2 = [mapping[np.sign(y)] for y in r[y2lab]]
    return [labels1, labels2]


def up_down_plot(h, fig, rpkm, i=1, ylabel='RPKM', xlabel='', ylim=[-600,500], #xlabel='strain',
                y1lab='rpkm_tel', y2lab='rpkm_cen', nrow=2, ncol=4, tel_color='b', cen_color='y'):
    strain_map = {'CK2308': 'CCC1\u0394',
    'CK6644': 'ezh2\u0394',
    'SL362': '4KR',
    'SL363': '6KR',
    'SL364': '10KR',
    'SL365': '4DE',
    'SL366': '5DE',
    'SL367': '9DE',
    'SL368': '6DE',
    'SL369': 'CC\u0394',
    'SL370': 'Y31A',
    'SL371': 'W52A',
    'SL372': 'WT (untagged)',
    'SL3': 'WT'}    
    r_temp = rpkm[rpkm['hist']==h]

    #change me to exclude custom samples from graph
    drop_list = ['SL368', 'SL370', 'SL371', 'SL3']
    r = r_temp.drop(r_temp.loc[r_temp['strain'].isin(drop_list)].index, inplace=False)
    
    #assign new variable to each row for plot order and sort
    #change me for custom reindexing
    plot_order = {'CK2308': 2,
    'CK6644': 1,
    'SL362': 4,
    'SL363': 5,
    'SL364': 6,
    'SL365': 7,
    'SL366': 8,
    'SL367': 9,
    'SL369': 3,
    'SL372': 0}    
    r['plot_order'] = r['strain'].apply(lambda x: plot_order[x])
    r = r.sort_values(by = 'plot_order', ascending = True)

    y1 = [abs(x) for x in r[y1lab]]
    y2 = [-abs(x) for x in -r[y2lab]]
    x = [strain_map[x] for x in r['strain']]


    ax = fig.add_subplot(nrow,ncol,i)
    r1 = ax.bar(x, y1, width=.8, color=tel_color)
    r2 = ax.bar(x, y2, width=.8, color=cen_color)
    #adds a (-) sign above "negative" rpkm value columns
    text1, text2 = get_labels(h, r, y1lab, y2lab)
    autolabel(r1, text1, 'above')
    autolabel(r2, text2, 'below')
    
    plt.xticks(rotation = 90)
    plt.ylim(ylim)
    plt.xlabel(xlabel)

    plt.ylabel(ylabel)
    plt.title(h)
    
    return fig

def up_down_plot_all(rpkm, fig, PLOT_DIR, tel_color='b', cen_color='y',
                     title='title', y1lab='rpkm_tel', y2lab='rpkm_cen', ylim=[-100,100],
                     nrow=2, ncol=4, legend_i=8):

    for i,h in enumerate(list(set(rpkm['hist']))):
        fig = up_down_plot(h, fig, rpkm, i+1, y1lab=y1lab, y2lab=y2lab, ylim=ylim, nrow=nrow, ncol=ncol)
    fig.suptitle(title)

    legend_elements = [Patch(facecolor=tel_color, edgecolor=tel_color, label='Telomeres'),
                       Patch(facecolor=cen_color, edgecolor=cen_color, label='Centromeres'),
                       Patch(facecolor='white',   edgecolor='white',   label='(-) indicates value is below background')]

    ax = fig.add_subplot(nrow,ncol,legend_i)
    plt.legend(handles=legend_elements, loc='upper left')
    plt.axis('off')
    t = PLOT_DIR + 'RPKM_' + title.split('\n')[1].replace(' ','_')
    plt.savefig(t + '.svg')
    plt.savefig(t + '.png')
    return fig

if __name__ == '__main__':
    execute()