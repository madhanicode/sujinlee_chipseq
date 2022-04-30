import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FormatStrFormatter
from scipy.interpolate import make_interp_spline, BSpline
import re, glob
import os
from tqdm import tqdm

#pipeline works as follows, reconstructed:
#1. sujin runs commands in jupyter notebook for jordan's pipeline (separate pipeline for chip-seq)
#   this generates .bedgraph files for alignment of reads to genome
#2. get_background.sh is run, using bedtools to pull reads mapping to non telomeric or centromeric
#   regions and get average read depth
#3. RPKM.sh? may have been run, this file contains lines that convert the WCE .bedgraph to _cen.bedgraph or _tel.bedgraph
#   calls bedtools map -o sum -c 4 -a ${refdir}/tel_windows_sorted.bed -b ${f} > ${o}/${b}_tel.bedgraph
#   or equivalent for cen. This depends on definitions in tel/cen_windows_sorted .bed files.
#   wrote pared down script, WCE_bedgr_to_cen_tel.sh to call command for all outputs on sujin bedgraph files

#usage notes:
#activate dani environment (conda activate ~/miniconda3/envs/dani_ipynb_dec_29/)
#check for lines with comment "change me to". fix directory
#ensure that text "cen" or "tel" does not appear in path to files!
#see dani_analysis_pt2 for up down bargraphs (Figure 4 blue + gold bargraphs)

def execute():
    kind ='SLHM25_26'
    BASE_DIR = '/home/manning/scripts/sujin_testing/'
    outdir   = kind +'/'

    PLOT_DIR = BASE_DIR + 'temp_graphs_subtract_line/'
    isExist = os.path.exists(PLOT_DIR)
    if not isExist:
        os.makedirs(PLOT_DIR)

    bedgraphdir = BASE_DIR + 'bedgraph_ct/'
    cen_files = glob.glob(bedgraphdir+'*_cen.bedgraph')
    tel_files = glob.glob(bedgraphdir+'*_tel.bedgraph')

    #no_spike_cen_files = glob.glob(DIR+'*cen_input.bedgraph')
    #no_spike_tel_files = glob.glob(DIR+'*tel_input.bedgraph')
    table = make_table(tel_files)

    ###README
    ###this chunk generates a graph for each chromosome, comparing wt vs a given bedgraph/strain

    # SLHM25_26 variables (no spike in control)
    #DATA_NO_SPIKE = '/data/sujin/SLHM25_26/MULTI/bedgraphs_WCE'
    WT = 'SL3'
    #no_spike_cen_files = glob.glob(DIR+'*_cen.bedgraph')
    #no_spike_tel_files = glob.glob(DIR+'*_tel.bedgraph')    

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
    #h = 'R1'
    #s = 'SL362'
    #wt, bg = get_files(table, WT, h, s, True)
    #draw_wt_v_bg(wt, bg, h, s, strain_map)
    
    ### README
    ### This chunk generates plots a la dumesic + homer

    REGIONS  = ['cen','tel']
    STRAINS = ['CK2308','CK6644','SL362','SL363','SL364','SL365',
        'SL366','SL367','SL368','SL369','SL370','SL371','SL372','SL3']
    HIST = ['R1','R2']    
    spike = {h:{s:{r:{} for r in REGIONS} for s in STRAINS} for h in HIST} # Spike in control
    spike = populate_dict(cen_files+tel_files, spike)
    #no_spike = {h:{s:{r:{} for r in REGIONS} for s in STRAINS} for h in HIST}   # input normalization
    #no_spike = populate_dict(no_spike_cen_files+no_spike_tel_files, no_spike)
    #spike_input_norm = {h:{s:{r:{} for r in REGIONS} for s in STRAINS} for h in HIST} # Spike in  + input normalization
    labels = {'tel':'Distance from Telomere (kb)',
          'cen':'Distance from Centromere (kb)'}

    settings = {'SLHM20_21':{'ncol':4, 'nrow':2, 'figsize':(20,10)},
                'SLHM17_18':{'ncol':4, 'nrow':2, 'figsize':(20,10)},
                'SLHM23':   {'ncol':4, 'nrow':3, 'figsize':(20,15)},
                'SLHM19':   {'ncol':4, 'nrow':3, 'figsize':(20,15)},
                'SLHM25_26':{'ncol':5, 'nrow':3, 'figsize':(20,15)}}

    yline_dict = {h:{s:{r:{} for r in REGIONS} for s in STRAINS} for h in HIST}
    for r in REGIONS:
        for h in HIST:
            yline_dict[h]['SL372'][r] = spike[h]['SL372'][r].subtracted.mean()

    # for r in REGIONS:
    for r in REGIONS:
        for h in HIST:
            fig = plt.figure(figsize=settings['SLHM25_26']['figsize'])
            for i,s in enumerate(tqdm(STRAINS)):
                indiv_plot_dir = PLOT_DIR + 'individual_plots/bg_sums/'
                isExist = os.path.exists(indiv_plot_dir)
                if not isExist:
                    os.makedirs(indiv_plot_dir)
                indiv_plot_name = f'{indiv_plot_dir}{r}_{h}_{s}.svg'
                #print(h,s,r)
                bg = spike[h][s][r]
                if len(bg)==0: 
                    continue
                
                xlim, ylim = find_limits2(h,r, REGIONS, HIST, STRAINS, spike)
                fig = smoothplot(bg, i+1, strain_map, r, yline_dict[h]['SL372'][r], 
                                 fig=fig, x_lim=xlim, y_lim=ylim, 
                                 title=s, 
                                 y_label='ChIP signal', x_label=labels[r], 
                                 x_metric='id',
                                 nrow=settings[kind]['nrow'], ncol=settings[kind]['ncol'], 
                                 save_individual_plots=indiv_plot_name)

            title = h+ '\nNormalized Bedgraph Sums'
            f  = PLOT_DIR+ 'bg_sums_input_'+h+'_'+r+'_'+title.split('\n')[1].replace(' ','_')
            fig.suptitle(title)
            plt.tight_layout()
            fig.subplots_adjust(top=0.9)

            fig.savefig(f + '.svg') 
            fig.savefig(f + '.png')
    
def populate_dict(files, d):
    #change me to directory for output from get_background.sh
    background_file = '/home/manning/scripts/sujin_testing/temp_holding/bg_averages.txt'
    background = pd.read_csv(background_file, sep=' ').dropna()
    #apply passes each filename in background['file'] to which_hist/get_strain
    background['hist'] = background.file.apply(lambda x: which_hist(x))
    background['strain'] = background.file.apply(lambda x: get_strain(x))

    for f in files:
        h = which_hist(f)
        r = cen_or_tel(f)
        s = get_strain(f)
        bg = readin_bedgraph(f, background)
        d[h][s][r] = bg
    return d

def readin_bedgraph(f, background, spikein=False, for_supplements=False):
    bg = pd.read_csv(f, names=['chr','start','end','id','sum_'], sep='\t').reset_index()

    #determine if a bedgraph is cen or tel file
    tel_is_true = re.findall('cen|tel', f)[0] == 'tel'
    if tel_is_true:
        #get last character of left split from bg_id, (L or R)
        bg['side'] = bg.id.apply(lambda x: x.split('_')[0][-1])
        #get an id field, i.e. 1L_1<, 1L_2<
        bg['id'] = bg.id.apply(lambda x: x.split('_')[1])
    else:
        bg['side'] = 'C'

    bg['id']    = pd.to_numeric(bg['id'], errors='coerce')
    bg['sum_']  = pd.to_numeric(bg['sum_'], errors='coerce')
    bg['fname'] = f
    h = which_hist(bg.fname[0])
    s = get_strain(bg.fname[0])
    bg['region'] = bg.fname.apply(lambda x: cen_or_tel(x))

    #background file has 2 fields at start, "file" and "average". "file" was split into "hist" and "strain"
    bg['background'] = float(background[(background['hist']==h) & (background['strain']==s)]['average'])

    bg['length'] = bg['end'] - bg['start']
    bg['average_coverage'] = bg['sum_'] / bg['length']
    bg['subtracted'] = bg['average_coverage'] - bg['background']
    #divide vs subtract background
    #bg['subtracted'] = bg['average_coverage'] / bg['background']
    bg['log_fold_change'] = np.log(bg['average_coverage'] / bg['background'])
    bg['hist'] = h
    bg['strain'] = s
    bg['chr_n'] = bg['chr'].apply(lambda x: int(x.split('chr')[1]))

    bin_group = bg.groupby('id').mean().reset_index()
    if for_supplements==True:
        return bg
    return bin_group

def which_hist(f):
    mini = f.split('_S')[0]
    hist = re.findall('WCE|R1|R2', mini)[0]
    return hist

def get_strain(f):
    # change me, listing all different strains in the experiment 
    STRAINS = ['CK2308','CK6644','SL362','SL363','SL364','SL365',
           'SL366','SL367','SL368','SL369','SL370','SL371','SL372','SL3']
    strn = re.findall('_|'.join(STRAINS)+'_', f)[0]
    strn = strn[:-1] # remove trailing _
    return strn

#ensure that text "cen" or "tel" does not appear in path to files!
def cen_or_tel(f):
    c_or_t = re.findall('cen|tel', f)[0]
    return c_or_t

def make_table(f):    
    """Returns table. f is either cen_files or tel_files"""
    table = pd.DataFrame(columns=['fname','hist','cen_tel','strain'])
    table['fname'] = f
    table['hist'] = table['fname'].apply(lambda x: which_hist(x))
    table['cen_tel'] = table['fname'].apply(lambda x: cen_or_tel(x))
    table['strain'] = table['fname'].apply(lambda x: get_strain(x))
    return table

def get_files(table, WT, h, s, for_supp = False):
    file_wt = table[(table['strain']==WT) & (table['hist']==h)]['fname'].values[0]
    file_bg = table[(table['strain']==s) & (table['hist']==h)]['fname'].values[0]

    wt = readin_bedgraph(file_wt, for_supplements = for_supp)
    bg = readin_bedgraph(file_bg, for_supplements = for_supp)
    return wt, bg

def make_spline(x,y):
    """Returns smoothed x,y"""
    x_ = np.linspace(x.min(), x.max(), 300) 
    spl = make_interp_spline(x, y, k=3)  # type: BSpline
    y_ = spl(x_)
    return(x_,y_)

#bg specified by a specific strain title
def draw_wt_v_bg(wt, bg, h, s, strain_map):
    fig = plt.figure(figsize=[10,30])
    i=1
    color='grey'
    color2='red'
    ncol=2
    nrow=14
    alpha=.5

    for c in set(wt['chr_n']):
        for si in set(wt['side']):
            
            ax = fig.add_subplot(nrow,ncol,i)
            
            g = wt[(wt['chr_n']==c) & (wt['side']==si)]
            g = g.sort_values(by='id')
            x = g['id']
            y = g['subtracted']
            y = y.fillna(0)
            x,y = make_spline(x,y)

            #plt.plot(x,y,color=color, alpha=.5)
            plt.fill_between(x,y,color=color, alpha=.3)
            plt.gca().set_xlim(left=1)
            plt.gca().set_ylim(bottom=0)

            g = bg[(bg['chr_n']==c) & (bg['side']==si)]
            g = g.sort_values(by='id')
            x = g['id']
            y = g['subtracted']
            y = y.fillna(0)
            x,y = make_spline(x,y)

            #plt.plot(x,y,color=color2, alpha=.5)
            plt.fill_between(x,y,color=color2, alpha=.3)
            plt.gca().set_xlim(left=1)
            plt.gca().set_ylim(bottom=0)
            
            y2 = max(max(wt.subtracted), max(bg.subtracted))*1.1
            x2 = max(max(wt.id), max(bg.id))
            x1 = min(min(wt.id), min(bg.id))
            
            plt.xlim([x1,x2])
            plt.ylim([0,y2])

            
            legend_properties = {'weight':'bold'}
            #label = ''.join(str(x) for x in n)
            label = ''.join(str(x) for x in [c,si])
            legend_elements = [Patch(facecolor='white', edgecolor='white',   label=label)]
            plt.legend(handles=legend_elements, loc='upper right', prop=legend_properties, 
                      frameon=False)

            i+=1     
    s = strain_map[s]
    fig.suptitle(f'{h} {s}')
    fig.subplots_adjust(top=0.96)
    plt.savefig('temp.png')
    return fig

def find_limits2(hist, region, REGIONS, HIST, STRAINS, spike, x_metric='id', y_metric='subtracted'):
    y_mins, y_maxs, x_mins, x_maxs = [], [], [], []
    for r in REGIONS:
        if r == region:
            for h in HIST:
                for s in STRAINS:
                    bg = spike[h][s][r]
#                     print(bg.columns)
                    if len(bg)==0: continue
                    y = bg[y_metric]
                    y_mins.append(y.min())
                    y_maxs.append(y.max())
                
    for s in STRAINS:
        bg = spike[hist][s][region]
        if len(bg)==0: continue
        x = bg[x_metric]
        x_mins.append(x.min())
        x_maxs.append(x.max())
        
    x_mins = [1] if region=='tel' else x_mins
#     return ( (min(x_mins), max(x_maxs)), (min(y_mins), max(y_maxs)) )
    return ( (min(x_mins), max(x_maxs)), (0, max(y_maxs)*1.1) )

def smoothplot(bin_group, i, strain_map, region, y_line, x_metric='id', y_metric='subtracted',
               x_lim=None, y_lim=None, fig=None, ncol=4, nrow=2,
               x_label='Insert x label', y_label='Insert y label',
               y_scaling_factor=100000, title='Insert title', save_individual_plots=None):    

    x = bin_group[x_metric]
    y = bin_group[y_metric] #/ y_scaling_factor
    #print(x)
    #print(y)

    x_, y_ = make_spline(x,y)
    y_[y_<0] = 0

    if x_lim==None: x_lim = [min(x), max(x)]
    if y_lim==None: y_lim = [min(y), max(y)]
    else: y_lim = y_lim #np.divide(y_lim, y_scaling_factor, dtype=float)

    ax = fig.add_subplot(nrow,ncol,i)
    color = 'grey'
    plt.plot(x_, y_, color=color)
    plt.fill_between(x_, y_, color=color)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if region == 'cen':
        plt.yticks(np.arange(min(y_lim), 5, 0.5))
        plt.xlim(x_lim)
        plt.ylim([0,5])
        plt.hlines(y= y_line, xmin=min(x_lim), xmax=max(x_lim), colors='purple', linestyles='--', lw=1)
    else:
        #ticks = np.insert(np.arange(5, max(y_lim)+5, 5.0), 0, 1.)
        #for use with subtraction
        ticks = np.arange(0, max(y_lim)+5, 5)
        plt.yticks(ticks)        
        plt.xlim(x_lim)
        plt.ylim(y_lim)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))


   
    legend_properties = {'weight':'bold'}
    label = strain_map[title]
    legend_elements = [Patch(facecolor='white', edgecolor='white',   label=label)]
    plt.legend(handles=legend_elements, loc='upper right', prop=legend_properties, 
              frameon=False)

    if save_individual_plots!=None:
        plt.savefig(save_individual_plots)
    
    return(fig)

if __name__ == '__main__':
    execute()