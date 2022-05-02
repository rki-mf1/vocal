import argparse as ap
import timeit as ti
import sys
import os
import pandas as pd
import jinja2


ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def render_template(template, **kwargs):
    ''' renders a Jinja template into HTML '''
    # check if template exists
    if not os.path.exists(template):
        print('No template file present: %s' % template)
        sys.exit()

    templateLoader = jinja2.FileSystemLoader(searchpath="/")
    templateEnv = jinja2.Environment(loader=templateLoader)
    templ = templateEnv.get_template(template)
    return templ.render(**kwargs)

def get_TOP10_clusters(vocal_alerts_clusters_summaries_all):
    cluster_df = vocal_alerts_clusters_summaries_all[vocal_alerts_clusters_summaries_all['cluster_ID_in_alert_level'].between(0, 10)]
    cluster_df= cluster_df.sort_values(["alert_level", "cluster_ID_in_alert_level"],ascending=[False,True])
    return cluster_df

def get_notECDC_PINKRED(_clusters_summaries_all,_ecdc_df ):
    '''
    Not in ECDC list and predicted with RED and PINK alert
    '''
    _result = _clusters_summaries_all[(~_clusters_summaries_all['Lineages'].isin( _ecdc_df['Pango lineage'].to_list())) & (_clusters_summaries_all['alert_level'].isin(['red', 'pink']))]
    return _result.sort_values(by=['n_samples'], ascending=False)

def get_ECDC_GREY( _samples_all, _ecdc_df):
    '''
    cluster of ECDC variants with GREY alert
    '''
    #varaint_ECDC_df  = _samples_all[_samples_all['LINEAGE.LATEST'].isin(ECDC_VOC['Pango lineage'].to_list())]
    #filtered_cluster_VOC = _clusters_summaries_all[_clusters_summaries_all['cluster_ID_in_alert_level'].isin(varaint_ECDC_df.cluster_ID_in_alert_level.to_list())]
    _samples_all = _samples_all[ (_samples_all['alert_level']== 'grey') & (_samples_all['LINEAGE'].isin(_ecdc_df['Pango lineage'].to_list())) ]
    #_samples_all = _samples_all.groupby(["cluster_ID_in_alert_level",'LINEAGE'], sort=True)["ID"].count()
    _samples_all = _samples_all.sort_values(by=['s_moc_roi_tot','nMutationsTotal_D','nMutationsTotal_M','nMutationsTotal_I','nLineageDefining_D','nLineageDefining_M','nLineageDefining_I'], ascending=False)
    return _samples_all

### Samples
def get_High_PM(_samples_all):    
    '''
    Samples that have high PM (high private mutaions, SNP > 20 or DEL > 5 or INS > 5)
    '''
    _samples_all=_samples_all[ (_samples_all['s_pm_M'] > 20 ) |(_samples_all['s_pm_I'] > 5 ) | (_samples_all['s_pm_D'] > 5 )]
    return _samples_all.sort_values(by=['s_pm_M','s_pm_I','s_pm_D'], ascending=False)

def get_samples_VOC_GREY(_samples_all, ECDC_VOC):
    '''
    VOC consider alert level samples that have grey 
    '''
    _result = _samples_all[(_samples_all['alert_level'].isin(['grey']) ) & (_samples_all['LINEAGE'].isin(ECDC_VOC['Pango lineage'].to_list()))]
    _result =  _result.where(pd.notnull(_result), None)
    return _result

def main(args):
    ## Arguments
    _clusters_df = pd.read_csv(args.vocal_alert_clusters)
    _samples = pd.read_csv( args.vocal_alert_samples)
    _vocal_annoDB_dir = args.vocal_annoDB_dir
    _short_form = args.short_sum
    if(args.from_date == ''):
        from_date = '-'
    else:
        from_date = args.from_date

    if(args.to_date == '' ):
        to_date = '-'
    else:
        to_date = args.to_date
    output = args.out_file 
    
    # Read Program Version
    vocal_version = os.path.join( ROOT_DIR,'.version')
    _sw_dict={}
    with open(vocal_version) as f:
        for line in f:
            sn, v = line.split('=')
            _sw_dict[sn] = v 
    # Read DB Version
    vocal_DB_version = os.path.join( ROOT_DIR,'data/.version')
    _db_dict={}
    with open(vocal_DB_version) as f:
        for line in f:
            if len(line.split('=')) == 0:
                continue
            sn, v = line.split('=')
            _db_dict[sn] = v 
        # Read Annonation DB
        _ecdc_df  = pd.read_csv(os.path.join(_vocal_annoDB_dir,'ECDC_assigned_variants.csv'))
        ECDC_not_VOC = _ecdc_df[_ecdc_df['Status']!='VOC'].reset_index(drop=True)
        ECDC_VOC = _ecdc_df[_ecdc_df['Status']== 'VOC'].reset_index(drop=True)
        # _linage_sublinage_df = pd.read_csv(os.path.join( ROOT_DIR,'data/lineage.all.tsv'))

        # Preprocessing DATA 

        ## Top 10  
        _top10_df = get_TOP10_clusters(_clusters_df)[['alert_level','n_samples','Lineages','first_seen_isolate','cluster_ID_in_alert_level','s_moc_M.avg','s_moc_D.avg','s_moc_I.avg','ListFrequentMutations_gt30perc']]
        ## All cluster
        _all_cluster = _clusters_df[['alert_level','n_samples','Lineages','first_seen_isolate','last_seen_isolate','date_range','cluster_ID_in_alert_level']]
        ## cluster not in ECDC but got predition  in PINK and RED
        _ECDC_notVOC_PINKRED = get_notECDC_PINKRED(_clusters_df,_ecdc_df)[['alert_level','n_samples','first_seen_isolate', 'Lineages','cluster_ID_in_alert_level']]
        
        ## New cluster today 
        ## -- soon --

        ## samples of ECDC variants with GREY alert.
        _ECDC_GREY = get_ECDC_GREY(_samples, _ecdc_df)[['ID','VariantType','nMutationsTotal_D','nMutationsTotal_M','nMutationsTotal_I','nLineageDefining_D','nLineageDefining_M','nLineageDefining_I', 'LINEAGE']]##[['alert_level','n_samples','first_seen_isolate', 'Lineages','cluster_ID_in_alert_level']]
        #_ECDC_GREY["ID"] = pd.util.hash_array(_ECDC_GREY["ID"].to_numpy())

        ## samples have high mutation. 
        _sample_high_mut = get_High_PM(_samples)[['ID','alert_level','s_pm_M','s_pm_D','s_pm_I', 'LINEAGE','VariantType']]
        #_sample_high_mut["ID"] = pd.util.hash_array(_sample_high_mut["ID"].to_numpy())
        
        ## VOC samples with GREY alert
        #_cluster_VOC_grey = get_samples_VOC_GREY(_samples, ECDC_VOC)[['ID','ListMutationsSelected_D','ListMutationsSelected_M','ListMutationsSelected_I','LINEAGE']]
    if(not _short_form):

        # Read the template file 
        html_template = os.path.join( ROOT_DIR,'email-template/email.html')
        #if not os.path.exists(tmp_dir):
        #    os.makedirs(tmp_dir)
    
        # clean tmp directory
        #try:
        #    shutil.rmtree(tmp_dir) 
        #except OSError as e:
        #    print("Error: %s - %s." % (e.filename, e.strerror))
    else:
        _top10_df.drop('ListFrequentMutations_gt30perc', axis=1, inplace=True) 
        # list only top 10 
        if(len( _ECDC_notVOC_PINKRED) > 10):
           _ECDC_notVOC_PINKRED = _ECDC_notVOC_PINKRED[0:10]

        if(len(_sample_high_mut) > 10):
            _sample_high_mut = _sample_high_mut[0:10]
        if(len(_ECDC_GREY) > 10):
            _ECDC_GREY = _ECDC_GREY[0:10]

        html_template = os.path.join( ROOT_DIR,'email-template/email.sum.html')
    
    # generate HTML from template
    html = render_template(html_template, **locals()) # locals means using all local var in this program

    # to save the results
    html = render_template(html_template, **locals()) # locals means using all local var in this program
    with open(output, "w") as fh:
        fh.write(html)


if __name__ == '__main__':

    parser = ap.ArgumentParser(description='Report Generator',
                                 epilog=''' Please let us know if there is any problem at https://github.com/rki-mf1/vocal/
                                        ''')
    parser.add_argument('-s', '--vocal_alert_samples',
                        help='vocal_alerts_samples_all.csv file',
                        required=True)

    parser.add_argument('-c', '--vocal_alert_clusters', 
                        help='vocal_alerts_clusters_summaries_all.csv file',
                        required=True)
                        # lambda s: datetime.datetime.strptime(s, '%Y-%m-%d')
    parser.add_argument('-a', '--vocal_annoDB_dir', 
                        default= ROOT_DIR+"/data",)
    
    parser.add_argument('--from_date', type=str,
                        default= '',
                        help='start date (format YYYY-MM-DD)')

    parser.add_argument('--to_date', type=str,
                        default= '',
                        help='end date (format YYYY-MM-DD)')

    parser.add_argument('-S', '--short_sum', 
                        help="Generate report in short form",
                        action="store_true")

    parser.add_argument('-o', '--out_file',  help='Ouput in HTML format',
                        required=True,)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()
    # print(args)
    t1 = ti.default_timer()
    main(args)
    t2 = ti.default_timer()
    print('Processing time: {} seconds'.format(round(t2 - t1), 2))