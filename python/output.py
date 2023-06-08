import time
import pandas as pd
import collections
import seaborn as sns
import matplotlib.pyplot as plt


def write_output(result):
    result = pd.DataFrame.from_dict(result, orient='index')
    result = result.reset_index().add_prefix('data')
    result.rename(columns=
                  {'dataindex': 'gene_symbol', 'data0': 'gene_id',
                   'data1': 'gene_short_name', 'data2': 'gene_full_name', 'data3': 'protein_id',
                   'data4': 'protein_function', 'data5': 'go_bp', 'data6': 'go_mf', 'data7': 'go_cc',
                   'data8': 'Sc_prot_id', 'data9': 'sc_prot_name', 'data10': 'query_cov',
                   'data11': 'e_value', 'data12': 'perc_identity', 'data13': 'acc_len'}, inplace=True)
    result.to_csv('result.csv')


def create_images():
    pd.options.display.max_colwidth = 2000
    df = pd.read_csv('result.csv')
    for go_name in ('go_bp', 'go_mf', 'go_cc'):
        go = df[go_name]
        go = go.dropna().to_string(index=False)
        go = go.split('\n')
        go = [i.split('\\n') for i in go]
        go = [item for inner_list in go for item in inner_list]
        go = [i.split('[')[0].strip() for i in go]
        go = list(filter(None, go))

        counter = collections.Counter(go)
        dic = dict(counter)

        res_go = pd.DataFrame.from_dict(dic, orient='index').add_prefix('data')
        res_go = res_go.reset_index()
        res_go.rename(columns={'index': 'go', 'data0': 'count'}, inplace=True)
        res_go = res_go.sort_values(by='count', ascending=False)
        res_go.to_csv(f'table_counts_for_{go_name}.csv')

        res_go = res_go[res_go['count'] > 1]
        if not res_go.empty:
            res_go = res_go[:13]
            plt.figure(figsize=(9, 7))
            ax = sns.barplot(x="count", y="go", data=res_go, palette='Set2', saturation=0.7, width=0.9)
            ax.set_title(f'{go_name}', weight='bold', fontsize=14)
            ax.tick_params(labelsize=10)
            sns.set()

            plt.savefig(f'{go_name}.png')
            time.sleep(1)
