"""
script to combine PPA and NPA accuracy tables
"""
import pandas as pd
import sys

bip_1_out_dir = sys.argv[1]
version1 = sys.argv[2]
bip_2_out_dir = sys.argv[3]
version2 = sys.argv[4]
out_dir = "/ghds/groups/bioinformatics/02_DEVELOPMENT/200623_BIP_VALIDATION_PIPLINES/accuracy_comparison/"

fileName = out_dir + 'accuracy1_bip_' + version1 + '_accuracy2_bip_' + version2 + '.pdf'

PPA_tab_f = "ppa_combined4.tsv"
NPA_tab_f = "npa_combined4.tsv"
ppa_colnames = ['variant category', 'variant type', 'acceptable LLCI PPA', 'CDx+ LBP70+', 'CDx- LBP70+', 'PPA', 'LLCI', 'ULCI']
npa_colnames = ['variant category', 'variant type', 'acceptable LLCI NPA', 'CDx+ LBP70-', 'CDx- LBP70-', 'NPA', 'LLCI', 'ULCI']

ppa_tab_df1 = pd.read_csv(bip_1_out_dir + PPA_tab_f, sep = '\t')
ppa_tab_df2 = pd.read_csv(bip_2_out_dir + PPA_tab_f, sep = '\t')
ppa_tab_df1.columns = ppa_colnames
ppa_tab_df2.columns = ppa_colnames

npa_tab_df1 = pd.read_csv(bip_1_out_dir + NPA_tab_f, sep = '\t')
npa_tab_df2 = pd.read_csv(bip_2_out_dir + NPA_tab_f, sep = '\t')
npa_tab_df1.columns = npa_colnames
npa_tab_df2.columns = npa_colnames

def roundVals(col):
    return [round(float(item), 1) for item in col]

#######################
##  Merged Accuracy PPA
#######################

merged_PPA_cols = ['variant category', 'variant type', 'acceptable PPA', 'bip:' + version1+ 'CDx+/LBP70+', 'bip:' + version2+ 'CDx+/LBP70+', 'bip:' + version1+ 'CDx-/LBP70+', 'bip:' + version2+ 'CDx-/LBP70+', 'bip:' + version1+ ' PPA', 'bip:' + version2+ ' PPA', 'LLCI', 'ULCI']
PPA_merged_df = pd.DataFrame(columns = merged_PPA_cols)

PPA_merged_df['variant category'], PPA_merged_df['variant type'], PPA_merged_df['acceptable PPA'] = ppa_tab_df1['variant type'], ppa_tab_df1['variant category'], ppa_tab_df1['acceptable LLCI PPA']
PPA_merged_df['bip:' + version1+ 'CDx+/LBP70+'], PPA_merged_df['bip:' + version2+ 'CDx+/LBP70+'] = ppa_tab_df1['CDx+ LBP70+'], ppa_tab_df2['CDx+ LBP70+']
PPA_merged_df['bip:' + version1+ 'CDx-/LBP70+'], PPA_merged_df['bip:' + version2+ 'CDx-/LBP70+'] = ppa_tab_df1['CDx- LBP70+'], ppa_tab_df2['CDx- LBP70+']
PPA_merged_df['bip:' + version1+ ' PPA'], PPA_merged_df['bip:' + version2+ ' PPA'] = roundVals(ppa_tab_df1['PPA']), roundVals(ppa_tab_df2['PPA'])
PPA_merged_df['LLCI'], PPA_merged_df['ULCI']  = roundVals(ppa_tab_df1['LLCI']), roundVals(ppa_tab_df2['ULCI'])

PPA_merged_df.to_csv(out_dir + 'accuracy_PPA_merged_table_1_bip_'+version1+'_2_bip_'+version2+'.tsv', sep = '\t', index = False)

#######################
##  Merged Accuracy NPA
#######################


merged_NPA_cols = ['variant category', 'variant type', 'acceptable NPA', 'bip:' + version1+ 'CDx+/LBP70-', 'bip:' + version2+ 'CDx+/LBP70-', 'bip:' + version1+ 'CDx-/LBP70-', 'bip:' + version2+ 'CDx-/LBP70-', 'bip:' + version1+ ' NPA', 'bip:' + version2+ ' NPA', 'LLCI', 'ULCI']
NPA_merged_df = pd.DataFrame(columns = merged_NPA_cols)

NPA_merged_df['variant category'], NPA_merged_df['variant type'], NPA_merged_df['acceptable NPA'] = npa_tab_df1['variant type'], npa_tab_df1['variant category'], npa_tab_df1['acceptable LLCI NPA']
NPA_merged_df['bip:' + version1+ 'CDx+/LBP70-'], NPA_merged_df['bip:' + version2+ 'CDx+/LBP70-'] = npa_tab_df1['CDx+ LBP70-'], npa_tab_df2['CDx+ LBP70-']
NPA_merged_df['bip:' + version1+ 'CDx-/LBP70-'], NPA_merged_df['bip:' + version2+ 'CDx-/LBP70-'] = npa_tab_df1['CDx- LBP70-'], npa_tab_df2['CDx- LBP70-']
NPA_merged_df['bip:' + version1+ ' NPA'], NPA_merged_df['bip:' + version2+ ' NPA'] = roundVals(npa_tab_df1['NPA']), roundVals(npa_tab_df2['NPA'])
NPA_merged_df['LLCI'], NPA_merged_df['ULCI']  = roundVals(npa_tab_df1['LLCI']), roundVals(npa_tab_df2['ULCI'])

NPA_merged_df.to_csv(out_dir + 'accuracy_NPA_merged_table_1_bip_'+version1+'_2_bip_'+version2+'.tsv', sep = '\t', index = False)


################
##  build report
################

from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.pagesizes import letter
from reportlab.platypus import Table
from reportlab.platypus import TableStyle
from reportlab.lib import colors
from reportlab.lib.units import mm, inch
from reportlab.lib.styles import getSampleStyleSheet



styles = getSampleStyleSheet()
styleN = styles['Heading1']

pagesize = (600 * mm, 400 * mm)

pdf = SimpleDocTemplate(
    fileName,
    pagesize=pagesize
)



def construct_table(data):
    
    table = Table(data)

    # add style
    style = TableStyle([
        #('BACKGROUND', (0,0), (3,0), colors.green),
        ('TEXTCOLOR',(0,0),(-1,0),colors.black),

        ('ALIGN',(0,0),(-1,-1),'CENTER'),

        ('FONTNAME', (0,0), (-1,0), 'Courier'),
        ('FONTSIZE', (0,0), (-1,0), 12),

        ('BOTTOMPADDING', (0,0), (-1,0), 12),

        ('BACKGROUND',(0,1),(-1,-1),colors.beige),
    ])
    table.setStyle(style)


    # 3) Add borders
    ts = TableStyle(
        [
        ('BOX',(0,0),(-1,-1),2,colors.black),

        #('LINEBEFORE',(2,1),(2,-1),2,colors.red),
        #('LINEABOVE',(0,2),(-1,2),2,colors.green),

        ('GRID',(0,0),(-1,-1),2,colors.black),
        ]
    )
    table.setStyle(ts)
    
    return table



data = list()
data.append(list(PPA_merged_df.columns))

for i, row in PPA_merged_df.iterrows():
    data.append(list(PPA_merged_df.loc[i]))
    
table = construct_table(data)

P = Paragraph("Table summarizing the comparison of accuracy (Positive Percent Agreement) across variant classes between two BIP versions", styleN)

elems = []
elems.append(P)
elems.append(table)


data = list()
data.append(list(NPA_merged_df.columns))

for i, row in NPA_merged_df.iterrows():
    data.append(list(NPA_merged_df.loc[i]))
    
table = construct_table(data)

P = Paragraph("Table summarizing the comparison of accuracy (Negative Percent Agreement) across variant classes between two BIP versions", styleN)

P2 = Spacer(1, 1.25*inch)

elems.append(P2)
elems.append(P)
elems.append(table)

pdf.build(elems)

