#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
#import lxml, html5lib
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import dash_table

from zipfile import ZipFile
from glob import glob

from time import sleep


# In[ ]:


batches = glob("Batches/*.zip")
availableBatches = [batch for batch in batches]

omim_map = pd.read_csv('assets/OMIM_map.txt', sep = '\t')

MGI_IDs = pd.read_csv('assets/MGI_IDs.csv', dtype = str)
MGI_IDs = dict(zip(MGI_IDs['Symbol'], MGI_IDs['MGI_ID']))

hpo = pd.read_csv('assets/HPO_dict.txt', sep = '\t')
hpo = hpo[hpo['Associated Genes'].str.len() > 1]

refSeqSummaries = pd.read_csv('assets/RefSeqSummaries.tsv', sep = '\t', dtype = str)

default_style = {'width' : '57%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}
hide = {'display' : 'None'}
full_display = {'width' : '100%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'textAlign' : 'center'}
info_style = {'width' : '100%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'left'}
break_line = html.Hr(style={'height' : '4px', 'width' : '60%', 'color' : '#111111','display' : 'inline-block', 'marginLeft':'auto', 'marginRight':'auto'})


# In[ ]:


def generateHeatTable(gene, Chromosome, Position):
    sleep(2)
    children = []
    try:
        url = 'https://www.ncbi.nlm.nih.gov/clinvar?term=(' + gene + '%5BGene%20Name%5D)%20AND%20pathogenic'
        df = pd.read_html(url)[0]
        if df.shape[0] > 100:
            children.append(html.P('Too many pathogenic variants for this gene!😵'))
        elif df.shape[0] == 0:
            children.append(html.P('Could not find pathogenic variants for this gene in ClinVar🤷 🤷‍♂️'))
        else:
            del df['VariationLocation']
            del df['Gene(s)']
            del df['Clinical significance (Last reviewed)']
            del df['Review status']
            del df['Protein change']
            df.fillna('', inplace = True)
            def getGRCh38(value):
                try:
                    location = value.split('GRCh38: Chr')[1]
                    if '-' in location:
                        location = location.split('-')[0]
                    return location
                except:
                    return 'delete'
            def getVariation(value):
                try:
                    variant = value.replace('GRCh38/hg38 ','').split('GRC')[0]
                    return variant
                except:
                    return 'delete'
            df['Coordinates'] = df['VariationLocation.1'].apply(getGRCh38)
            df['Variation'] = df['VariationLocation.1'].apply(getVariation)
            del df['VariationLocation.1']
            df = df[df['Coordinates'] != 'delete']
            df = df[df['Variation'] != 'delete']
            df['Condition(s)'] = df['Condition(s)'].str.replace('See cases', 'Check link ➜')
            df['Position'] = df['Coordinates'].str.split(':').str[1]
            df['Position'] = df['Position'].astype(int)
            coordinates = str(Chromosome) + ':' + str(Position)
            if coordinates in df['Coordinates']:
                False
            else:
                df = df.append({
                    'Coordinates' : coordinates,
                    'Position' : int(Position),
                    'Condition(s)' : '👨‍🔬👩‍🔬',
                    'Accession' : 'Your variant',
                    'Variation' : '👩‍⚕️ 👨‍⚕️'
                },
            ignore_index=True)
            df.sort_values('Position', inplace = True)
            if df['Position'].max() - int(Position) >= int(Position) - df['Position'].min():
                maxSpace = df['Position'].max() - int(Position)
            else:
                maxSpace = int(Position) - df['Position'].min()
            df['gradient'] = (1 - abs(int(Position) - df['Position']) / maxSpace) * 255 ** 0.7
            df['gradient'] = df['gradient'].astype(int)
            df.at[df.index[df['Coordinates'] == coordinates].tolist(), 'gradient'] = 255
            df['gradient'] = df['gradient'].apply(lambda x : 'rgb(255,' + str(255 - x) + ',' + str(255 - x)+')')
            del df['Position']
            df = df[['Coordinates','Condition(s)','Accession','Variation', 'gradient']]
            cols = ('Coordinates','Condition(s)','Accession','Variation')
            df['link'] = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + df['Accession']
            print(df)
            children.append(html.Table(
                [html.Tr([html.Th(col) for col in cols], style = {'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'})] +
                [html.Tr([
                    html.Td(df.iloc[i][col], style = {'backgroundColor' : df.iloc[i]['gradient'], 'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'}) if col != 'Accession' else html.Td(html.A(href=df.iloc[i]['link'], children=df.iloc[i][col], target='_blank'), style = {'backgroundColor' : df.iloc[i]['gradient'], 'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'}) for col in cols
                ]) for i in range(len(df))],
                style = {'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'}
            )
        )
    except:
        children.append(html.P('Could not find pathogenic variants for this gene in ClinVar🤷 🤷‍♂️'))
    return children


# In[ ]:


external_stylesheets = ['assets/BAMdelbee_main.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "BAMdelbee | Ohad Birk's Lab"
app.config['suppress_callback_exceptions'] = False


# In[ ]:


app.layout = html.Div([
    html.Div([
        html.Img(id = 'header_image',
            src = 'assets/logo_birklab.png',
            style = {
                'height' : '25%',
                'width' : '25%',
                'float' : 'left'
            }
        ),
        html.Span(''),
        html.Img(id = 'header_image2',
            src = 'assets/logo_BAMdelbee.png',
            style = {
                'height' : '20%',
                'width' : '20%',
                'float' : 'right',
                'marginLeft' : 4
            }
        ),
    ]),
    html.Span(),
    html.H2(
        'Select a batch',
        id = 'instructions',
        style = {'width' : '100%', 'height' : '100%', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}

    ),
    dcc.Dropdown(
        id = 'BatchesDropdown',
        options = [{'label': v[8:-4], 'value': v} for v in batches],
        multi = False,
        placeholder = 'Select batch',
        style = default_style
    ),
    dcc.RadioItems(
        id = 'referenceRadioButtons',
        options=[
            {'label': 'hg38', 'value': 'hg38'},
            {'label': 'hg19', 'value': 'hg19'},
        ],
        value='hg38',
        labelStyle={'display': 'inline-block'},
        style = default_style
    ),
    html.Div([
        html.Div([
            html.Div(id = 'checkListContainer'),
        ], id = 'checklist_div',
        style = {'display' : 'None'}
        )],
        id = 'parameters_div',
        style = default_style
    ),
    html.Div(
        dcc.Loading(
            id = 'loading_results',
            fullscreen = True,
            children = [html.Div([html.Div(id = 'loading_table_results')])],
            type = 'circle',
            style = {'marginTop' : '120px'}
        ),
        id = 'loading_results_div',
        style={'width': '100%', 'height': '10%', 'borderWidth': '1px', 'borderStyle': 'grooved', 'borderRadius': '5px', 'textAlign': 'center', 'marginLeft':'auto', 'marginRight':'auto'}
    ),
    html.Div(
        dcc.Loading(
            id = 'loading_heatTable',
            fullscreen = True,
            children = [html.Div([html.Div(id = 'loading_heattable_results')])],
            type = 'circle',
            style = {'marginTop' : '120px'}
        ),
        id = 'loading_heatTable_div',
        style={'width': '100%', 'height': '10%', 'borderWidth': '1px', 'borderStyle': 'grooved', 'borderRadius': '5px', 'textAlign': 'center', 'marginLeft':'auto', 'marginRight':'auto'}
    ),
    html.Div(id = 'table_div', style = hide)
])


# In[ ]:


@app.callback(
    [Output('BatchesDropdown', 'style'),
     Output('referenceRadioButtons', 'style'),
     Output('checklist_div', 'style'),
     Output('instructions', 'children'),
     Output('checkListContainer', 'children')],
    [Input('BatchesDropdown', 'value')]
)
def chooseBatch(batch):
    if batch != None:
        zBatch = ZipFile(batch)
        samples = zBatch.read('results/samples.txt').decode('utf-8', 'ignore').split('\n')
        samples = [{'label' : str(sample), 'value' : str(sample)} for sample in samples]
        checkList = dcc.Checklist(id = 'samplesCheckList', options = samples, labelStyle={'display': 'block'}) 
        btn = html.Button('Analyze', id = 'Analysis_starting_button', disabled = True, n_clicks = 0, style = {'text-transform' : 'none', 'font-size' : '30px'})
        return [hide, hide, default_style, batch[8:-4], [checkList, btn]]
    else:
        return [default_style, default_style, hide, "Select batch:", None]

@app.callback(
    [Output('Analysis_starting_button', 'disabled')],
    [Input('samplesCheckList', 'value')]
)
def chooseSamples(samples):
    if samples not in [None, []]:
        return [False]
    else:
        return [True]

@app.callback(
    [Output('parameters_div', 'style'),
     Output('table_div', 'style'),
     Output('table_div', 'children'),
     Output('loading_results_div', 'children')],
    [Input('Analysis_starting_button', 'n_clicks'),
     Input('BatchesDropdown', 'value'),
     Input('samplesCheckList', 'value'),
     Input('referenceRadioButtons', 'value')]
)
def update_figure(clicks, batch, affected, referenceGenome):
    if clicks in (0, None):
        return [default_style, hide, [], None]
    else:
        if referenceGenome == 'hg38':
            geneMap = pd.read_csv('assets/NCBI_Genes_hg38.map')
        else:
            geneMap = pd.read_csv('assets/NCBI_Genes_hg19.map')
        geneMap['Chr'] = geneMap['Chr'].astype(str)
        
        zBatch = ZipFile(batch)
        deletions = []
        sizes = []
        genes = []
        
        for chromosome in [str(c) for c in range(1,23)] + ['X', 'Y']:
            try:
                df = pd.read_csv(zBatch.open('results/chr' + chromosome + '.coverage'), sep = '\t')
            except:
                continue
            samples = df.columns[1:]
            controls = [sample for sample in samples if sample not in affected]
            for sample in affected:
                df = df[df[sample] == 0]
            for control in controls:
                df = df[df[control] > 5]
            if df.shape[0] < 2:
                continue
            coordinates = df['Coordinate'].tolist()
            while len(coordinates) != 0:
                for coordinate in coordinates:
                    start = coordinate
                    end = coordinate
                    while end + 50 in coordinates:
                        end = end + 50
                    coordinates = [coordinate for coordinate in coordinates if coordinate > end]
                    if end != start:
                        deletions.append(chromosome + ':' + str(start) + '-' + str(end))
                        sizes.append(end-start)
                        genes_in_deletion = geneMap[(geneMap['Chr'] == chromosome) & (geneMap['Start'] <= start) & (geneMap['End'] >= end)]['Gene'].tolist()
                        genes.append(', '.join(genes_in_deletion))
        df = pd.DataFrame()
        df['Coordinates'] = deletions
        df['Chromosome'] = df['Coordinates'].apply(lambda x : int(x.split(':')[0].upper().replace('X', '23').replace('Y', '24').replace('M', '25').replace('MT', '25')))
        df['~Size (bp)'] = sizes
        df['Genes'] = genes
        df['Start'] = df['Coordinates'].apply(lambda x : int(x.split(':')[1].split('-')[0]))
        df.sort_values(['Chromosome', 'Start'], inplace = True)
        df.reset_index(drop = True, inplace = True)
        df['End'] = df['Coordinates'].apply(lambda x : int(x.split(':')[1].split('-')[1]))
        del df['Chromosome']
        del df['Start']
        del df['End']
        df = df[df['~Size (bp)'] > 50]

        ddt = dash_table.DataTable(
            id = 'table',
            data = df.to_dict('records'),
            columns = [{'id': c, 'name': c} for c in df.columns],
            page_action = 'none',
            row_selectable='single',
            fixed_rows={'headers' : True},
            style_cell = {
                'textAlign': 'left',
                'overflow': 'hidden',
                'textOverflow': 'ellipsis',
                'maxWidth': 0
            },
            tooltip_delay=0,
            tooltip_duration=None,

            style_data_conditional=[
                {
                'if': {'row_index': 'odd'},
                'backgroundColor': 'rgb(247, 247, 255)'
                },
                {
                'if': {'row_index': 'even'},
                'backgroundColor': 'rgb(252, 252, 255)'
                },
            ],
            style_data={
                'whiteSpace': 'normal',
                'height': 'auto',
                },
            style_header = {
                'backgroundColor': 'rgb(220, 220, 255)',
                'fontWeight': 'bold',
                'whiteSpace' : 'normal'
                },
            style_table={
                'maxHeight': '250px',
                'maxwidth' : '150%',
                'border': 'thin lightgrey solid'
            },
            style_as_list_view = False
        )
        info = html.Div('Click on a deletion for more information', id = 'info')
        return [hide, full_display, [ddt, info], None]

@app.callback(
    [Output('info', 'children'),
     Output('info', 'style')],
    [Input('table', 'derived_virtual_selected_rows'),
     Input('table', 'derived_virtual_data'),
     Input('referenceRadioButtons', 'value')]
)
def getVariantData(selected_row_index, data, referenceGenome):
    if selected_row_index in [[], None]:
        return [html.P('Select row for more information'), default_style]
    else:
        elements = []
        coordinates = data[selected_row_index[0]]['Coordinates']
        genes = data[selected_row_index[0]]['Genes']
        
        #IGV & UCSC & gnomAD 
        igv_link = 'http://localhost:60151/goto?locus=chr' + coordinates
        igvA = html.A('IGV', href = igv_link, target = '_blank')
        ucsc_link = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=' + referenceGenome + '&position=' + coordinates.replace(':','%3A')
        ucscA = html.A('UCSC', href = ucsc_link, target = '_blank')
        if referenceGenome == 'hg38':
            gnomAD_Url = 'https://gnomad.broadinstitute.org/region/' + coordinates + '?dataset=gnomad_r3'
        else:
            gnomAD_Url = 'https://gnomad.broadinstitute.org/region/' + coordinates + '?dataset=gnomad_r2_1'
        gnomADLink = html.A('gnomAD',target='_blank', href = gnomAD_Url)
        elements.append(html.Div([igvA, html.Span('    '), ucscA, html.Span('    '), gnomADLink]))
        
        if len(genes) > 0:
            
            #RefSeq
            for gene in genes.split(', '):
                try:
                    refseq = refSeqSummaries[refSeqSummaries['Symbol'] == gene]
                    fullName = refseq['Full Name'].tolist()[0]
                    elements.append(html.Div(fullName, style = {'fontStyle' : 'oblique'}))
                    summary = refseq['Summary'].tolist()[0]
                    elements.append(html.Div(summary))
                    elements.append(break_line)
                except: False
        
            #HPO
            hpo_terms = hpo[hpo['Associated Genes'].str.contains('|'.join(genes.split(', ')))]['HPO'].tolist()
            phenotypes = ', '.join(hpo_terms)
            if len(phenotypes.strip()) > 0:
                elements.append(html.Div('HPO associated phenotypes:', style = {'fontStyle' : 'oblique'}))
                elements.append(html.Div(phenotypes))
                elements.append(break_line)
        
            #GeneCards
            geneCards = [html.Span('GeneCards: ')]
            for gene in genes.split(', '):
                geneCards_link = 'https://www.genecards.org/cgi-bin/carddisp.pl?gene=' + gene
                geneCardsA = html.A(gene, href = geneCards_link, target = '_blank')
                geneCards.append(html.Span('    '))
                geneCards.append(geneCardsA)
            elements.append(html.Div(geneCards))
        
            #OMIM
            omims = [html.Span('OMIM: ')]
            for gene in genes.split(', '):
                try:
                    omim = str(omim_map[omim_map['Approved Gene Symbol (HGNC)'] == gene]['MIM Number'].tolist()[0]).replace('.0','')
                    omim_link = 'https://omim.org/entry/' + omim
                    omimA = html.A(gene, href = omim_link, target = '_blank')
                    omims.append(html.Span('    '))
                    omims.append(omimA)
                except: False
            if len(omims) > 1:
                elements.append(html.Div(omims))
        
            #GTEx
            GTExs = [html.Span('GTEx: ')]
            for gene in genes.split(', '):
                GTEx_link = 'https://gtexportal.org/home/gene/' + gene
                GTExA = html.A(gene, href = GTEx_link, target = '_blank')
                GTExs.append(html.Span('    '))
                GTExs.append(GTExA)
            elements.append(html.Div(GTExs))

            #MGI
            mgis = [html.Span('MGI: ')]
            for gene in genes.split(', '):
                try:
                    mgi_link = 'http://www.informatics.jax.org/marker/phenotypes/MGI:' + MGI_IDs[gene]
                    mgiA = html.A(gene, href = mgi_link, target = '_blank')
                    mgis.append(html.Span('    '))
                    mgis.append(mgiA)
                except: False
            if len(mgis) > 1:
                elements.append(html.Div(mgis))
        
        #ClinVar HeatTable
        if len(genes) > 0:
            elements.append(html.Div(break_line))
            clinVarBtn = btn = html.Button('ClinVar HeatTable', id = 'HeatTableBtn', n_clicks = 0, style = {'text-transform' : 'none'})
            heatTableDiv = html.Div(id = 'heatTableDiv')
            elements.append(clinVarBtn)
            elements.append(heatTableDiv)
        
        return [elements, info_style]

@app.callback(
    [Output('heatTableDiv', 'children'),
     Output('loading_heatTable_div', 'children'),
     Output('HeatTableBtn', 'style')],
    [Input('HeatTableBtn', 'n_clicks'),
     Input('table', 'derived_virtual_selected_rows'),
     Input('table', 'derived_virtual_data')]
)
def update_figure(clicks, selected_row_index, data):
    if clicks in [0, None]:
        return [[None], None]
    else:
        coordinates = data[selected_row_index[0]]['Coordinates']
        gene = data[selected_row_index[0]]['Genes']
        if ', ' in gene:
            gene = gene.split(', ')[0]
        return [generateHeatTable(gene, coordinates.split(':')[0], coordinates.split(':')[1].split('-')[0]), None, {'display' : 'None'}]


# In[ ]:


if __name__ == '__main__':
    import webbrowser
    webbrowser.open('http://127.0.0.1:6161/')
    app.run_server(port=6161)

