#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd

import dash
from dash.dependencies import Input, Output
from dash import dcc, html, dash_table

from zipfile import ZipFile
from glob import glob
import urllib


# In[ ]:


batches = glob("assets/*.zip")

omim_map = pd.read_csv('assets/OMIM_map.txt', sep = '\t')

MGI_IDs = pd.read_csv('assets/MGI_IDs.csv', dtype = str)
MGI_IDs = dict(zip(MGI_IDs['Symbol'], MGI_IDs['MGI_ID']))

hpo = pd.read_csv('assets/HPO_dict.txt', sep = '\t')
hpo = hpo[hpo['Associated Genes'].str.len() > 1]

refSeqSummaries = pd.read_csv('assets/RefSeqSummaries.tsv', sep = '\t', dtype = str, low_memory = False)

default_style = {'width' : '57%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}
default_style = {'width' : '100%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}

hide = {'display' : 'None'}
full_display = {'width' : '100%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'textAlign' : 'center'}
info_style = {'width' : '100%', 'height' : '100%', 'margin' : '10px', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'left'}
break_line = html.Hr(style={'height' : '4px', 'width' : '60%', 'color' : '#111111','display' : 'inline-block', 'marginLeft':'auto', 'marginRight':'auto'})


# In[ ]:


def generateHeatTable(gene, Chromosome, Position):
    df = pd.read_csv('assets/CleanerVar.tsv', sep = '\t', low_memory = False) 
    df = df[df['Gene'] == gene]
    if df.shape[0] == 0:
        return html.P('Could not find pathogenic variants for this gene in ClinVar🤷 🤷‍♂️')
    del df['Gene']
    del df['Chromosome']
    df.fillna('', inplace = True)
    df['Reference'] = df['Reference'].astype(str).str.replace('na','')
    df['Variant'] = df['Variant'].astype(str).str.replace('na','')
    df['Variation'] = df.apply(lambda x : x['Reference'] + '>' + x['Variant'], axis = 1)
    del df['Reference']
    del df['Variant']
    df['Variation'] = df.apply(lambda x : x['Mutation Type'] if x['Variation'] == '>' else x['Variation'], axis = 1)
    del df['Mutation Type']
    df['VariationID'] = df['VariationID'].apply(lambda x : str(int(x)))
    df['Link'] = 'https://www.ncbi.nlm.nih.gov/clinvar/variation/' + df['VariationID']
    del df['VariationID']
    df = df.append({
        'Significance' : '👩‍⚕️👨‍⚕️',
        'Phenotypes' : '👩‍⚕️👨‍⚕️',
        'Position' : int(Position),
        'Variation' : 'Your variant',
        'Link' : ''
    }, ignore_index=True)
    df.sort_values('Position', inplace = True)
    if df['Position'].max() - int(Position) >= int(Position) - df['Position'].min():
        maxSpace = df['Position'].max() - int(Position)
    else:
        maxSpace = int(Position) - df['Position'].min()
    df['gradient'] = (1 - abs(int(Position) - df['Position']) / maxSpace) * 255 ** 0.7
    df['gradient'] = df['gradient'].astype(int)
    df.at[df.index[df['Position'] == Position].tolist(), 'gradient'] = 255
    df['gradient'] = df['gradient'].apply(lambda x : 'rgb(255,' + str(255 - x) + ',' + str(255 - x)+')')
    df['Position'] = df['Position'].apply(lambda x : str(Chromosome) + ':' + str(x))
    cols = ('Position', 'Phenotypes', 'Variation', 'Link')
    df['Variation'] = df['Variation'].apply(lambda x : x[:25] + '(...)' if len(x) > 25 else x)
    return html.Table(
        [html.Tr([html.Th(col) for col in cols], style = {'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'})] +
        [html.Tr([
            html.Td(df.iloc[i][col], style = {'backgroundColor' : df.iloc[i]['gradient'], 'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'}) if col != 'Link' else html.Td(html.A(href=df.iloc[i]['Link'], children=df.iloc[i][col], target='_blank'), style = {'backgroundColor' : df.iloc[i]['gradient'], 'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'}) for col in cols
        ]) for i in range(len(df))],
        style = {'textAlign' : 'center', 'border': '1px solid black', 'border-collapse' : 'collapse'}
    )

def sortBatches(batches):
    options = [{'label': v[7:-4], 'value': v} for v in sorted(batches)]
    return options


# In[ ]:


external_stylesheets = ['assets/BAMdelbee_main.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "BAMdelbee | Ohad Birk's Lab"
app.config['suppress_callback_exceptions'] = False


# In[ ]:


app.layout = html.Div([
    html.Div([
        html.Img(id = 'header_image',
            src = 'assets/logo_BAMdelbee.png',
            style = {
                'height' : '16%',
                'width' : '16%',
                'float' : 'left',
                'marginLeft' : 4
            }
        ),
        html.Span(''),
        html.Img(id = 'header_image2',
            src = 'assets/logo_birklab.png',
            style = {
                'height' : '16%',
                'width' : '16%',
                'float' : 'right',
                'margin' : 4
            }
        ),
    ]),
    html.P(html.Br()),
    html.H2(
        '',
        id = 'selected_batch',
        style = {'width' : '100%', 'height' : '100%', 'font-family' : 'gisha', 'margin' : 'auto', 'textAlign' : 'center'}
    ),
    html.P(
        'Select sample/s:',
        id = 'selected_samples',
        style = {'display' : 'none'}
    ),
    html.Div([
        dcc.Dropdown(
            id = 'BatchesDropdown',
            options = sortBatches(batches),
            multi = False,
            placeholder = 'Select a batch',
            style = {'width' : '70%', 'height' : '100%', 'font-family' : 'gisha', 'margin' : 'auto', 'textAlign' : 'center'}
        ),
        html.Br(id = 'spacer1'),
        dcc.RadioItems(
            id = 'referenceRadioButtons',
            options=[
                {'label': 'hg38', 'value': 'hg38'},
                {'label': 'hg19', 'value': 'hg19'},
            ],
            value='hg38',
            labelStyle={'display': 'inline-block'},
            style = {'width' : '100%', 'height' : '100%', 'font-family' : 'gisha', 'margin' : 'auto', 'textAlign' : 'center'}
        ),
        ],
        style = {'width' : '100%', 'height' : '100%', 'font-family' : 'gisha', 'margin' : 'auto', 'textAlign' : 'center'}
    ),
    
    
    html.P(html.Br(), id = 'temp_break'),
    
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
     Output('selected_batch', 'children'),
     Output('checkListContainer', 'children'),
     Output('temp_break', 'style'),
     Output('selected_samples', 'style'),
     Output('spacer1', 'style')],
    [Input('BatchesDropdown', 'value')]
)
def chooseBatch(batch):
    if batch != None:
        zBatch = ZipFile(batch)
        samples = zBatch.read('samples.txt').decode('utf-8', 'ignore').split('\n')
        samples = sorted(samples)
        samples = [{'label' : str(sample), 'value' : str(sample)} for sample in samples]
        checkList = dcc.Checklist(id = 'samplesCheckList', options = samples, labelStyle={'display': 'block'}) 
        btn = html.Button('Analyze', id = 'Analysis_starting_button', disabled = True, n_clicks = 0, style = {'text-transform' : 'none', 'font-size' : '30px'})
        checklist_style = {'width' : '100%', 'height' : '100%', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}
        return [hide, hide, checklist_style, batch[7:-4], [checkList, btn], {'display' : 'none'}, default_style, {'display' : 'none'}]
    else:
        dropdown_style = {'width' : '70%', 'height' : '100%', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}
        ref_style = {'width' : '100%', 'height' : '100%', 'font-family' : 'gisha', 'marginLeft' : 'auto', 'marginRight' : 'auto', 'textAlign' : 'center'}
        return [dropdown_style, ref_style, hide, '', None, default_style, {'display' : 'none'}, {'display' : 'inline-block'}]

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
     Output('loading_results_div', 'children'),
     Output('selected_samples', 'children')],
    [Input('Analysis_starting_button', 'n_clicks'),
     Input('BatchesDropdown', 'value'),
     Input('samplesCheckList', 'value'),
     Input('referenceRadioButtons', 'value')]
)
def update_figure(clicks, batch, affected, referenceGenome):
    if clicks in (0, None):
        return [default_style, hide, [], None, 'Select sample\s:']
    else:
        if referenceGenome == 'hg38':
            geneMap = pd.read_csv('assets/NCBI_Genes_hg38.map', low_memory = False)
        else:
            geneMap = pd.read_csv('assets/NCBI_Genes_hg19.map', low_memory = False)
        geneMap['Chr'] = geneMap['Chr'].astype(str)
        
        zBatch = ZipFile(batch)
        deletions = []
        sizes = []
        genes = []
        isCoding = []
        exons = pd.read_csv('assets/' + referenceGenome + '_exons.csv.gz', low_memory = False)
        exons['Chromosome'] = exons['Chromosome'].astype(str)
        
        for chromosome in [str(c) for c in range(1,23)] + ['X', 'Y']:
            try:
                df = pd.read_csv(zBatch.open('chr' + chromosome + '.coverage'), sep = '\t', low_memory = False)
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
                        coding = exons[exons['Gene'].isin(genes_in_deletion)]
                        coding = coding[coding['Chromosome'] == chromosome]
                        coding = coding[(coding['Start'] >= start) & (coding['Start'] <= end)]
                        if coding.empty:
                            isCoding.append('')
                        else:
                            isCoding.append('Yes')
        df = pd.DataFrame()
        df['Coordinates'] = deletions
        df['Chromosome'] = df['Coordinates'].apply(lambda x : int(x.split(':')[0].upper().replace('X', '23').replace('Y', '24').replace('M', '25').replace('MT', '25')))
        df['~Size (bp)'] = sizes
        df['Genes'] = genes
        df['Start'] = df['Coordinates'].apply(lambda x : int(x.split(':')[1].split('-')[0]))
        df.sort_values(['Chromosome', 'Start'], inplace = True)
        df.reset_index(drop = True, inplace = True)
        df['End'] = df['Coordinates'].apply(lambda x : int(x.split(':')[1].split('-')[1]))
        df['Coding'] = isCoding
        del df['Chromosome']
        del df['Start']
        del df['End']
        df = df[df['~Size (bp)'] > 50]
        df = df[df['Genes'] != '']
        
        if df.empty:
            return [hide, full_display, [html.H2('🤷 No deletion were found 🤷‍♂️')], None, ', '.join(affected)]
        
        ddt = dash_table.DataTable(
            id = 'table',
            data = df.to_dict('records'),
            columns = [{'id': c, 'name': c} for c in df.columns],
            page_action = 'naive',
            page_size = 25,
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
        csv_report = df.fillna('').to_csv(na_rep = '', index = False).replace(',nan,', ',,').replace('[','').replace(']','').replace("'",'')
        download_href = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_report)
        download_image = html.A('Download table ⬇️', href = download_href, download = 'BAMdelbee_' + '_'.join(affected) + '.csv', title = 'Download table', style = {'text-decoration' : 'none', 'fontSize' : '150%', 'float' : 'right', 'marginTop' : 1})
        info = html.Div('Click on a deletion for more information', id = 'info')
        return [hide, full_display, [ddt, download_image, info], None, ', '.join(affected)]

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
        
        #IGV & UCSC % gnomAD
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
                hpos = html.Details(
                    [
                    html.Summary('HPO associated phenotypes', style = {'fontStyle' : 'oblique'}),
                    html.Div(html.Div(phenotypes))
                    ],
                )
                elements.append(hpos)
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
            clinVarBtn = html.Button('ClinVar HeatTable', id = 'HeatTableBtn', n_clicks = 0, style = {'text-transform' : 'none'})
            heatTableDiv = html.Div(id = 'heatTableDiv', style = {'display' : 'inline-block', 'margin':'auto'})
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
        return [[None], None, {'display' : 'inline-block', 'text-transform' : 'none'}]
    else:
        coordinates = data[selected_row_index[0]]['Coordinates']
        gene = data[selected_row_index[0]]['Genes']
        if ', ' in gene:
            gene = gene.split(', ')[0]
        return [generateHeatTable(gene, coordinates.split(':')[0], coordinates.split(':')[1].split('-')[0]), None, {'display' : 'None'}]


# In[ ]:


if __name__ == '__main__':
    import webbrowser
    webbrowser.open('http://127.0.0.1:6261/')
    app.run_server(port=6261)

