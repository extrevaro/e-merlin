import cobra
import pandas as pd
from reframed import to_cobrapy
from reframed.cobra.transcriptomics import gene_to_reaction_expression, GIMME
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from numpy import percentile
from pymodulon.enrichment import *

supported_species = ['escherichia_coli', 'pseudomonas_putida', 'bacillus_subtilis']

def set_medium(model, media_definition, carbon_source):
    model_with_media = model.copy()
    
    try:
        for source, uptake in carbon_source.items():
            media_definition[source] = uptake
            
        model_with_media.medium = media_definition
            
    except Exception as e:
        print(e)
        if type(carbon_source) != dict or type(media_definition) != dict:
            print("Both 'media_definition'  and 'carbon_source' must be a dict with pattern {<exchange_reaction> : <uptake_rate>}")
    
    return model_with_media


def reformat_bounds(model):
    working_model= model.copy()
    for reaction in working_model.reactions:
        reaction.bounds = tuple( b if abs(b)<1000 else (b/abs(b))*1000 for b in reaction.bounds )
    return working_model


def get_active_rxns(rf_media_model, gene_exp_replicates, cutoff, parsimonious=False):
    replicate_list = []
    for replicate in gene_exp_replicates:
        gene_exp = gene_exp_replicates[replicate]
        gimmes_sol = GIMME(rf_media_model, gene_exp, cutoff=cutoff, parsimonious=parsimonious)
        gimmes_result = gimmes_sol.to_dataframe()
        gimmes_result = gimmes_result.rename(columns={'value': replicate+'_value'})
        replicate_list.append(gimmes_result)

        active_rxns = pd.concat(replicate_list, axis=1)
        #delete rows whose all values are 0
        active_rxns = active_rxns[~(active_rxns == 0).all(axis=1)]                                       
        active_rxns = active_rxns.assign(Mean=active_rxns.mean(1), Stddev=active_rxns.std(1))

    return active_rxns

#To apply this function you need to have a fully populated ICA structure, 
#including the iModulon and Annotation tables. Moreover you also need a 
#dictionary with this pattern : {<reaction> : <GPR(s)>}
def get_reaction_functional_annotation(ica_data, rxn_gene_dict):
    annot_table = ica_data.gene_table
    rxn_iM_cog_data = {'Reaction': [],
                       'Gene' : [],
                       'iModulon' : [],
                       'COG' : []}

    for r in rxn_gene_dict.keys():
        for g in rxn_gene_dict[r]:
            try:
                rxn_iM_cog_data['iModulon'].append(ica_data.imodulons_with(g))
                rxn_iM_cog_data['Reaction'].append(r)
                rxn_iM_cog_data['Gene'].append(g)
                rxn_iM_cog_data['COG'].append(annot_table.loc[g].cog)

            except Exception as e:
                print(e)
                rxn_iM_cog_data.pop(r, None)
                
    rxn_iM_cog_df = pd.DataFrame.from_dict(rxn_iM_cog_data)
    
    return rxn_iM_cog_df


def get_network_data(cobra_model, active_rxns, interest_condition):
    sustrates = set([s.id for rxn_id in active_rxns.index for s in cobra_model.reactions.get_by_id(rxn_id).reactants])
    s_l = []
    p_l = []
    r_l = []
    f_l = []

    for s in sustrates:
        for r in cobra_model.metabolites.get_by_id(s).reactions:
            if r.id in active_rxns.index:
                for p in r.products:
                    s_l.append(s)
                    p_l.append(p.id)
                    r_l.append(r.id)
                    #it is correct to use the mean or it is better to use a sample value
                    f_l.append(active_rxns.loc[r.id, interest_condition+'_value'])

    all_metabolites = s_l+p_l
    network_df = pd.DataFrame.from_dict({'sustrate':s_l, 
                                         'product':p_l, 
                                         'reaction':r_l,
                                         'flux':f_l})
    
    #Delete all reactions whose flux is 0 for this condition
    network_df = network_df.loc[~(network_df['flux'] == 0)]
    
    return network_df

class Media:
    '''
    This a media class where different media definitions are stored
    but without the carbon source, which should be specified in the
    'set_medium' function.
    '''
    species = supported_species
    
    media_dict = { 'pseudomonas_putida_m9' : { 'EX_ca2_e' : 1000,
                                               'EX_cl_e' : 1000,
                                               'EX_co2_e' : 1000,
                                               'EX_cobalt2_e' : 1000,
                                               'EX_cu2_e' : 1000,
                                               'EX_fe2_e' : 1000,
                                               'EX_fe3_e' : 1000,
                                               'EX_h_e' : 1000,
                                               'EX_h2o_e' : 1000,
                                               'EX_k_e' : 1000,
                                               'EX_mg2_e' : 1000,
                                               'EX_mn2_e' : 1000,
                                               'EX_mobd_e' : 1000,
                                               'EX_na1_e' : 1000,
                                               'EX_tungs_e' : 1000,
                                               'EX_zn2_e' : 1000,
                                               'EX_ni2_e' : 1000,
                                               'EX_sel_e' : 1000,
                                               'EX_so4_e' : 1000,
                                               'EX_nh4_e' : 1000,
                                               'EX_pi_e' : 1.024, #calculated by Blas in minimal media
                                               'EX_cbl1_e' : .01,
                                               'EX_o2_e' : 20 },
                  
                   'escherichia_coli_m9' : {   'EX_ca2_e' : 1000,
                                               'EX_cl_e' : 1000,
                                               'EX_co2_e' : 1000,
                                               'EX_cobalt2_e' : 1000,
                                               'EX_cu2_e' : 1000,
                                               'EX_fe2_e' : 1000,
                                               'EX_fe3_e' : 1000,
                                               'EX_h_e' : 1000,
                                               'EX_h2o_e' : 1000,
                                               'EX_k_e' : 1000,
                                               'EX_mg2_e' : 1000,
                                               'EX_mn2_e' : 1000,
                                               'EX_mobd_e' : 1000,
                                               'EX_na1_e' : 1000,
                                               'EX_tungs_e' : 1000,
                                               'EX_zn2_e' : 1000,
                                               'EX_ni2_e' : 1000,
                                               'EX_sel_e' : 1000,
                                               'EX_slnt_e' : 1000,
                                               'EX_so4_e' : 1000,
                                               'EX_nh4_e' : 1000,
                                               'EX_pi_e' : 1000,
                                               'EX_cbl1_e' : .01,
                                               'EX_o2_e' : 20}
                 }

    def __init__(self, s=str(), n=str()):
        self.species = s
        self.name = n
        self.definition = self.media_dict['_'.join([s,n])]

    def get_predefined_media_definition(self):
        return self.definition
    

class OD_coefficient:
    '''
    This a coefficient class where different species-level coefficients 
    are stored to enable programatic conversion from OD to biomass flux
    that is a required input for the function GIMMES.
    '''
    species = supported_species    

    coeff_dict = { 'pseudomonas_putida_m9' : 0.6,
                   'escherichia_coli_m9' : 0.54,
                   'bacillus_subtilis_m9' : 0.48
                 }
    
    def __init__(self, s=str(), m=str()):
        self.species = s
        self.media = m
        self.definition = self.coeff_dict['_'.join([s,m])]

    def get_predefined_media_definition(self):
        return self.definition

def plot_enrichments_results(enrichment_result, functional_class):
    p_value_fig = px.bar(enrichment_result, x=functional_class, y=['p-Value_BAR', 'p-Value_NBR'],
                 barmode='group',
                 height=250,
                 opacity = [l/max(enrichment_result.target_lenght_BAR.tolist()) 
                             for l in enrichment_result.target_lenght_BAR.tolist()] +
                           [l/max(enrichment_result.target_lenght_NBR.tolist())
                             for l in enrichment_result.target_lenght_NBR.tolist()])

    recall_fig = px.bar(enrichment_result, x=functional_class, y=['recall_BAR', 'recall_NBR'],
                 barmode='group',  color_discrete_sequence=px.colors.qualitative.D3,
                 height=250,
                 opacity = [l/max(enrichment_result.target_lenght_BAR.tolist()) 
                             for l in enrichment_result.target_lenght_BAR.tolist()] +
                           [l/max(enrichment_result.target_lenght_NBR.tolist())
                             for l in enrichment_result.target_lenght_NBR.tolist()])


    # For as many traces that exist per Express figure, get the traces from each plot and store them in an array.
    # This is essentially breaking down the Express fig into it's traces
    figure1_traces = []
    figure2_traces = []
    for trace in range(len(p_value_fig["data"])):
        figure1_traces.append(p_value_fig["data"][trace])

    for trace in range(len(recall_fig["data"])):
        figure2_traces.append(recall_fig["data"][trace])

    #Create a 2x1 subplot
    this_figure = make_subplots(rows=2, cols=1, vertical_spacing=0.02 )

    # Get the Express fig broken down as traces and add the traces to the proper plot within in the subplot
    for traces in figure1_traces:
        this_figure.append_trace(traces, row=1, col=1)

    for traces in figure2_traces:
        this_figure.append_trace(traces, row=2, col=1)

    this_figure.update_xaxes(visible=False, showticklabels=False)
    this_figure['layout']['yaxis'].update(title='p Value')
    this_figure['layout']['yaxis2'].update(title=functional_class+' recall')
    this_figure['layout']['xaxis2'].update(title=functional_class, visible=True, showticklabels=True, tickangle=45)

    this_figure.update_layout(
        title=dict(
            text='<b>'+functional_class+' enrichment analysis in BARs & NBRs</b>',
            x=0.5,
            y=0.95,
            font=dict(
                family="Arial",
                size=20,
                color='#000000'
            )
        )
    )

    return this_figure


def nbr_sensitivity_to_cutoff(cutoff_range, gene_exp_replicates, rf_media_model, carbon_source):
    '''
    This function computes the set of NBRs for each cutoff value for the given expression
    data and GEM.
    INPUTS:
            cutoff_range : range with the cutoff interval to test
            
            gene_exp_replicates : dict with the replicate ids extracted from expression data
                                  as keys and the expression values for each gene as values
                                  
            rf_media_model : the GEM as a refremed model
            
            carbon_source : dict with carbon source exchange reaction as key and flux as value
            
    OUTPUT:
            display a plot representing the variation of the number of NBRs respect the cutoff
            
            returns a list containing the NBRs for each cutoff value            
    '''
    nbr_reaction_set_list = []
    number_of_nbr_list=[]
    for cutoff in cutoff_range:
        print('Applying method for cutoff value of %s' % str(cutoff))
        expression_set = set()
        #See the active reactions in at least one of the replicates according to expression data:
        for replicate in gene_exp_replicates:
            rxn_exp = gene_to_reaction_expression(rf_media_model, gene_exp_replicates[replicate], or_func=max)
            threshold = percentile([abs(exp) for exp in list(rxn_exp.values())], cutoff)
            print(threshold)
            rxn_exp = {r_id for r_id, exp in rxn_exp.items() if exp > threshold}
            expression_set |= rxn_exp

        expression_set

        #See the active reactions in at least one of the replicates according to GIMME:
        try:
            gimme_active_rxns = get_active_rxns(rf_media_model, gene_exp_replicates, cutoff, parsimonious=True)
            
        except:
            print('GIMME gives error with a cutoff value of %s' % str(cutoff))
            number_of_nbr_list.append(0)
            continue
            
        active_but_not_in_gimme_model  = set([rxn for rxn in expression_set if rxn not in gimme_active_rxns.index])
        print('%s reactions are being expressed according RNAseq but not active according GIMME' % 
              len(active_but_not_in_gimme_model))

        active_and_in_gimme_model = set([rxn for rxn in expression_set if rxn in gimme_active_rxns.index])
        print('%s reactions are being expressed according RNAseq and GIMME' % 
              len(active_and_in_gimme_model))

        not_rxn_exp = set(rf_media_model.reactions)-expression_set

        not_active_but_in_gimme = set([rxn for rxn in not_rxn_exp if rxn in gimme_active_rxns.index])
        print('%s reactions are not being expressed according RNAseq but active according GIMME' % 
              len(not_active_but_in_gimme))

        not_active_and_not_in_gimme_model = set([rxn for rxn in not_rxn_exp if rxn not in gimme_active_rxns.index])
        print('%s reactions are not being expressed according RNAseq and GIMME' % 
              len(not_active_and_not_in_gimme_model))
        #Generate the context model by deletio¡ng those reactions that are inactive according to GIMME and expression data
        GIMME_model = to_cobrapy(rf_media_model)
        GIMME_model.remove_reactions(list(not_active_and_not_in_gimme_model))
        GIMME_model.optimize()
        #Compute NBRs & BARs by performin a FVA with fraction_of_optimum set to 80% of wt
        fva_gimme = cobra.flux_analysis.flux_variability_analysis(GIMME_model, fraction_of_optimum=0.8)

        biomass_reactions = set(fva_gimme.loc[(fva_gimme['minimum']!=0) | (fva_gimme['maximum']!=0)].index)
        print('%s reactions are related with biomass component using %s as carbon source' %
              (len(biomass_reactions), list(carbon_source.keys())[0]))

        not_biomass_reactions= set(fva_gimme.loc[~((fva_gimme['minimum']!=0) | (fva_gimme['maximum']!=0))].index)
        print('%s reactions are no related with biomass component using %s as carbon source' %
              (len(not_biomass_reactions), list(carbon_source.keys())[0]))
        #Populate 'nbr_reaction_set_list' by adding NBRs foun for this cutoff
        nbr_reaction_set_list.append(not_biomass_reactions)
        number_of_nbr_list.append(len(not_biomass_reactions))
        
    #Plot the nº of NBRs respect the cutoff value
    fig = px.scatter(x=[c for c in cutoff_range], y=number_of_nbr_list)
    fig.update_layout(title=dict(text='<b>NBR Sensitivity to cutoff value<b>',
                                 x=0.5),
                      xaxis=dict(title='Cutoff Value'),
                      yaxis=dict(title='Number of NBRs in GIMME model')
                     )
    display(fig)
    
    return nbr_reaction_set_list


def function_sensitivity_to_cutoff(reaction_set_list, gene_presence, ica_data, rf_media_model):
    '''
    This function computes the set of enriched functions for each cutoff value for a given
    reaction set list, ica data, iModulon genes and a GEM.
    INPUTS:
            reaction_set_list : list containing the NBRs for each cutoff value
            
            gene_presence : dataframe containing the genes associated to each iModulon
            
            ica_data : ica object generated with all available data as described in
                       pymodulon documentation
                                  
            rf_media_model : the GEM as a refremed model
            
    OUTPUT:       
            returns a dataframe containing the presence of each function across the cutoff values            
    '''
    all_genes = ica_data.gene_names
    im_list = ica_data.imodulon_names
    im_table = ica_data.imodulon
    media_model = to_cobrapy(rf_media_model)
    enriched_functions_list = []
    for g_s in reaction_set_list:
        nbr_gene_dict = {rxn : [g.id for g in media_model.reactions.get_by_id(rxn).genes]
                      for rxn in g_s}
        nbr_functional_annotation = get_reaction_functional_annotation(ica_data, nbr_gene_dict)
        nbr_functional_annotation['Reaction_type'] = ['NBR']*len(nbr_functional_annotation)
        nbr_genes = list(nbr_functional_annotation['Gene'].unique())

        nbr_pv_list = []
        nbr_recall_list = []

        for im in im_list:
            target_set = gene_presence.loc[gene_presence['iModulon']==im, 'Gene'].values.tolist()
            nbr_enrichment_result = compute_enrichment(nbr_genes, target_set, all_genes)
            nbr_pv_list.append(nbr_enrichment_result.pvalue)
            nbr_recall_list.append(nbr_enrichment_result.recall)

        enrichment_data = {'iModulon': im_list,
                           'p-Value_NBR' : nbr_pv_list,
                           'recall_NBR' : nbr_recall_list}

        enrichment_df = pd.DataFrame.from_dict(enrichment_data).sort_values(by='p-Value_NBR')
        #Consider enriched iModulons those with p-value less than 0.05:
        nbr_enriched_im = enrichment_df.loc[enrichment_df['p-Value_NBR'] <= 0.05, 'iModulon'].tolist()
        #For each iModulon that our set is enriched in, compute its function
        nbr_enriched_functions = [ f for im in nbr_enriched_im for f in im_table.loc[[im]].function.values ]
        enriched_functions_list.append(nbr_enriched_functions)

    #For each iModulon function compute its presence among the NBR sets computed with different cutoff values
    sensitivity_data = { 'Function' : [i_f for i_f in im_table.function.unique().tolist()],
                         'Presence' : [ [f for f_s in enriched_functions_list for f in set(f_s) ].count(i_f)/len(reaction_set_list)
                                        for i_f in im_table.function.unique().tolist()] }

    sensitivity_result = pd.DataFrame.from_dict(sensitivity_data)
    return sensitivity_result


        