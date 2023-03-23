import cobra
import pandas as pd
from reframed import to_cobrapy
from reframed.cobra.transcriptomics import gene_to_reaction_expression, GIMME
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px
from numpy import percentile
from pymodulon.enrichment import *
from itertools import compress
import itertools


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


def reaction_sensitivity_to_cutoff(cutoff_range, gene_exp_replicates, rf_media_model, carbon_source, reaction_set='both'):
    '''
    This function computes the set of NBRs for each cutoff value for the given expression
    data and GEM.
    INPUTS:
            cutoff_range : range with the cutoff interval to test
            
            gene_exp_replicates : dict with the replicate ids extracted from expression data
                                  as keys and the expression values for each gene as values
                                  
            rf_media_model : the GEM as a refremed model
            
            carbon_source : dict with carbon source exchange reaction(s) as key and flux(es) as value
            
            reaction_set : str indicating type of reaction set, could be 'NBR', 'BAR' or 'both'
            
    OUTPUT:
            display and saves a plot representing the variation of the number of reactions respect the cutoff
            and returns:
            
            - GIMME_model : a GIMME model containing BAR and NBR set of the most restrictive threshold

            - reaction_set_list : a dict containing the reaction set for each cutoff value and with the pattern:
            
                { <reaction_set> : <reaction_set_list> }
            
            - number_reactions_list : the same as reaction_set_list but containing the lenght of the sets instead
                
    '''
    reaction_set_list = dict()
    number_reactions_list= dict()
    if reaction_set == 'both':
        reaction_set_list['NBR'] = list()
        reaction_set_list['BAR'] = list()
        number_reactions_list['NBR'] = list()
        number_reactions_list['BAR'] = list()
    
    else:
        reaction_set_list[reaction_set] = list()
        number_reactions_list[reaction_set] = list()        
        
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
 
        #See the active reactions in at least one of the replicates according to GIMME:
        try:
            gimme_active_rxns = get_active_rxns(rf_media_model, gene_exp_replicates, cutoff, parsimonious=True)
            
        except:
            print('GIMME gives error with a cutoff value of %s' % str(cutoff))
            
            for reaction_class in number_reactions_list:
                number_reactions_list[reaction_class].append(0)
                
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
        #Generate the context model by deletion those reactions that are inactive according to GIMME and expression data
        GIMME_model = to_cobrapy(rf_media_model)
        GIMME_model.remove_reactions(list(not_active_and_not_in_gimme_model))
        #Compute NBRs & BARs by performing a FVA with fraction_of_optimum set to 80% of wt
        fva_gimme = cobra.flux_analysis.flux_variability_analysis(GIMME_model, fraction_of_optimum=0.8)

        biomass_reactions = set(fva_gimme.loc[(fva_gimme['minimum']!=0) | (fva_gimme['maximum']!=0)].index)
        print('%s reactions are related with biomass component using %s as carbon source(s)' %
              (len(biomass_reactions), ', '.join(list(carbon_source.keys()))) )

        not_biomass_reactions= set(fva_gimme.loc[~((fva_gimme['minimum']!=0) | (fva_gimme['maximum']!=0))].index)
        print('%s reactions are no related with biomass component using %s as carbon source(s)' %
              (len(not_biomass_reactions), ', '.join(list(carbon_source.keys()))) )
        #Populate 'reaction_set_list' by adding NBRs found for this cutoff
        if reaction_set == 'both':         
            reaction_set_list['NBR'].append(not_biomass_reactions)
            reaction_set_list['BAR'].append(biomass_reactions)
            number_reactions_list['NBR'].append(len(not_biomass_reactions))
            number_reactions_list['BAR'].append(len(biomass_reactions))
        
        elif reaction_set != 'both':
            if reaction_set == 'NBR':
                target_set = not_biomass_reactions
            else:
                target_set = biomass_reactions

            reaction_set_list[reaction_set].append(target_set)
            number_reactions_list[reaction_set].append(len(target_set))
        
        else:
            raise Exception("reaction_set can only be 'NBR', 'BAR' or 'both', but %s was given." % reaction_set)
            break
    
    return GIMME_model, reaction_set_list, number_reactions_list

def gene_sensitivity_to_cutoff(reaction_set_list, rf_media_model):
    gene_list = []
    media_model = to_cobrapy(rf_media_model)
    for r_s in reaction_set_list:
        gene_dict = {rxn : [g.id for g in media_model.reactions.get_by_id(rxn).genes]
                              for rxn in r_s}

        gene_list.append(set([g for g_s in list(gene_dict.values()) for g in g_s ]))

    gene_sensitivity_data =  { 'Gene' : [gene.id for gene in media_model.genes],
                               'Presence' : [[g for g_s in gene_list for g in g_s].count(gene.id)/len(reaction_set_list)
                                               for gene in media_model.genes] }
        
    gene_sensitivity_result = pd.DataFrame.from_dict(gene_sensitivity_data)
        
    return gene_sensitivity_result

def get_enrichment_result(query_set, functional_data, rf_media_model, ica_data, functional_class_list, functional_class='Subsystem', input_type='Gene'):
    
    '''
    This function generates a dataframe with enriched functions for a given reaction/gene
    set within a ica data and GEM.
    INPUTS:
            query_set : set containing reactions/genes           
            
            functional_data : dict containing dataframes for all the functional classes whose
                              analysis is permited (iModulon and Subsystem). 
                              Each dataframe has the genes associated to each functional class. 

            rf_media_model : the GEM as a reframed model
            
            ica_data : ica object generated with all available data as described in
                       pymodulon documentation
                       
            functional_class_list : list containing the cathegories of a given functional class
            
            functional_class : string; can be either 'Subsystem' or 'iModulon'
            
            input_type : string; can be either 'Gene' or 'Reaction'
                   
    OUTPUT:       
            returns a dataframe consisting on the enrichment result, which has the columns:
            
                        | <functional_class> | p-Value | recall | target_length |
    '''
    all_genes = ica_data.gene_names

    if input_type == 'Reaction':
        media_model = to_cobrapy(rf_media_model)
        gene_dict = {rxn : [g.id for g in media_model.reactions.get_by_id(rxn).genes]
                     for rxn in query_set}
            
        functional_annotation = get_reaction_functional_annotation(ica_data, gene_dict)
        genes = list(functional_annotation['Gene'].unique())

    pv_list = []
    recall_list = []
    len_list = []
                
    for c in functional_class_list:
        print('Computing %s %s enrichment' % (c, functional_class, ))
        composition_df = functional_data[functional_class]
        if functional_class.startswith('iModulon'):
            function_to_search = functional_class.split('_')[0]
            
        else:
            function_to_search = functional_class

        target_set = composition_df.loc[composition_df[function_to_search]==c, 'Gene'].values.tolist()
        if input_type == 'Reaction':
            enrichment_result = compute_enrichment(genes, target_set, all_genes)
            
        if input_type == 'Gene':
            try:
                enrichment_result = compute_enrichment(query_set, target_set, all_genes)
            except:
                query_set -= query_set-set(all_genes)
                enrichment_result = compute_enrichment(query_set, target_set, all_genes)
                
        pv_list.append(enrichment_result.pvalue)
        recall_list.append(enrichment_result.recall)
        len_list.append(enrichment_result.TP/enrichment_result.target_set_size)
        
    enrichment_data = {functional_class: functional_class_list,
                               'p-Value' : pv_list,
                               'recall' : recall_list,
                               'target_lenght' : len_list}

    enrichment_df = pd.DataFrame.from_dict(enrichment_data).sort_values(by='p-Value')
    
    return enrichment_df

def function_sensitivity_to_cutoff(reaction_set_list, rf_media_model, ica_data, functional_data):
    '''
    This function computes the set of enriched functions for each cutoff value for a given
    reaction set list, ica data, iModulon genes and a GEM.
    INPUTS:
            reaction_set_list : list containing a reaction set for each cutoff value
            
            rf_media_model : the GEM as a reframed model
            
            ica_data : ica object generated with all available data as described in
                       pymodulon documentation            
            
            functional_data : dict containing dataframes for all the functional classes whose
                              analysis is permited (iModulon_<functional class column> and Subsystem). 
                              Each dataframe has the genes associated to each functional class.                          
                   
    OUTPUT:       
            returns a list of dataframes containing the presence of each function across the 
            cutoff values. A element of the list consist on th dateframe for a given functional class           
    '''
    sensitivity_result = []
    media_model = to_cobrapy(rf_media_model)
    enriched_functions_list = []
    for functional_class in functional_data.keys():
        if functional_class.startswith('iModulon'):
            function_to_search, function_category = functional_class.split('_')
            class_list = ica_data.imodulon_names
            im_table = ica_data.imodulon
            
        elif functional_class == 'Subsystem':
            function_to_search = functional_class
            class_list = functional_data[function_to_search].Subsystem.unique().tolist()
         
        else:
            raise Exception( "Keys of functional data need to be either iModulon or Subsystem, but %s was given"
                            % functional_class)
            break
            
        all_genes = ica_data.gene_names
        media_model = to_cobrapy(rf_media_model)
        for g_s in reaction_set_list:
            enrichment_df = get_enrichment_result(g_s, functional_data, rf_media_model, ica_data, class_list, functional_class=functional_class, input_type='Reaction')
            #Consider enriched functional classes those with p-value less than 0.05:
            enriched_fc = enrichment_df.loc[enrichment_df['p-Value'] <= 0.05, functional_class].tolist()
            #For each functional class that our set is enriched in, compute its function
            if function_to_search == 'iModulon':
                enriched_functions = [ f for im in enriched_fc for f in im_table.loc[[im]][function_category].values ]
                enriched_functions_list.append(enriched_functions)
                
            if function_to_search == 'Subsystem':
                enriched_functions_list.append(enriched_fc)

        #For each functional class compute its presence among the reaction sets computed with different cutoff values
        if function_to_search == 'iModulon':
            rxn_sensitivity_data = { 'Function' : [i_f for i_f in im_table[function_category].unique().tolist()],
                                     'Presence' : [ [f for f_s in enriched_functions_list for f in set(f_s) ].count(i_f)/len(reaction_set_list)
                                                for i_f in im_table[function_category].unique().tolist()] }
        
        if function_to_search == 'Subsystem':
            rxn_sensitivity_data = { 'Function' : class_list,
                                     'Presence' : [ [f for f_s in enriched_functions_list for f in set(f_s) ].count(ss)/len(reaction_set_list)
                                                for ss in class_list] } 

        sensitivity_result.append(pd.DataFrame.from_dict(rxn_sensitivity_data))
        
    return sensitivity_result

def get_functional_class_composition(functional_data, imodulon_function, reaction_sets):
    '''                                                     
    This function computes the set of enriched functions for each cutoff value for a given
    reaction set list, ica data, iModulon genes and a GEM.
    INPUTS:
            functional_data : dict containing dataframes for all the functional classes whose
                              analysis is permited (iModulon and Subsystem). Each dataframe has
                              the genes associated to each functional class. 
            
            imodulon_function : dict containing imodulon names as keys and a list of associated 
                                functions as values.
            
            reaction_sets : dict containing reaction class as keys ('NBR' or 'BAR' only for now)
                            and a set of the genes corresponding to the reactions of each one of 
                            the classes.
            
    OUTPUT:       
            returns a dataframe of 3 columns:
            
                    * <Functional class>_functions : names of the functions associated with the class
                    
                    * BAR_Set_percentage : % of each of the functions in the BAR set
                    
                    * NBR_Set_percentage : % of each of the functions in the NBR set
                    
           
            The % of the functions is defined as the ratio between the number of function-specific 
            genes and the number of total genes.
            
    COMMENT:
            maybe a recall value could be interesting here in order to account for functions highly
            present but with few associated genes.            
    '''

    fc_composition_result = {}
    category = None

    for functional_class in functional_data.keys():
        fc_presence_bar = []
        fc_presence_nbr = []
        df = functional_data[functional_class]
        
        if functional_class.startswith('iModulon'):
            function_label, category = functional_class.split('_')
            
        elif functional_class == 'Subsystem':
            function_label = functional_class
            
        else:
            raise Exception("'functional_data' keys can be only 'iModulon' or 'Subsystem' but %s was given" % functional_class)
            break
            
        
        class_list = df[function_label].unique().tolist()
        class_g_dict = { c : df.loc[df[function_label]==c, 'Gene'].values.tolist()
                         for c in class_list }

        if functional_class.startswith('iModulon'):
            function_list = [f for im in class_list for f in imodulon_function[im]]

        elif functional_class == 'Subsystem':
            function_list = class_list

        for c in class_list:
            bar_count=0
            nbr_count=0
            for g in class_g_dict[c]:
                if g in reaction_sets['BAR']:
                    bar_count += 1
                if g in reaction_sets['NBR']:
                    nbr_count += 1

            fc_presence_bar.append((bar_count/len(reaction_sets['BAR']))*100)
            fc_presence_nbr.append((nbr_count/len(reaction_sets['NBR']))*100)

        fc_composition_data = { function_label+'_function' : function_list,
                                'BAR_Set_percentage' : fc_presence_bar,
                                'NBR_Set_percentage' : fc_presence_nbr }

        fc_composition_df = pd.DataFrame.from_dict(fc_composition_data)
        fc_composition_df.set_index(function_label+'_function', inplace=True)
        fc_composition_result[function_label] = fc_composition_df

    return fc_composition_result

def plotly_upset_plot(df):
    # an array of dimensions d x d*2^d possible subsets where d is the number of columns
    subsets = []
    # the sizes of each subset (2^d array)
    subset_sizes = [ ]
    d = len(df.columns)
    for i in range(1, d + 1):
        subsets = subsets + [list(x) for x in list(itertools.combinations(df.columns, i))]

    for s in subsets:
        curr_bool = [1]*len(df)
        for col in df.columns:
            if col in s: curr_bool = [x and y for x, y in zip(curr_bool, list(df.loc[:, col].copy()))]
            else: curr_bool = [x and not y for x, y in zip(curr_bool, list(df.loc[:, col].copy()))]
        subset_sizes.append(sum(curr_bool))


    plot_df = pd.DataFrame({'Intersection': subsets, 'Size':subset_sizes})
    #exclude all combinations not sharing any element 
    plot_df = plot_df.loc[plot_df['Size']!=0]
    plot_df = plot_df.sort_values(by = 'Size', ascending = False)
    #if there are lots of combinations (>25) just show the top25
    if len(plot_df) > 25:
        plot_df = plot_df[:25]

    max_y = max(plot_df['Size'])+0.1*max(plot_df['Size'])

    subsets = list(plot_df['Intersection'])
    scatter_x = []
    scatter_y = []
    for i, s in enumerate(subsets):
        for j in range(d):
            scatter_x.append(i)
            scatter_y.append(-j*max_y/d-0.1*max_y)

    fig = go.Figure()
    #     fig.add_trace(go.Scatter(x=[-1.2,len(subsets)],y= [max_y+0.1*max_y,max_y+0.1*max_y],fill='tozeroy'))
    template =  ['' for x in scatter_x]
    fig_width = 1000
<<<<<<< HEAD
    marker_size = 3*(fig_width/len(scatter_x)) if len(scatter_x) >= 5 else 80
=======
    marker_size = 3*(fig_width/len(scatter_x)) if len(scatter_x) > 7 else 100
>>>>>>> 38e3a675274839b4d5fdf3e6a34e3d58d40ed2eb
        
    fig.add_trace(go.Scatter( x = scatter_x, y = scatter_y,
                                      mode = 'markers', showlegend=False,
                                      marker=dict(size=marker_size,color='#C9C9C9'),
                                      hovertemplate = template
                            )
                 )

    fig.update_layout( xaxis=dict(showgrid=False, zeroline=False),
                               yaxis=dict(showgrid=True, zeroline=False),
                               plot_bgcolor = "#FFFFFF",
                               margin=dict(t=40, l=150)
                     ) 

    for i, s in enumerate(subsets):
        scatter_x_has = []
        scatter_y_has = []
        for j in range(d):
            if df.columns[j] in s:
                scatter_x_has.append(i)
                scatter_y_has.append(-j*max_y/d-0.1*max_y)
                fig.add_trace(go.Scatter( x = scatter_x_has, y = scatter_y_has,
                                                  mode = 'markers', showlegend=False,
                                                  marker=dict(size=marker_size,color='#000000',showscale=False),
                                                  hovertemplate = template
                                                )
                                     )

    fig.update_xaxes(showticklabels=False) # Hide x axis ticks 
    fig.update_yaxes(showticklabels=False) # Hide y axis ticks
    fig.update_traces(hoverinfo=None)

    plot_df['Intersection'] = ['+'.join(x) for x in plot_df['Intersection']]
    template =  [f'<extra><br><b>{lab}</b><br><b>N-Count</b>: {n}</extra>' for  lab, n in zip(plot_df['Intersection'], plot_df['Size'])]
    bar = go.Bar(x = list(range(len(subsets))), y = plot_df['Size'], marker = dict(color='#000000'),  text = plot_df['Size'], hovertemplate = template, textposition='outside', hoverinfo='none')
    fig.add_trace(bar)

    template =  ['' for x in range(d)]
    max_string_len = max([len(x) for x in df.columns])
    fig_lab = go.Scatter(x = [-0.01*max_string_len]*d, y = scatter_y, text = df.columns, mode = 'text', textposition='middle left',showlegend=False, hovertemplate = template)
    fig_lab = go.Scatter(x = [-0.01*max_string_len]*d, y = scatter_y, text = df.columns, mode = 'text', textposition='middle left',showlegend=False, hovertemplate = template)
    fig.add_trace(fig_lab)
    fig.update_layout(title = '<b>Intersections<b>',height=800, width=fig_width, yaxis_range=[-max_y-0.1*max_y-1,max_y+0.1*max_y], xaxis_range = [-0.13*max_string_len, len(subsets)], showlegend = False, title_x=0.5)

    return fig 
