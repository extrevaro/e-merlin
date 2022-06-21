import cobra
import pandas as pd
from reframed.cobra.transcriptomics import GIMME

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
        gimmes_sol = GIMME(rf_media_model, gene_exp, parsimonious=parsimonious)
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

    