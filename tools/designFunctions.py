import cobra

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


class Media:
    '''
    This a media class where different media definitions are stored
    but without the carbon source, which should be specified in the
    'set_medium' function.
    '''
    species = ['escherichia_coli', 'pseudomonas_putida', 'bacillus_subtilis']
    
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
                                               'EX_pi_e' : 1000,
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


    