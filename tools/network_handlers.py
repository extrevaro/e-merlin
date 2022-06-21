import networkx as nx
from math import floor, sqrt
import plotly.graph_objects as go
from .designFunctions import get_active_rxns, get_network_data

def get_subnetwork_centers(subnetworks):
    pos = []
    d = 0
    for sg in subnetworks:
        pos.append((d,0))
        d += sqrt(len(sg.nodes()))
            
    return pos

def metabolites_to_reactions(metabolite_pair, network_df):
    rxn_names = network_df.loc[ (network_df['sustrate'] == metabolite_pair[0]) &
                                (network_df['product'] == metabolite_pair[1]),
                                'reaction' ].unique().tolist()
    return rxn_names


def draw_subnetworks(network_list, flux_threshold=0.05):
    data = []
    #determine the center for each subgraph:
    centers = get_subnetwork_centers(network_list)
    max_adjacency = max([len(a[1]) for SG in network_list for a in SG.adjacency()])
    #Set layout for get tho coordenates of each node:
    for SG in network_list:
        pos = nx.spring_layout(SG,
                               center=centers[network_list.index(SG)], scale=len(SG.nodes())/max([sqrt(len(sg.nodes())) for sg in network_list]))
        #Set edge position and text
        edge_x = []
        edge_y = []
        edge_text = []
        edge_flux_list = []
        for edge in SG.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)
            edge_text.append('_'.join(edge))
            edge_flux_list.append(SG[edge[0]][edge[1]]['flux'])

        n_of_coordinates = len(edge_x)/len(edge_text)
        low_flux_edge_trace = go.Scatter(
                                x= [x for x in edge_x if edge_flux_list[floor(edge_x.index(x)/n_of_coordinates)]<=flux_threshold], 
                                y= [y for y in edge_y if edge_flux_list[floor(edge_y.index(y)/n_of_coordinates)]<=flux_threshold],
                                mode='lines',
                                text = [t for t in edge_text if edge_flux_list[edge_text.index(t)]<flux_threshold],
                                line=dict(width=0.5, color= '#808888'),
                                hoverinfo='text')
        data.append(low_flux_edge_trace)


        high_flux_edge_trace = go.Scatter(
                                x= [x for x in edge_x if abs(edge_flux_list[floor(edge_x.index(x)/n_of_coordinates)])>flux_threshold], 
                                y= [y for y in edge_y if abs(edge_flux_list[floor(edge_y.index(y)/n_of_coordinates)])>flux_threshold],
                                mode='lines',
                                text = [t for t in edge_text if edge_flux_list[edge_text.index(t)]>flux_threshold],
                                line=dict(width=2, color= '#990099'),
                                hoverinfo='text')
        data.append(high_flux_edge_trace)

        #Set node position and text
        node_x = []
        node_y = []
        node_text = []
        node_adjacencies = []
        for node in SG.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            adjacenty = [len(a[1]) for a in SG.adjacency() if a[0]== node][0]
            node_adjacencies.append(adjacenty)
            node_text.append(' \n '.join([node, 'Connections: '+str(adjacenty)]))

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            text = node_text,
            hoverinfo='text',
            marker=dict(
                showscale=True,
                #to assure each graph has the same scale:
                cmax=max_adjacency,
                cmin=0,
                # colorscale options
                #'Greys' | 'YlGnBu' | 'Greens' | 'YlOrRd' | 'Bluered' | 'RdBu' |
                #'Reds' | 'Blues' | 'Picnic' | 'Rainbow' | 'Portland' | 'Jet' |
                #'Hot' | 'Blackbody' | 'Earth' | 'Electric' | 'Viridis' |
                colorscale='RdBu',
                reversescale=True,
                color=[],
                size=6,
                opacity=0.8,
                colorbar=dict(
                    thickness=15,
                    title='Node Connections',
                    xanchor='left',
                    titleside='right'
                ),
                line_width=2))

        #Set custom visualization
        node_trace.marker.color = node_adjacencies
        data.append(node_trace)

    fig = go.Figure(data=data,
                 layout=go.Layout(
                    title = 'Subnetworks originated from GIMMES processing',
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                    )

    fig.show()
    

def extract_subnetworks(rf_media_model, gene_exp_replicates, node_threshold = 20, cutoff_range = range(1,34), parsimonious=False):
    # The function has as input:
    #
    #    Required:
    #        'rf_media_model', and 'gene_exp_replicates'
    #    Optional:
    #        'node_threshold' (default -> 20) 'cutoff_range' (default -> range(1,34))
    #     
    #
    #   Needs to use pandas ,reframed and Networkx
    # 
    # It returns:
    #    'subnetwork_dict' : a dict containing all the different subnetworks found with its total counts 
    #    
    import pandas as pd
    from reframed import to_cobrapy
    from reframed.cobra.transcriptomics import GIMMES

    subnetwork_dict = {}
    for i in cutoff_range:
        #Get the active reaction predicted by GIMMES for each replicate,
        #Only reactions whose flux is 0 in all the replicates are excluded
        active_rxns = get_active_rxns(rf_media_model, gene_exp_replicates, i, parsimonious)

        #generate the network(s) dataframe
        cobra_model = to_cobrapy(rf_media_model)
        for replicate in gene_exp_replicates:
            #Compute network data for each condition. Here all reactions
            #whose value is 0 for the given condition will be excluded
            network_df = get_network_data(cobra_model, active_rxns, replicate)
            
            #Find the different subnetworks and plot them:
            G=nx.from_pandas_edgelist(network_df, 'sustrate', 'product', edge_attr=['reaction', 'flux'], create_using=nx.DiGraph)
            UG=G.to_undirected()
            S = [UG.subgraph(c).copy().to_directed() for c in nx.connected_components(UG)]
            subgraphs = [sg for sg in S if len(sg.nodes())<=node_threshold]
            print('Using fluxes from %s, %s network(s) were found for cutoff %s' % (replicate, str(len(S)), str(i)))
            if len(subgraphs) >= 1:
                print('Their numbers of nodes are:')
                print('%s | '*len(S) % tuple(str(len(sg.nodes())) for sg in S) )
                for sg in subgraphs:
                    graph_id = []
                    print(sg.edges())
                    for m_p in sg.edges():
                        rxns = metabolites_to_reactions(m_p, network_df)
                        if '-'.join(rxns) not in graph_id:
                            graph_id.append('-'.join(rxns))
                            
                        print(graph_id)
                    
                    graph_id = '-'.join(graph_id)
                    
                    if graph_id not in subnetwork_dict.keys():
                        new_dict = { graph_id : [None, {r : [i] if r == replicate else [] for r in gene_exp_replicates}]}
                    else:
                        new_dict[graph_id][1].update({replicate : subnetwork_dict[graph_id][1][replicate]+[i]})

                print(new_dict)
                subnetwork_dict = new_dict
                
    unique_cutoff_values = set()
    for sg in subnetwork_dict.keys():
        for c_l in subnetwork_dict[sg][1].values():
            unique_cutoff_values |= set(c_l)
            
        subnetwork_dict[sg][0] = len(unique_cutoff_values)
    
    return subnetwork_dict