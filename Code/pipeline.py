import os
import cobra
from cobra import *
import metquest
from metquest import *
import networkx as nx
from networkx import *
import re
from re import *
import shutil
import pandas as pd
from itertools import permutations
from itertools import combinations
import argparse
from argparse import ArgumentParser


global artifacts
# biomass, DNA replication, RNA transcription, Protein biosynthesis
artifacts = [
    "cpd17041",
    "cpd17042",
    "cpd17043",
    "cpd11416"]


def format_seed(
        seed_dir,
        site_name,
        model_compartments_dict,
        model_id_list=[]):
    seed_list = open(
        os.path.join(
            seed_dir,
            site_name +
            ".txt"),
        'r').read().splitlines()
    cofactor_list=[]
    if(os.path.isfile(os.path.join(seed_dir, "cofactors.txt"))):
        cofactor_list = open(
            os.path.join(
                seed_dir,
                "cofactors.txt"),
            'r').read().splitlines()

    total = set()
    total = set(seed_list).copy().union(set(cofactor_list).copy())
    seed_set = set()
    for model_id in model_id_list:
        for compartment in model_compartments_dict[model_id]:
            for met in total:
                seed_set.add(model_id + " " + met + "_" + compartment)
    return seed_set


def get_model_info(model_dir):
    model_id_dict = {}
    model_compartments_dict = {}
    model_e0_dict={}
    for model_file in os.listdir(model_dir):
        model = cobra.io.read_sbml_model(os.path.join(model_dir, model_file))
        model_id_dict[model_file.replace(".xml", "")] = model.id
        model_compartments_dict[model.id] = list(model.compartments.keys())
        model_e0_dict[model.id]=[]

        for reaction in model.reactions:
            if(reaction.id.endswith("e0") and not(reaction.id.startswith("EX"))):
                '''
                flag=1
                for metabolite in list(reaction.metabolites.keys()):
                    if(not(metabolite.id.endswith("e0"))):
                        flag=0
                        break
                if(not(flag==0)):
                '''
                model_e0_dict[model.id].append(reaction.id)

    return model_id_dict, model_compartments_dict, model_e0_dict




def remove_artifacts(graph, namemap, model_compartments_dict, model_id_list=[]):
    global artifacts
    remove = set()
    for key, values in namemap.items():
        if("bio1" in values):
            remove.add(key)

    artifacts_set = set()
    for i in artifacts:
        for model_id in model_id_list:
            for j in model_compartments_dict[model_id]
                if(graph.has_node(i+"_"+j)):
                    artifacts_set.add(i+"_"+j)
                if(graph.has_node(model_id+"_"+j)):
                    artifacts_set.add(model_id + " " + i+"_"+j)
            

    for i in artifacts_set:
        for in_node in graph.out_edges(i):
            remove = remove.copy().union(set(in_node))
        for out_node in list(graph.in_edges(i)):
            remove = remove.copy().union(set(out_node))

    for remove_node in remove:
        graph.remove_node(remove_node)
        #print("----removing artifact: " + remove_node + "----")
    return graph


def remove_files(temp_dir):
    for file in os.listdir(temp_dir):
        os.remove(os.path.join(temp_dir, file))


def get_stuck(graph, namemap, status_dict, model_id_list,model_e0_dict):
    i = ''
    stuck = []
    for model_id in model_id_list:
        e0_list=model_e0_dict[model_id]
        for i in graph.nodes:
            if(i.startswith("Org_" + model_id) and (i not in status_dict.keys())):
                if((re.search("Org_(.*) (.*)", i).group(2).startswith('I') or re.search("Org_(.*) (.*)", i).group(2).startswith('R')) and (namemap[i] not in e0_list)):
                    stuck.append(i)
    return stuck


def get_individual_graph(model_dir, model_id_dict):
    print("----------Creating individual graphs----------")
    individual_graph_dict = {}
    individual_namemap_dict = {}

    for model_file in os.listdir(model_dir):
        model_file_name = re.search("(.*).xml", model_file).group(1)
        shutil.copyfile(
            os.path.join(
                model_dir, model_file), os.path.join(
                temp_dir, model_file))

        num = 1
        graph, namemap = create_graph(temp_dir, num)
        graph = remove_artifacts(
            graph, namemap,model_compartments_dict, [
                model_id_dict[model_file_name]])

        individual_graph_dict[model_id_dict[model_file_name]] = [graph]
        individual_namemap_dict[model_id_dict[model_file_name]] = [namemap]

        remove_files(temp_dir)

    return individual_graph_dict, individual_namemap_dict


def get_pair_graph(model_dir, model_id_dict, list_dir):
    print("----------Creating pair graphs----------")
    pair_graph_dict = {}
    pair_namemap_dict = {}
    done = []
    for list_file in os.listdir(list_dir):
        model_list = open(
            os.path.join(
                list_dir,
                list_file),
            'r').read().splitlines()

        combinations_list = combinations(model_list, 2)
        for comb in combinations_list:
            if(comb[0]+"_"+comb[1] not in done and comb[1] + "_" + comb[0] not in done):
                done.append(comb[0]+"_"+comb[1])
                pair = model_id_dict[comb[0]] + "_" + model_id_dict[comb[1]]
                shutil.copyfile(
                    os.path.join(
                        model_dir,
                        comb[0] +
                        ".xml"),
                    os.path.join(
                        temp_dir,
                        comb[0] +
                        ".xml"))
                shutil.copyfile(
                    os.path.join(
                        model_dir,
                        comb[1] +
                        ".xml"),
                    os.path.join(
                        temp_dir,
                        comb[1] +
                        ".xml"))
                num = 2
                graph, namemap = create_graph(temp_dir, num)
                graph = remove_artifacts(
                    graph, namemap, [
                        model_id_dict[comb[0]],model_id_dict[comb[1]]])
                pair_graph_dict[pair] = [graph]
                pair_namemap_dict[pair] = [namemap]

                remove_files(temp_dir)
    return pair_graph_dict, pair_namemap_dict


def get_community_graph(model_dir, model_id_dict, list_dir):
    print("----------Creating community graphs----------")
    community_graph_dict={}
    community_namemap_dict={}
    for list_file in os.listdir(list_dir):
        site_name = list_file.replace(".txt", "")

        model_list = open(
            os.path.join(
                list_dir,
                list_file)).read().splitlines()
        for model_file_name in model_list:
            shutil.copyfile(
                os.path.join(
                    model_dir,
                    model_file_name +
                    ".xml"),
                os.path.join(
                    temp_dir,
                    model_file_name +
                    ".xml"))

        # create a bi-partite graph, use the site name as graph name
        graph, namemap = create_graph(temp_dir, len(model_list))
        graph = remove_artifacts(
            graph, namemap,model_compartments_dict, [
                model_id_dict[model_file_name] for model_file_name in model_list])
        community_graph_dict[site_name] = [graph]
        community_namemap_dict[site_name] = [namemap]
        remove_files(temp_dir)
    return community_graph_dict, community_namemap_dict


def metabolic_support_index(root_dir, list_dir, seed_dir, model_id_dict, model_compartments_dict,model_e0_dict, individual_graph_dict, individual_namemap_dict, pair_graph_dict, pair_namemap_dict):
    print("----------Calculating metabolic support indices----------")
    msi_list=[]
    msi_list = [["Microorganism", "In the presence of", "MSI", "Site"]]
    sites = []

    for list_file in os.listdir(list_dir):
        sites.append(list_file)

    for list_file in sites:
        seed = set()
        model_list = open(
            os.path.join(
                list_dir,
                list_file)).read().splitlines()

        permutations_list = permutations(model_list, 2)

        for perm in permutations_list:

            individual_graph = individual_graph_dict[model_id_dict[perm[0]]][0]
            individual_namemap = individual_namemap_dict[model_id_dict[perm[0]]][0]
            if(model_id_dict[perm[0]] + "_" + model_id_dict[perm[1]] in pair_graph_dict):
                pair_graph = pair_graph_dict[model_id_dict[perm[0]
                                                           ] + "_" + model_id_dict[perm[1]]][0]
                pair_namemap = pair_namemap_dict[model_id_dict[perm[0]
                                                           ] + "_" + model_id_dict[perm[1]]][0]
            else:
                pair_graph = pair_graph_dict[model_id_dict[perm[1]
                                                           ] + "_" + model_id_dict[perm[0]]][0]
                pair_namemap = pair_namemap_dict[model_id_dict[perm[1]
                                                           ] + "_" + model_id_dict[perm[0]]][0]
            seed_pair = set()
            seed_pair = format_seed(seed_dir,
                                    list_file.replace(".txt",
                                                      ""),
                                    model_compartments_dict,
                                    [model_id_dict[perm[0]],
                                     model_id_dict[perm[1]]])

            seed_individual = set()
            seed_individual = format_seed(seed_dir, list_file.replace(
                ".txt", ""), model_compartments_dict, [model_id_dict[perm[0]]])

            lower_bound_individual, status_dict_individual, scope_individual = forward_pass(
                individual_graph, seed_individual)

            stuck_individual = get_stuck(individual_graph, individual_namemap, status_dict_individual, [
                                         model_id_dict[perm[0]]], model_e0_dict)

            lower_bound_pair, status_dict_pair, scope_pair = forward_pass(
                pair_graph, seed_pair)
            stuck_pair = get_stuck(pair_graph,pair_namemap, status_dict_pair, [
                                   model_id_dict[perm[0]]],model_e0_dict)

            MSI = (1 - (len(stuck_pair) / len(stuck_individual))) * 100
            print(MSI)

            msi_list.append([perm[0], perm[1], MSI, list_file.replace(".txt", "")])

    df = pd.DataFrame(msi_list)
    df.columns = df.iloc[0]
    df=df[1:]
    df.to_csv(root_dir + "/Metabolic Support Index.csv")
    print("Done")


def community_support_index_individual(
        root_dir,
        list_dir,
        seed_dir,
        model_id_dict,
        model_compartments_dict,
        model_e0_dict,
        individual_graph_dict,
        individual_namemap_dict,
        community_graph_dict,
        community_namemap_dict):
    print("----------Calculating community support indices to estimate the dependence of an individual on the community----------")
    csi_dict = {}
    for list_file in os.listdir(list_dir):
        site_name = list_file.replace(".txt", "")
        csi_dict[site_name] = {}
        model_list = open(
            os.path.join(
                list_dir,
                list_file)).read().splitlines()
        model_id_list = []
        for i in model_list:
            model_id_list.append(model_id_dict[i])

        community_graph = community_graph_dict[site_name][0]
        community_namemap = community_namemap_dict[site_name][0]

        seed_community = set()
        seed_community = format_seed(
            seed_dir,
            site_name,
            model_compartments_dict,
            model_id_list)
        lower_bound_community, status_dict_community, scope_community = forward_pass(
            community_graph, seed_community)


        for model in model_list:
            csi_dict[site_name][model] = []
            individual_graph=individual_graph_dict[model_id_dict[model]][0]
            individual_namemap=individual_namemap_dict[model_id_dict[model]][0]
            seed_individual = set()
            seed_individual = format_seed(
                seed_dir, site_name, model_compartments_dict, [
                    model_id_dict[model]])

            lower_bound_individual, status_dict_individual, scope_individual = forward_pass(
                individual_graph, seed_individual)

            stuck_individual = get_stuck(
                individual_graph, individual_namemap, status_dict_individual, [
                    model_id_dict[model]],model_e0_dict)
            stuck_community = []
            stuck_community = get_stuck(
                community_graph, community_namemap, status_dict_community, [
                    model_id_dict[model]],model_e0_dict)
            csi_dict[site_name][model] = (1 - (len(stuck_community) / len(stuck_individual))) * 100

    df = pd.DataFrame(csi_dict)
    df.to_csv(
        root_dir +
        "/Community Support Index_Dependence on the community.csv")
    print("Done")


def community_support_index_community(
        root_dir,
        list_dir,
        seed_dir,
        model_dir,
        model_id_dict,
        model_compartments_dict,
        model_e0_dict,
        community_graph_dict,
        community_namemap_dict):

    print("----------Calculating community support indices to estimate the dependence of the community on the individual----------")
    csi_dict = {}
    for list_file in os.listdir(list_dir):
        site_name = re.search("(.*).txt", list_file).group(1)
        csi_dict[site_name] = {}

        model_list = open(
            os.path.join(
                list_dir,
                list_file)).read().splitlines()
        model_id_list = []
        for i in model_list:
            model_id_list.append(model_id_dict[i])

        for model_out in model_list:
            csi_dict[site_name][model_out] = []
            model_in_list = []
            for model_in in model_list:
                if(model_in != model_out):
                    model_in_list.append(model_id_dict[model_in])
                    shutil.copyfile(
                        os.path.join(
                            model_dir,
                            model_in +
                            ".xml"),
                        os.path.join(
                            temp_dir,
                            model_in +
                            ".xml"))
            sub_graph, sub_namemap = create_graph(temp_dir, len(model_in_list))
            sub_graph = remove_artifacts(sub_graph, sub_namemap,model_compartments_dict, model_in_list)

            remove_files(temp_dir)

            sub_seed = format_seed(
                seed_dir, site_name, model_compartments_dict, model_in_list)
            lower_bound_sub, status_dict_sub, scope_sub = forward_pass(
                sub_graph, sub_seed)
            stuck_sub = []
            stuck_sub = get_stuck(sub_graph,sub_namemap,status_dict_sub, model_in_list,model_e0_dict)

            seed_community = set()
            seed_community = format_seed(
                seed_dir, site_name, model_compartments_dict, model_id_list)
            community_graph = community_graph_dict[site_name][0]
            community_namemap = community_namemap_dict[site_name][0]
            lower_bound_community, status_dict_community, scope_community = forward_pass(
                community_graph, seed_community)
            stuck_community = []
            stuck_community = get_stuck(
                community_graph, community_namemap, status_dict_community, model_in_list,model_e0_dict)

            csi_dict[site_name][model_out] = (1 - (len(stuck_community) / len(stuck_sub))) * 100

    df = pd.DataFrame(csi_dict)
    df.to_csv(root_dir + "/Community Support Index_Leave one out.csv")
    print("Done")


def cluster(
        root_dir,
        list_dir,
        seed_dir,
        model_dir,
        model_id_dict,
        model_compartments_dict,
        model_e0_dict,
        community_graph_dict,
        community_namemap_dict):
    print("----------Calculating community support indices for family-level analyses----------")
    cluster_path=os.path.join(root_dir,"cluster.csv")
    cluster_dict = {}
    for key, value in pd.read_csv(
        cluster_path,
        index_col=0,
        header=0).to_dict()['Cluster'].items():
        if(value in cluster_dict):
            cluster_dict[value].append(key)
        else:
            cluster_dict[value] = [key]

    cluster_loo_dict = {}

    for site in os.listdir(list_dir):
        #get F-L name
        site_name = re.search("(.*).txt", site).group(1)
        #get community graph for site
        community_graph = community_graph_dict[site_name][0]
        community_namemap = community_namemap_dict[site_name][0]
        #Initialise dictionary
        cluster_loo_dict[site_name] = {}
        #get all model names in site
        model_list = open(os.path.join(list_dir, site)).read().splitlines()
        #get model id
        model_id_list = []
        for i in model_list:
            model_id_list.append(model_id_dict[i])

        # ................................................

        #get list of models to keep in for each cluster
        loo_dict = {}
        loo_cluster_no_set=set()
        for i in model_list:
            loo_cluster_no_set.add(pd.read_csv(cluster_path,index_col=0,header=0).to_dict()['Cluster'][i])

        for cluster_no, cluster_list in cluster_dict.items():
            if(cluster_no in loo_cluster_no_set):
                loo_dict[cluster_no] = []
                for model in model_list:
                    if model not in cluster_list:
                        loo_dict[cluster_no].append(model)

        for cluster_in, cluster_list_in in loo_dict.items():
            for model in cluster_list_in:
                shutil.copyfile(
                    os.path.join(
                        model_dir,
                        model + ".xml"),
                    os.path.join(
                        temp_dir,
                        model + ".xml"))

            model_id_list_loo = []
            for i in cluster_list_in:
                model_id_list_loo.append(model_id_dict[i])

            graph_loo, namemap_loo = create_graph(
                temp_dir, len(cluster_list_in))

            graph_loo = remove_artifacts(
                graph_loo, namemap_loo,model_compartments_dict, model_id_list_loo)
            seed_loo = format_seed(
                seed_dir,
                site_name,
                model_compartments_dict,
                model_id_list_loo)
            lower_bound_loo, status_dict_loo, scope_loo = forward_pass(
                graph_loo, seed_loo)
            stuck_loo = []
            stuck_loo = get_stuck(
                graph_loo, namemap_loo,status_dict_loo, model_id_list_loo,model_e0_dict)

            seed_community = format_seed(
                seed_dir, site_name, model_compartments_dict, model_id_list)
            lower_bound_community, status_dict_community, scope_community = forward_pass(
                community_graph, seed_community)
            stuck_community_loo = []
            stuck_community_loo = get_stuck(
                community_graph,community_namemap ,status_dict_community, model_id_list_loo,model_e0_dict)

            cluster_loo_dict[site_name][cluster_in] = (
                1 - (len(stuck_community_loo) / len(stuck_loo))) * 100
            # print(stuck_community,len(stuck_loo))
            remove_files(temp_dir)

    df = pd.DataFrame(cluster_loo_dict)
    # df.columns = df.iloc[0]
    df.to_csv(root_dir + "/Community Support Index leave family out.csv")
    print("Done")


def maincall(root_dir, model_dir, seed_dir, list_dir, msi, csi):
    global temp_dir
    temp_dir = os.path.join(root_dir, "temp/")
    if(not(os.path.isdir(temp_dir))):
        os.mkdir(temp_dir)

    print("----------Please comment lines 212 to 216 in the create_graph function of the construct_graph module----------")

    model_id_dict, model_compartments_dict,model_e0_dict = get_model_info(model_dir)
    if(msi and not(csi)):
        individual_graph_dict, individual_namemap_dict = get_individual_graph(
            model_dir, model_id_dict)
        pair_graph_dict, pair_namemap_dict = get_pair_graph(
            model_dir, model_id_dict, list_dir)
        msi(root_dir,
            list_dir,
            seed_dir,
            model_id_dict,
            model_compartments_dict,
            individual_graph_dict,
            pair_graph_dict)

    if(csi and not(msi)):
        community_graph_dict, community_namemap_dict = get_community_graph(
            model_dir, model_id_dict, list_dir)
        community_support_index_individual(
            root_dir,
            list_dir,
            seed_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            individual_graph_dict,
            individual_namemap_dict,
            community_graph_dict,
            community_namemap_dict)
        community_support_index_community(
            root_dir,
            list_dir,
            seed_dir,
            model_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            community_graph_dict,
            community_namemap_dict)
        cluster(
            root_dir,
            list_dir,
            seed_dir,
            model_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            community_graph_dict,
            community_namemap_dict)

    if(msi and csi):
        individual_graph_dict, individual_namemap_dict = get_individual_graph(
            model_dir, model_id_dict)
        pair_graph_dict, pair_namemap_dict = get_pair_graph(
            model_dir, model_id_dict, list_dir)
        metabolic_support_index(root_dir,
            list_dir,
            seed_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            individual_graph_dict,
            individual_namemap_dict,
            pair_graph_dict,
            pair_namemap_dict)
        community_graph_dict, community_namemap_dict = get_community_graph(
            model_dir, model_id_dict, list_dir)
        community_support_index_individual(
            root_dir,
            list_dir,
            seed_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            individual_graph_dict,
            individual_namemap_dict,
            community_graph_dict,
            community_namemap_dict)
        community_support_index_community(
            root_dir,
            list_dir,
            seed_dir,
            model_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            community_graph_dict,
            community_namemap_dict)
        cluster(
            root_dir,
            list_dir,
            seed_dir,
            model_dir,
            model_id_dict,
            model_compartments_dict,
            model_e0_dict,
            community_graph_dict,
            community_namemap_dict)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-r",
        "--rootdir",
        help="Enter the directory for analyses")
    parser.add_argument(
        "-m",
        "--modeldir",
        help="Enter the directory with models")
    parser.add_argument(
        "-s",
        "--seeddir",
        help="Enter the directory to the seed file")
    parser.add_argument("-l", "--listdir", help="Enter list directory")
    parser.add_argument(
        "-csi",
        "--community",
        help="Community support indices",
        action='store_true')
    parser.add_argument(
        "-msi",
        "--pair",
        help="Metabolic support indices",
        action='store_true')

    args = parser.parse_args()

    maincall(
        root_dir=args.rootdir,
        model_dir=args.modeldir,
        seed_dir=args.seeddir,
        list_dir=args.listdir,
        msi=args.pair,
        csi=args.community)


if __name__ == '__main__':
    main()
