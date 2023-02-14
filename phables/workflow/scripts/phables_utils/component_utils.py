def get_components(
    assembly_graph,
    unitig_names,
    smg_unitigs,
    unitig_phrogs,
    circular,
    edges_lengths,
    cicular_len,
):
    """
    Get connected components with PHROGs and no SMGs.
    """

    pruned_vs = {}

    i = 0

    for component in assembly_graph.components():

        if len(component) > 1:

            head_present = False
            connector_present = False
            tail_present = False
            lysis_present = False

            for unitig in component:

                if unitig_names[unitig] in smg_unitigs:
                    break
                elif unitig_names[unitig] in unitig_phrogs:

                    for protein in unitig_phrogs[unitig_names[unitig]]:
                        if "head and packaging" in protein:
                            head_present = True
                        if "connector" in protein:
                            connector_present = True
                        if "portal" in protein:
                            tail_present = True
                        if "lysis" in protein:
                            lysis_present = True

            if (head_present or connector_present or tail_present or lysis_present):
                pruned_vs[i] = component
                i += 1

        if len(component) == 1:

            if (
                unitig_names[component[0]] in unitig_phrogs
                and unitig_names[component[0]] in circular
                and edges_lengths[unitig_names[component[0]]] > cicular_len
            ):
                pruned_vs[i] = component
                i += 1

    return pruned_vs
