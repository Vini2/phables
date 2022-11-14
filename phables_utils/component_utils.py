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

            has_phrog = False

            for unitig in component:

                if unitig_names[unitig] in smg_unitigs:
                    break
                elif unitig_names[unitig] in unitig_phrogs:
                    if (
                        "terminase large subunit head and packaging"
                        in unitig_phrogs[unitig_names[unitig]]
                        or "capsid assembly protein tail"
                        in unitig_phrogs[unitig_names[unitig]]
                        or "portal protein head and packaging"
                        in unitig_phrogs[unitig_names[unitig]]
                    ):
                        has_phrog = True

            if has_phrog:
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
