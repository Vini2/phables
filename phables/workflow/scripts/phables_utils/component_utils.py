def get_components(
    assembly_graph,
    unitig_names,
    smg_unitigs,
    unitig_phrogs,
    circular,
    edges_lengths,
    cicular_len,
    phrog_dict,
):
    """
    Get connected components with PHROGs and no SMGs.

    Category checks are exact-match against phrog_dict[phrog] (gene_utils.get_phrog_unitigs
    now stores just the category, not concatenated with the product/annot text). Used to be
    a substring test, which false-matched real data: e.g. phrog_5's product is "tail
    completion or Neck1 protein" with category "connector", and the word "tail" in that
    product text satisfied `"tail" in phrog_dict[phrog]` even though the PHROG isn't a tail
    gene. 110 PHROGs hit this in the current annotation table (103 false "tail" matches, 7
    false "connector" matches) -- mostly connector-category genes whose product names
    reference the tail apparatus they connect to.
    """

    pruned_vs = {}
    likely_complete = {}

    i = 0

    comp_phrogs = {}

    for component in assembly_graph.components():
        phrogs_found = set()

        head_present = False
        connector_present = False
        tail_present = False
        lysis_present = False

        # Check every unitig in the component for a bacterial single-copy marker
        # gene before scoring any PHROG evidence. This has to be a separate,
        # unconditional pass over the whole component: doing the SMG check and
        # the PHROG-category scan in the same loop with a `break` on SMG only
        # stops iterating -- it doesn't undo category flags a unitig earlier in
        # the same component already set, so whether a mixed component got
        # excluded ended up depending on iteration order, not on whether an SMG
        # was actually present. Skip PHROG scoring entirely once any SMG is
        # found; there's no point computing it for a component that's excluded
        # either way.
        has_smg = any(unitig_names[unitig] in smg_unitigs for unitig in component)

        if has_smg:
            continue

        if len(component) > 1:
            for unitig in component:
                if unitig_names[unitig] in unitig_phrogs:
                    for phrog in unitig_phrogs[unitig_names[unitig]]:
                        if phrog_dict[phrog] == "head and packaging":
                            head_present = True
                        if phrog_dict[phrog] == "connector":
                            connector_present = True
                        if phrog_dict[phrog] == "tail":
                            tail_present = True
                        if phrog_dict[phrog] == "lysis":
                            lysis_present = True

                        phrogs_found.add(phrog)

            if head_present or connector_present or tail_present or lysis_present:
                pruned_vs[i] = component
                comp_phrogs[i] = phrogs_found
                i += 1

        if len(component) == 1:
            unitig = component[0]
            phrogs_present = False

            if unitig_names[unitig] in unitig_phrogs:
                for phrog in unitig_phrogs[unitig_names[unitig]]:
                    if phrog_dict[phrog] == "head and packaging":
                        head_present = True
                    if phrog_dict[phrog] == "connector":
                        connector_present = True
                    if phrog_dict[phrog] == "tail":
                        tail_present = True
                    if phrog_dict[phrog] == "lysis":
                        lysis_present = True

                    phrogs_found.add(phrog)

            # Check PHROG categories in unitig (should contain at least one)
            if head_present or connector_present or tail_present or lysis_present:
                phrogs_present = True

            if phrogs_present and edges_lengths[unitig_names[unitig]] > cicular_len:
                pruned_vs[i] = component
                comp_phrogs[i] = phrogs_found

                # Check if all PHROG categories are present in unitig
                if (
                    head_present
                    and connector_present
                    and tail_present
                    and lysis_present
                ):
                    likely_complete[i] = 1
                else:
                    likely_complete[i] = 0

                i += 1

    return pruned_vs, comp_phrogs, likely_complete
