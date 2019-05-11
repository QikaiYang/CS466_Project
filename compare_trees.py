"""
Comparison of two trees on different leaf sets
Copyright (c) 2018 NJMerge Developers
Erin K. Molloy <molloy.erin.k@gmail.com>
All rights reserved.
License: 3-Clause BSD,
see https://opensource.org/licenses/BSD-3-Clause
"""
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import sys
 
def compare_trees(tr1, tr2):
    # Find leaf labels that are in both trees
    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])
    com = lb1.intersection(lb2)
 
    # Restrict trees to shared leaf set
    if com != lb1 or com != lb2:
        com = list(com)
        tns = dendropy.TaxonNamespace(com)
 
        tr1.retain_taxa_with_labels(com)
        tr1.migrate_taxon_namespace(tns)
 
        tr2.retain_taxa_with_labels(com)
        tr2.migrate_taxon_namespace(tns)
    com = list(com)
 
    # Update tree bipartitions
    tr1.update_bipartitions()
    tr2.update_bipartitions()
 
    # Compute number of leaves and number of internal edges
    nl = len(com)
    ei1 = len(tr1.internal_edges(exclude_seed_edge=True))
    ei2 = len(tr2.internal_edges(exclude_seed_edge=True))
 
    # Compute number of false positives and false negatives
    [fp, fn] = false_positives_and_negatives(tr1, tr2)
 
    # Compute symmetric difference rate
    sd = float(fp + fn) / (ei1 + ei2)
 
    # Compute Robinson-Foulds error rate
    rf = float(fp + fn) / (2 * nl - 6)
 
    return(nl, ei1, ei2, fp, fn, sd, rf)
 
 
if __name__ == "__main__":
    argc = len(sys.argv)
    assert (argc == 3), "Usage: compare_trees.py tree1 tree2"
 
    tax = dendropy.TaxonNamespace()
    tr1 = dendropy.Tree.get(path=sys.argv[1],
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
 
    tr2 = dendropy.Tree.get(path=sys.argv[2],
                            schema='newick',
                            rooting='force-unrooted',
                            taxon_namespace=tax)
 
    # Unroot trees
    tr1.collapse_basal_bifurcation(set_as_unrooted_tree=True)
    tr2.collapse_basal_bifurcation(set_as_unrooted_tree=True)
 
    # Compute RF distance
    [nl, ei1, ei2, fp, fn, sd, rf] = compare_trees(tr1, tr2)
    sys.stdout.write('%d %d %d %d %d %f %f' % (nl, ei1, ei2, fp, fn, sd, rf))
    sys.stdout.flush()