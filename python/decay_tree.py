# -*- coding: utf-8 -*-

import numpy as np

def build_decay_tree(parent, all_decays):
    my_decays = filter(lambda d: abs(d[2]) == abs(parent), all_decays)
    if len(my_decays) == 0:
        return (parent, None)
    else:
        return (parent, [(label, br, [build_decay_tree(child, all_decays) for child in children])
                         for (label, br, _, children) in my_decays])

def rescale_branch(branch, scaling_factor):
    label, br, subtrees = branch
    return (label, br / scaling_factor, subtrees)

def rescale_tree(tree, scaling_factor):
    top, branches = tree
    return (top, [rescale_branch(branch, scaling_factor) for branch in branches])

def optimize_branch(branch):
    label, initial_br, subtrees = branch
    # Rescale subtrees individually from each other
    scaling_factors, optimized_subtrees = zip(*[optimize_tree(subtree) for subtree in subtrees])
    # Rescale the branch by the maximum scaling factor of the trees,
    # to make sure all branching ratios remain below one
    branch_scaling_factor = max(scaling_factors)
    # Rescale the subtrees to a same scale, to preserve their relative proportions
    # (unless the weight is already zero, it which case it is not needed)
    scaled_subtrees = [rescale_tree(subtree, branch_scaling_factor / tree_scaling_factor)
                       if (subtree[1] is not None and tree_scaling_factor > 0) else subtree
                       for (tree_scaling_factor, subtree) in zip(scaling_factors, optimized_subtrees)]
    return (label, branch_scaling_factor * initial_br, scaled_subtrees)

def optimize_tree(tree):
    top, branches = tree
    # Leaf node: return the tree with a scaling factor of 1
    if branches is None:
        return (1, tree)
    # Otherwise, optimize each branch and rescale the tree
    else:
        optimized_branches = [optimize_branch(branch) for branch in branches]
        optimized_tree = (top, optimized_branches)
        branching_ratios = [branch[1] for branch in optimized_branches]
        total_branching_ratio = sum(branching_ratios)
        # Rescale the tree such that the branching ratios sum up to one
        # (unless it is zero, in which case we don't do anything because it would not make sense)
        rescaled_tree = (rescale_tree(optimized_tree, total_branching_ratio)
                         if total_branching_ratio > 0 else optimized_tree)
        # Return the tree and the scaling factor
        return (total_branching_ratio, rescaled_tree)

def branch_map(f, tree):
    _, branches = tree
    if branches is None:
        return f(None)
    else:
        return [(f((label, br, subtrees)), branch_map(f, subtree))
                for (label, br, subtrees) in branches
                for subtree in subtrees]

def node_map(f, tree):
    node, branches = tree
    if branches is None:
        return f(node)
    else:
        return [node_map(f, subtree)
                for (_, _, subtrees) in branches
                for subtree in subtrees]

def flatten_branching_ratios(tree):
    br, branches = tree
    lst = []
    if branches is None:
        lst.append(br)
    else:
        lst.extend(sum([[br * x for x in flatten_branching_ratios(branch)]
                        for branch in branches], []))
    return lst

def check_trees(original_tree, optimized_tree):
    original_branching_ratios = branch_map(
        lambda branch: branch[1] if branch is not None else None, original_tree)
    optimized_branching_ratios = (optimized_tree[0], branch_map(
        lambda branch: branch[1] if branch is not None else None, optimized_tree[1]))
    flattened_original_br = flatten_branching_ratios((1, original_branching_ratios))
    flattened_optimized_br = flatten_branching_ratios(optimized_branching_ratios)
    return all([np.isclose(x, y, rtol=1e-10) for (x, y) in list(zip(flattened_original_br, flattened_optimized_br))])
