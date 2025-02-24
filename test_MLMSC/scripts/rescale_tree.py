"""Rescale trees for some c value for dupcoal"""

import argparse
import dendropy

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=str, help="Path to input tree")
    parser.add_argument("--scale", required=True, type=float, help="scaling factor")
    parser.add_argument("--output", required=True, type=str, help="Path to write results")
    return parser.parse_args()

def scale_tree(args):
    t = dendropy.Tree.get(path=args.input, schema="newick")
    for edge in t.preorder_edge_iter():
        if(edge.length == None):
            continue
        edge.length /= args.scale
    t.write(path=args.output, schema="newick")


def main():

    args = parse_arguments()
    scale_tree(args)

if __name__ == "__main__":
    main()
