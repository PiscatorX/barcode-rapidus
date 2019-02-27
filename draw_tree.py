#! /home/drewx/anaconda_ete/bin/python

from ete3 import PhyloTree, faces, TreeStyle
import subprocess
import argparse
import pprint

__author__ = "Andrew Ndhlovu"
__copyright__ = "Copyright 2018"
__license__ = "GPL"
__version__ = "3"
__maintainer__ = "Andrew Ndhlovu"
__email__ = "drewxdvst@outlook.com"
global code2name
code2name = None



class DrawTree(object):

    def __init__(self):

        parser = argparse.ArgumentParser(description="Draw phylogenetic tree")
        parser.add_argument('tree_file', action='store', type=str)
        parser.add_argument('-f','--img-format', dest='img_format', action='store', required=False,
                            default="png", type=str)
        parser.add_argument('-r','--ref', dest='ref_file', action='store',required=True, type=file)
        self.args = parser.parse_args()
        self.tree_bname =  self.args.tree_file.split('.')[0]
        self.tree =  PhyloTree(self.args.tree_file)
        self.img_format = self.args.img_format

        
    def get_ref(self):
        
        code2name = {}
        for line in self.args.ref_file.read().splitlines():
            data = tuple(line.split('\t'))[:2]
            if len(data) == 2:
                code2name[data[0]] = data[1]
        print code2name
        return code2name


    def draw(self, Ts):
        
        circular_style = TreeStyle()
        #circular_style.mode = "c" 
        circular_style.scale = 20
        img_fname =  '.'.join([self.tree_bname, self.img_format])
        self.tree.render(img_fname, tree_style = Ts, w=400, units="mm")
        #, tree_style=circular_style)
        try:
            subprocess.call(['gpicview',img_fname])
        except:
            pass
        
    
def node_aesthetic(node):
    
    #nameFace = faces.AttrFace("name", fsize=12, fgcolor="#009000")
    if node.is_leaf and node.name:
        if node.name in code2name:
            node.name = "[{0}] {1}".format(node.name, code2name.get(node.name, node.name))
        pprint.pprint(code2name)
        #print code2name[node.name]
        # faces.add_face_to_node(nameFace, node, column=0)
        # #We can also create faces on the fly
        # print node.name
        # longNameFace = faces.TextFace(code2name[node.name])
        # faces.add_face_to_node(longNameFace, node, column=0)

        # #text faces support multiline. We add a text face
        # #with the whole description of each leaf.
        #descFace = faces.TextFace(code2desc[node.name], fsize=10)
        # descFace.margin_top = 10
        # descFace.margin_bottom = 10
        # descFace.border.margin = 1

        # #Note that this faces is added in "aligned" mode
        # faces.add_face_to_node(descFace, node, column=0, aligned=True)

        #Sets the style of leaf nodes
        node.img_style["size"] = 12
        node.img_style["shape"] = "circle"


if __name__ == '__main__':
    phylo = DrawTree()
    code2name = phylo.get_ref()
    Ts = TreeStyle()
    Ts.layout_fn = node_aesthetic
    phylo.draw(Ts)
