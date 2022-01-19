#!/usr/bin/env python3
# #!/Users/asanto/anaconda3/bin/python
# #!/usr/bin/python
import sys, argparse
import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import CoordinationAnalysisModifier


   
def main(argv=None):
    # Parse in command-line arguments, and create the user help instructions
    parser = argparse.ArgumentParser(description='Average meshed-field density'
                                                 'and convert to 1-d density profile')
    parser.add_argument("--configfiles", type=str, required=True, nargs='+', help='configuration filename')
    parser.add_argument("--outputfile", type=str, required=True, help='filename with crystallinity as a function of time')
    parser.add_argument("--cutoff", type=float, required=True, help='filename with crystallinity as a function of time')
    parser.add_argument("--nbins", type=int, required=True, help='filename with crystallinity as a function of time')
    parser.add_argument("--ntypes", type=int, required=True, help='filename with crystallinity as a function of time')

    pipeline = import_file(parser.parse_args().configfiles)
    # Select type:
    
    # Coordination analysis:
    pipeline.modifiers.append(CoordinationAnalysisModifier(
        cutoff = parser.parse_args().cutoff, 
        number_of_bins = parser.parse_args().nbins, 
        partial = True))

    # Access the output DataTable:
    rdf_table = pipeline.compute().tables['coordination-rdf']
    
    # The y-property of the data points of the DataTable is now a vectorial property.
    # Each vector component represents one partial RDF.
    rdf_names = rdf_table.y.component_names
    
    ofile = open(parser.parse_args().outputfile, 'w')
    # Print a list of partial g(r) functions.
    ofile.write("# r g(r)")
    for component, name in enumerate(rdf_names):
        ofile.write("_%s " % name)
    ofile.write("\n")
    for ir in range(parser.parse_args().nbins):
        ofile.write("%e " % rdf_table.xy()[ir,0])
        for component, name in enumerate(rdf_names):
            ofile.write("%e " % rdf_table.y[ir,component])
            #ofile.write("%e " % rdf_table.xy()[ir,component-1])
        ofile.write("\n")
    ofile.close()

if __name__ == '__main__':
    sys.exit(main())
