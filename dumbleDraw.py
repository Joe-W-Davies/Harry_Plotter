#!/usr/bin/env python

from optparse     import OptionParser
import HarryPlotter as HP

#get options
def get_options():
    parser = OptionParser()
    parser.add_option('-p', '--plotcard', dest='plot_card',default='plotcard.py',
                      help=''' Load the plot card in .py/txt format. ''')

    parser.add_option('-o', "--outdir", dest='out_dir', default='plots',
                      help='Path you want to write the output histograms to')

    parser.add_option("-s", "--sampledir", dest="sample_dir",default='./data/',
                      help='''
                      Specify the detrectory where the trees are.
                      e.g.: --sampledir /data/trees
                      ''')
    return parser.parse_args()




#wont work since it needs to be invoken on an object

# NB reader uses objects from Sample, Options, Variable, Systematics

if __name__ == '__main__':
    (opt, args) = get_options()
  #STEPS:
  #get plot card: from options
    with open(opt.plot_card) as card: 
    #FIXME: what is the best way to fill as class attributes using a dict?? 
    #FIXME: esp if we want some attributes to just be their deafult values

    # HP.Reader.read
        #should: 
        #        - read variables options [function done]
        #        - read other options [function done]
        reader = HP.Reader(card.read())
        reader.read()
        #        - get trees for all samples and put them in df's []. takes account of concatting df's across years if we want to merge samples
        reader.trees_to_dfs()
        # get trees for systematics and hold them in dicts. Takes account of concatting df's across years if we want to merge samples
        reader.systs_to_dfs()


    #HP.plotting
        #shoud:  - perform the plotting based on the options
 
