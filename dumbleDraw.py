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

    parser.add_option('-a', "--all", dest='draw_all', default=False,
                      help='draw all variable specified on plot card')

    parser.add_option('-v', "--variable", dest='variable', default='',
                      help='name of a single variable to be drawn')
   
    return parser.parse_args()




#wont work since it needs to be invoken on an object

# NB Plotter uses objects from Sample, Options, Variable, Systematics

if __name__ == '__main__':
    (opt, args) = get_options()
  #STEPS:
  #get plot card: from options
    with open(opt.plot_card) as card: 
    # HP.Plotter.read
        #should: 
        #        - read variables options [function done]
        #        - read other options [function done]
        plotter = HP.Plotter(card.read(), output_dir= opt.out_dir )
        plotter.read()
        #        - get trees for all samples and put them in df's []. takes account of concatting df's across years if we want to merge samples
        plotter.trees_to_dfs()
        # get trees for systematics and hold them in dicts. Takes account of concatting df's across years if we want to merge samples
        plotter.systs_to_dfs()


    #Plotter.plot
        #shoud:  - perform the plotting based on the options
        if opt.draw_all:
            for var_name in plotter.var_names:
                plotter.draw(var_name)
            
        else: 
            plotter.draw(opt.variable)

