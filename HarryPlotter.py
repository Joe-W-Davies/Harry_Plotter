import json
import scipy.stats
from os import system, path    
import math
import numpy as np
import pandas as pd
import uproot as upr
import matplotlib.pyplot as plt
import warnings
try:
    plt.style.use("cms10_6_HP")
except IOError:
    warnings.warn('Could not import user defined matplot style file. Using default style settings...') 

                               
class Variable(object):              
    '''
    hold plotting information inside variable object
    '''

    def __init__(self, name, options={}):
        self.template = {
            'name'       : '',
            'bin_range'  : '(0,100)',
            'num_bins'   : '30',
            'cuts'       : '',
            'lumi'       : 1.0,
            'xlabel'     : '',
            'log'        : False,
            'norm'       : False,
            'ratio'      : False
        }
        self.__dict__ = self.template
        self.__dict__.update(options)
        self.name     = name
        self.xlabel   = self.xlabel.replace('#', "\\")

    #def __repr__(): pass

#--------------------------------------------------------------------------------------
    
        
class Sample(object):
    '''
    object for sample attributes
    '''

    def __init__(self, name, options={}):
        self.template={
        'file_path'      : '',
        'tree'           : '',
        'order'          : None,
        'label'          : '',
        'lumi'           : 1.0,
        'scale'          : 1.0,
        'year'           : None,
        'file_ext'       : '.root',
        'sample_set'     : '',
        'colour'         : 'red'
        }
        self.__dict__    = self.template
        self.__dict__.update(options)
        self.name        = name
        self.label       = self.label.lower()
        self.systematics = {} #dict of {syst:syst_obj}

        if self.file_path.endswith('/'): self.file_path = self.file_path[:-1]
        
    #def __repr__(): pass



#--------------------------------------------------------------------------------------

class Systematic(object):
    '''Object containing atrtibues related to systematic variations
       One object is created per systematic uncertainty, so that one
       sample that may hold multiple variations contained in separate objects.
    '''
    def __init__(self, name, down_tree='', up_tree=''):
            self.name                  = name
            self.up_tree               = up_tree
            self.down_tree             = down_tree
            self.up_frame              = None
            self.down_frame            = None
            self.up_tree_binned        = [] #[true_downs, true_ups]
            self.down_tree_binned      = [] #[true_downs, true_ups]
            self.true_up               = {}
            self.true_down             = {}
    #def __repr__(): pass

#--------------------------------------------------------------------------------------

class Options(object):
    '''Object containing the plotting options'''
    def __init__(self, options={}):
        self.template = {
            'ratio_plot'         : False,
            'ratio_range'        : (0,2),
            'global_labels'      : '',
            'cms_label'          : 'Preliminary',
            'energy_label'       : '(13 TeV)',
            'concat_years'       : False,
            'reweight_var'       : '',
            'reweight_samp_mc'   : '',
            'reweight_samp_data' : '',
            'stack_bkgs'         : False

        }
        self.__dict__  = self.template
        self.__dict__.update(options)

    #def __repr__(): pass

    def print_options(self, title_length, divider):
        print '|' + 'Options'.center(title_length) + '|'
        print divider
        for option, value in self.__dict__.iteritems():
    	    print '|' + (option+' = '+str(value)).ljust(title_length) + '|'

        #for option, value in self.__dict__.iteritems()
        #    print 

      
#--------------------------------------------------------------------------------------

class Plotter(object):
    '''Class to read the plot card'''

    def __init__(self, card, output_dir=''):
        self.card                = card 
        self.variables           = {}
        self.samples             = {} #all samples in { key:df } format
        self.options             = Options()
        self.ratio_range         = (0,2) #FIXME: could make changeable for each variable? (add att to var)
        self.cut_map             = {} 
        self.output_dir          = output_dir
        self.legend              = []
        self.var_names           = []
        self.sample_sets         = set() #unique set of name that span years i.e. processes
        self.sample_years        = set() 
        self.sample_labels       = set() 
        self.sample_lumis        = set() 
        self.signal_df_map       = {} 
        self.bkg_df_map          = {} 
        self.data_df_map         = {} 
        self.rew_scale_factors   = []
        self.rew_bin_edges       = np.array([])
        self.colour_sample_map   = {} 
        #FIXME: can probably put some of these attribues inside one systematic class object for each sample?
        self.syst_tree_name_map  = {}
        self.mc_totals           = {} 
        self.mc_stat_uncs        = {} 

        self.check_dir(output_dir)


    def read(self):
        config = json.loads(self.card)
        for key in config: 
            if 'variables' in key:
                  for var_name, var_options in config[key].iteritems():
                      var = Variable(var_name, var_options)
                      self.variables[var.name] = var
                      self.var_names.append(var.name)
                      self.cut_map[var.name] = var.cuts
            if 'options' in key:
                opts = Options(config[key])
                self.options = opts
            if 'systematics' in key:
                for syst_name, syst_options in config[key].iteritems():
                    self.syst_tree_name_map[syst_name] = [syst_options['down_tree'], syst_options['up_tree']]
            if 'samples' in key:
                samples = {}
                for sample_name, sample_opts in config[key].iteritems():
                    sample = Sample(sample_name, sample_opts)
                    samples[sample_name]  = sample
                    self.sample_sets.add(sample.sample_set)
                    self.sample_years.add(sample.year)
                    self.colour_sample_map[sample.sample_set] = sample.colour
                    self.sample_lumis.add(float(sample.lumi))
                self.samples = samples

        print self.sample_lumis
    def trees_to_dfs(self):
        """
        Function to read in files for each sample and convert to a DataFrame.
        If concat across years, hold dataframes in dict that looks like this: 
        {'sample_set_1' : [df_for_year_1, df_for_year_2, ... ], 'sample_set_2' : [...], ...}
        The lists will then be concattenated across years i.e. flattened out
        Otherwise just store with one key per sample:
        {'sample_id_1' : df_1, 'sample_id_2' : ..., ...}
        Dictionaries are split into signal, background, and Data
        """

        sig_df_map  = {}
        bkg_df_map  = {}
        data_df_map = {}

        if self.options.concat_years: 
            for sample_name, sample_obj in self.samples.iteritems():
                #dict keys are unique so dont really need nested "if" statement but nice for thinking about it
                if   sample_obj.label == 'signal'     :
                    if sample_obj.sample_set not in sig_df_map.keys() : sig_df_map[sample_obj.sample_set]  = [] 
                elif sample_obj.label == 'background' : 
                    if sample_obj.sample_set not in bkg_df_map.keys() : bkg_df_map[sample_obj.sample_set]  = []
                elif sample_obj.label == 'data'       : 
                    if sample_obj.sample_set not in bkg_df_map.keys() : data_df_map[sample_obj.sample_set] = []

                else: raise Exception ("Got incorrect label '%s', for sample: %s. Accepted labels are: 'signal', 'background', or 'data'" ) % (sample_obj.label, sample_name)
        
        for sample_name, sample_obj in self.samples.iteritems():
	    if not path.isdir(sample_obj.file_path): 
              raise Exception('%s directory does not exist!'%sample_obj.file_path)
            input_file = upr.open( '%s/%s%s' % (sample_obj.file_path, sample_obj.name, sample_obj.file_ext) )
            input_tree = input_file[sample_obj.tree]
            input_df   = input_tree.pandas.df(self.var_names+['weight'])

            if sample_obj.label.lower() != 'data' : input_df['weight']*=float(sample_obj.lumi)
            if len(self.options.reweight_var)==0: input_df['weight']*=float(sample_obj.scale)
          
            
            if self.options.concat_years:
                if sample_obj.label.lower()== 'signal'     : sig_df_map[sample_obj.sample_set].append(input_df)
                elif sample_obj.label.lower()=='background': bkg_df_map[sample_obj.sample_set].append(input_df)
                elif sample_obj.label.lower()=='data'      : data_df_map[sample_obj.sample_set].append(input_df)
                else: raise Exception ("Got incorrect label '%s', for sample: %s. Accepted labels are: 'signal', 'background', or 'data'" ) % (sample_obj.label, sample_name)
            else:
                if sample_obj.label.lower() == 'signal'      : sig_df_map[sample_obj.name] = input_df
                elif sample_obj.label.lower() == 'background': bkg_df_map[sample_obj.name] = input_df
                elif sample_obj.label.lower() == 'data'      : data_df_map[sample_obj.name] = input_df
                else: raise Exception ("Got incorrect label '%s', for sample: %s. Accepted labels are: 'signal', 'background', or 'data'" ) % (sample_obj.label, sample_name)

        if self.options.concat_years: 
            for df_dict in [sig_df_map, bkg_df_map, data_df_map]:
                for sample_set_name, sample_set_list in df_dict.iteritems():
                    try:  assert( len(self.sample_years) == len(sample_set_list) )
                    except: raise Exception('Number of years (%i) not equal to number of sub-samples (%i) in sample set: %s' %(len(self.sample_years), len(sample_set_list), sample_set_name) )
                    df_dict[sample_set_name] = pd.concat(sample_set_list)

        #apply reweighting (bkg->data) once all inputs have been read in
        #FIXME: accomodate stacking of bkgs i.e. derive the re-weighting from the stacked bkgs

        if len(self.options.reweight_var)!=0: 
            #derive SFs in samples and variable specified in options:
            scale_factors, bin_edges = self.derive_scale_factors(self.options.reweight_var,
                                                                 bkg_df_map[self.options.reweight_samp_mc], 
                                                                 data_df_map[self.options.reweight_samp_data])
            self.rew_scale_factors = scale_factors
            self.rew_bin_edges     = bin_edges

            #apply scale factors to all bkgs in the dict
            for sample_set_name, sample_df in bkg_df_map.iteritems():
                bkg_df_map[sample_set_name] = self.apply_reweighting( self.options.reweight_var,
                                                                      self.rew_scale_factors, 
                                                                      self.rew_bin_edges,
                                                                      sample_df)

        self.sig_df_map  = sig_df_map
        self.bkg_df_map  = bkg_df_map
        self.data_df_map = data_df_map

    def systs_to_dfs(self):
        """
        Function to read in systs for each sample and convert to a DataFrame. 
        Only reads systematics for processes labelled with background.


        The function creates final data structure of the samples attribute for this Reader class, as:
        self.samples = { 'sample_name_1': sample_object,
                         'sample_name_2': sample_object, ... } 

        Each sample object contains a "systematics" attribute, which is another dictionary with structure:

        sample.systematics = {'syst_1_name': syst_object, 'syst_2_name': syst_object } 
         
        where each syst_object includes both up/down frames for a single source of systematic variation.
        """

        for sample_obj in self.samples.values():
            if sample_obj.label=='background':
                for syst in self.syst_tree_name_map.keys():
                    sample_obj.systematics[syst] = Systematic(syst, self.syst_tree_name_map[syst][0], self.syst_tree_name_map[syst][1])

        for sample_name, sample_obj in self.samples.iteritems():
            if sample_obj.label=='background':
                input_file = upr.open( '%s/%s%s' % (sample_obj.file_path, sample_obj.name, sample_obj.file_ext) )
                for syst, syst_obj in sample_obj.systematics.iteritems():
                    print 'reading systematic %s for sample %s' %(sample_name, syst)
                    up_tree            = input_file[syst_obj.up_tree]
                    up_df              = up_tree.pandas.df(self.var_names+['weight'])
                    up_df['weight']    *= float(sample_obj.lumi) 
                    if len(self.options.reweight_var)!=0: 
                        up_df = self.apply_reweighting( self.options.reweight_var,
                                                        self.rew_scale_factors, 
                                                        self.rew_bin_edges,
                                                        up_df)
                    else: up_df['weight'] *= float(sample_obj.scale)
                    syst_obj.up_frame = up_df

                    down_tree          = input_file[syst_obj.down_tree]
                    down_df            = down_tree.pandas.df(self.var_names+['weight'])
                    down_df['weight'] *= float(sample_obj.lumi) 
                    if len(self.options.reweight_var)!=0: 
                        down_df = self.apply_reweighting(self.options.reweight_var,
                                                         self.rew_scale_factors,
                                                         self.rew_bin_edges,
                                                         down_df)
                    else: down_df['weight'] *= float(sample_obj.scale)
                    syst_obj.down_frame = down_df


        #idea here, to account for merging systematics, is to loop through nominal sample objects
        # get the systematics for each, and concat across sample sets. This is done by putting the indiidual
        #sample object in nested dictioanries, concatting, then re-filling the Plotter.samples attribute with
        # a new sample object for each set, taking care to delete the individual ones

        if self.options.concat_years:
            # set up dict to hold indiv dataframes for sets corresponding to "background" sample types
            set_to_syst_map = {}
            for sample_obj in self.samples.values(): 
                    if (sample_obj.label=='background') and (sample_obj.sample_set not in set_to_syst_map.keys()):
                        set_to_syst_map[sample_obj.sample_set] = {} #systematics dict for each set
                        for syst in self.syst_tree_name_map.keys():
                            set_to_syst_map[sample_obj.sample_set][syst+'Up'] = []
                            set_to_syst_map[sample_obj.sample_set][syst+'Down'] = []

            #fill data frames lists based off sample sets
            for sample_obj in self.samples.values():
                if sample_obj.label=='background':
                    for syst in self.syst_tree_name_map.keys():
                        set_to_syst_map[sample_obj.sample_set][syst+'Down'].append( sample_obj.systematics[syst].down_frame )
                        set_to_syst_map[sample_obj.sample_set][syst+'Up'].append( sample_obj.systematics[syst].up_frame )
                    del self.samples[sample_obj.name]

            #concat the list of df's for each systematic i.e. merge across years and put into a new systematic object
            for set_name, syst_dict in set_to_syst_map.iteritems():
                #only bkg objects in the dict atm so can set label to background
                set_object = Sample(set_name)
                set_object.label = 'background'
                for syst in self.syst_tree_name_map.keys():
                    set_object.systematics[syst] = Systematic(syst)
                    set_object.systematics[syst].down_frame = pd.concat(syst_dict[syst+'Down'])
                    set_object.systematics[syst].up_frame = pd.concat(syst_dict[syst+'Up'])

                self.samples[set_name] = set_object

    #so the final data structure of the samples attribute for this Reader class is:
    #self.samples = { 'sample_name_1': sample_object },
    #                 'sample_name_2': sample_object } 

    #the sample object contains an "systematics" attribute, which is a dictionary with structure:

    #FIXME: can maybe merge sample objects before r
    #sample.systematics = {'syst_1_name': syst_object, 'syst_2_name': syst_object } 
     
    # where each syst_object includes both up/down frames for a single source of systematic variation
    # and may be merged across years if specified (i.e. sample-> set where both are a Sample() object )

    #--------------------------------------------------------------

    def draw(self, var_key):
        '''
        Main function to do the drawing, called once per variable.
        '''

        variable   = self.variables[var_key]
        print 'drawing var: %s' % variable.name
        cut_string = self.assemble_cut_string(variable.name)
        bins       = np.linspace( eval(variable.bin_range)[0], eval(variable.bin_range)[1],
                                 variable.num_bins)

        if self.options.ratio_plot: 
            fig, axes = plt.subplots(nrows=2, ncols=1, dpi=200, sharex=True,
                                          gridspec_kw ={'height_ratios':[3,0.8], 'hspace':0.08})    
            axes = self.set_canv_style(axes, variable, bins)
        else:
            fig  = plt.figure()
            axes = fig.gca()

        if len(self.sig_df_map) != 0: 
            self.plot_signals(cut_string, axes, variable, bins)

        if len(self.data_df_map) != 0: 
            data_binned, bin_centres, data_stat_down_up = self.plot_data(cut_string, axes, variable, bins)

        if len(self.bkg_df_map) != 0: 
            self.plot_bkgs(cut_string, axes, variable, bins, data_binned, bin_centres, data_stat_down_up)

        if len(self.syst_tree_name_map.keys()) != 0: 
            self.plot_systematics(cut_string, axes, variable, bins)

        #default axes dont have enough whitespace above plots, so add more. 
        #also set logy here and not when filling hists, or else you get weird behaviour if y~0
        if self.options.ratio_plot: 
            axes = axes[0]
        axes.legend(loc='upper right', bbox_to_anchor=(0.94,0.94)) 
        canv_bottom, canv_top = axes.get_ylim()
        axes.legend(loc='upper right', bbox_to_anchor=(0.94,0.94)) 
        if variable.log:
            axes.set_ylim(top= canv_top*100, bottom=canv_bottom*0.01) #log plot
            axes.set_yscale('log', nonposy='clip')    
        else: axes.set_ylim(top= canv_top*1.25)

        fig.savefig('%s/%s.pdf' % (self.output_dir, variable.name))

    #--------------------------------------------------------------------

    #FIXME: first parts of each of the 3 functions below can be
    #FIXME  put into a single function?
    #FIXME  can probably refactor to avoid duplicating code 3 times

    def plot_signals(self, cut_string, axes, variable, bins):
      
        for sample_id, sig_frame in self.sig_df_map.iteritems():
            sig_frame_cut  = sig_frame.query(cut_string).copy()
            var_to_plot    = sig_frame_cut[variable.name].values
            var_weights    = sig_frame_cut['weight'].values

            sumw, _            = np.histogram(var_to_plot, bins=bins, weights=var_weights)
            if variable.norm: var_weights /= np.sum(sumw)

            #FIXME: check colours are consistent for the sample!.. or do this earlier
            #FIXME when being read in, and append it as sampe attribute 

            #normed option is depracated so do manually
            if variable.norm: var_weights /= np.sum(sumw)

            if self.options.ratio_plot:
                axes[0].hist(var_to_plot, bins=bins, label=sample_id, weights=var_weights, 
                             color=self.colour_sample_map[sample_id], histtype='step')
            else:
                axes.hist(var_to_plot, bins=bins, label=sample_id, weights=var_weights, 
                          color=self.colour_sample_map[sample_id], histtype='step')
            del sig_frame_cut

    #--------------------------------------------------------------

    def plot_data(self, cut_string, axes, variable, bins):
        for sample_id, data_frame in self.data_df_map.iteritems():
            data_frame_cut     = data_frame.query(cut_string).copy()
            var_to_plot        = data_frame_cut[variable.name].values
            var_weights        = data_frame_cut['weight'].values

            data_binned, bin_edges = np.histogram(var_to_plot, bins=bins, weights=var_weights)
            print '--> Integral of hist: %s, for sample: %s, is: %.2f' % (variable.name,sample_id,np.sum(data_binned))
            bin_centres = (bin_edges[:-1] + bin_edges[1:])/2
            data_stat_down, data_stat_up = self.poisson_interval(data_binned, data_binned)

            if variable.norm:
                data_binned     /= np.sum(data_binned)
                data_stat_down  /= np.sum(data_binned)
                data_stat_up    /= np.sum(data_binned)

            #FIXME: sort this niche issue out
            #dataUp[dataUp==np.inf] == 0

            data_binned[data_binned==0] = np.nan #removes markers at zero i.e. no entries

            if self.options.ratio_plot:
                axes[0].errorbar(bin_centres, data_binned, 
                                 yerr=[data_binned-data_stat_down, data_stat_up-data_binned],
                                 label=sample_id, fmt='o', ms=4, color='black', capsize=0)

            else: 
                axes.errorbar(bin_centres, data_binned, 
                              yerr=[data_binned-data_stat_down, data_stat_up-data_binned],
                              label=sample_id, fmt='o', ms=4, color='black', capsize=0)

            return data_binned, bin_centres, (data_stat_down, data_stat_up)

            
    #--------------------------------------------------------------

    def plot_bkgs(self, cut_string, axes, variable, bins, data_binned, bin_centres, data_stat_down_up):
        if self.options.ratio_plot:
            ratio_canv = axes[1]
            axes = axes[0]

        if self.options.stack_bkgs:
            var_stack       = []
            weight_stack    = []
            label_stack     = []
            colour_stack    = []
            stat_down_stack = []
            stat_up_stack   = []

        for sample_id, bkg_frame in self.bkg_df_map.iteritems():
            bkg_frame_cut     = bkg_frame.query(cut_string).copy()
            var_to_plot       = bkg_frame_cut[variable.name].values
            var_weights       = bkg_frame_cut['weight'].values 

            #FIXME: store this in dict, inside systematics object
            sumw, _                   = np.histogram(var_to_plot, bins=bins, weights=var_weights)
            self.mc_totals[sample_id] = sumw
            sumw2, _                  = np.histogram(var_to_plot, bins=bins, weights=var_weights**2)
            stat_down, stat_up        = self.poisson_interval(sumw, sumw2)
            self.mc_stat_uncs[sample_id]   = [(sumw-stat_down), (stat_up-sumw)]

            print '--> Integral of hist: %s, for sample: %s, is: %.2f' % (variable.name,sample_id,np.sum(sumw))

            #FIXME: if: add options to stack backgrounds!
 

            if self.options.stack_bkgs:
                var_stack.append(var_to_plot)
                weight_stack.append(var_weights)
                label_stack.append(sample_id)
                colour_stack.append(self.colour_sample_map[sample_id])

            else: 
                if variable.norm:
                    var_weights /= np.sum(sumw)
                    stat_down   /= np.sum(sumw)
                    stat_up     /= np.sum(sumw)
                axes.hist(var_to_plot, bins=bins, label=sample_id, weights=var_weights, 
                          color=self.colour_sample_map[sample_id], histtype='stepfilled')

                #required for comparing to each systematic flucatuation

                #if we have multiple bkgs that we dont want to stack, we can't define a ratio plot i.e. Data/ which bkg? so build this in below
                if self.options.ratio_plot:
                    if len(self.bkg_df_map.keys()) == 1:
                        data_bkg_ratio   = data_binned/sumw
                        ratio_canv.errorbar( bin_centres, data_bkg_ratio, yerr=[(data_binned-data_stat_down_up[0])/sumw,(data_stat_down_up[1] - data_binned)/sumw], fmt='o', ms=4, color='black', capsize=0, zorder=3)
                    else: raise Exception("Requested a ratio plot, but have given more than one background. Did you mean to set option 'stack_backgrounds' = True? ")

            del bkg_frame_cut

        #draw stack if exists
        if self.options.stack_bkgs:
            if variable.norm:
                normed_weights   = [ normed_weights / np.sum(sumw_stack) for normed_weights in weight_stack]
                normed_stat_down = [ normed_stat_down / np.sum(sumw_stack) for normed_stat_down in stat_down_stack]
                normed_stat_up   = [ normed_stat_up / np.sum(sumw_stack) for normed_stat_down in stat_up_stack]

            axes.hist(var_stack, bins=bins, label=label_stack, weights=weight_stack, 
                            color=colour_stack, histtype='stepfilled', stacked=True)
            
            mc_stack_binned, _   = np.histogram(np.concatenate(var_stack), bins=bins, weights=np.concatenate(weight_stack))

            
            data_bkg_ratio     = data_binned/mc_stack_binned
            ratio_canv.errorbar( bin_centres, data_bkg_ratio, yerr=[(data_binned-data_stat_down_up[0])/sumw,(data_stat_down_up[1] - data_binned)/sumw], fmt='o', ms=4, color='black', capsize=0, zorder=3)


            
    #--------------------------------------------------------------

    def set_canv_style(self, axes, variable, bins):
        x_err = abs(bins[-1] - bins[-2])
        axes[0].set_ylabel('Events/%.2f' %(2*x_err) , size=14, ha='right', y=1)
        if variable.norm: axes[0].set_ylabel('1/N dN/d(%s) /%.2f' % (variable.xlabel,x_err), size=14)
        axes[0].text(0, 1.005, r'\textbf{CMS}', ha='left', va='bottom', transform=axes[0].transAxes, size=16)           
        axes[0].text(0.13, 1.005, r'\emph{%s}'%self.options.cms_label, ha='left', va='bottom', transform=axes[0].transAxes, size=14)           
        axes[0].text(1, 1.005, r'%.1f~fb$^{-1}$ %s'%(sum(self.sample_lumis), self.options.energy_label), ha='right', va='bottom', transform=axes[0].transAxes, size=14)
        #if variable.norm: axes[0].set_ylabel('1/N dN/d(%s) /%.2f' % (variable.xlabel,x_err, ha='right', y=1)
       
       
        axes[1].set_xlabel(variable.xlabel, size=14, ha='right', x=1)
        axes[1].set_ylim(self.ratio_range[0],self.ratio_range[1])
        axes[1].grid(True, linestyle='dotted')

        return axes
    

    #--------------------------------------------------------------

    def plot_systematics(self, cut_string, axes, variable, bins):

        #the systematics info can be accessed through the self.samples object which has format:

        #self.samples = { 'sample_name_1': sample_object },
        #                 'sample_name_2': sample_object } 

        #the sample object contains an "systematics" attribute, which, if we are reading a backgroud sample,
        # is a dictionary with structure:

        #sample.systematics = {'syst_1_name': syst_object, 'syst_2_name': syst_object } 


        for sample, sample_obj in self.samples.iteritems():
            if sample_obj.label == 'background':
                for syst_name, syst_obj in sample_obj.systematics.iteritems():
                    both_syst_frames = {}
                    both_syst_frames['Down'] = syst_obj.down_frame.query(cut_string)
                    both_syst_frames['Up']   = syst_obj.up_frame.query(cut_string)
                     
                    for syst_type, i_frame in both_syst_frames.iteritems():
                        var_to_plot      = i_frame[variable.name].values
                        weight           = i_frame['weight'].values
                        i_syst_binned, _ = np.histogram(var_to_plot, bins=bins, weights=weight)
                        #compare variation to the nominal for given sample and fill bin list
                        #can probably do using numpy and not loop through the arrays i.e. arr1>arr2 mask and then fill
                        true_up_variations    = []
                        true_down_variations  = []

                        #compare the systematic change to the nominal bin entries for that proc.
                        for ybin_syst, ybin_nominal in zip(i_syst_binned, self.mc_totals[sample]):
                          if ybin_syst > ybin_nominal: 
                            true_up_variations.append(ybin_syst - ybin_nominal)
                            true_down_variations.append(0)
                          elif ybin_syst < ybin_nominal:
                            true_down_variations.append(ybin_nominal - ybin_syst)
                            true_up_variations.append(0)
                          else: #sometimes in low stat cases we get no change either way wrt nominal
                            true_up_variations.append(0)
                            true_down_variations.append(0)
                          

                        #FIX: we need a new atriibute for the syst class linking the true fluctuations
                        #FIX: with the name and syst type (up/down) of the tree e.g. :
 
                        # systematic.up_tree_binned:   [true_down, true_up]
                        # systematic.down_tree_binned: [true_down, true_up]
                     
                        # because EACH up/down tree has its own true_up and true_down flucs associated to it
                        # i.e. we can't just make a single true_up/true_down inclusively for a given syst
                         
                        if syst_type=='Down':
                            syst_obj.down_tree_binned = [np.asarray(true_down_variations), 
                                                         np.asarray(true_up_variations)]
                        else:
                            syst_obj.up_tree_binned   = [np.asarray(true_down_variations), 
                                                         np.asarray(true_up_variations)]

        #add all the up/down variations (separately) for each systematic in quadrature for each bin, 
        #for each proc (not adding procs together yet though

        for sample, sample_obj in self.samples.iteritems():
            if sample_obj.label == 'background':
                down_squares = [] 
                up_squares   = [] 

                for syst_name, syst_obj in sample_obj.systematics.iteritems():
                    down_squares.append( syst_obj.down_tree_binned[0]**2 )
                    down_squares.append( syst_obj.up_tree_binned[0]**2 )

                    up_squares.append( syst_obj.down_tree_binned[1]**2 )
                    up_squares.append( syst_obj.up_tree_binned[1]**2 )

                #now add up the each bin that has been squared (will add element wise since np array)
                syst_merged_downs = np.zeros(len(bins)-1)
                syst_merged_ups   = np.zeros(len(bins)-1)

                for down_array in down_squares:
                    syst_merged_downs += down_array
                for up_array in up_squares:
                    syst_merged_ups   += up_array

                #clear individual background systematic entries and merge create new object merged across systs
                sample_obj.systematics.clear()
                merged_syst_obj = Systematic('merged_systs')

                #combined with correpsonding stat error. note that if we are considering a sample set, the name and set attributes are identical now
                #NOTE: syst have already been squared above in prep for this step!
                syst_merged_downs = np.sqrt( syst_merged_downs + self.mc_stat_uncs[sample_obj.name][0]**2) 
                syst_merged_ups   = np.sqrt( syst_merged_ups   + self.mc_stat_uncs[sample_obj.name][1]**2) 

                merged_syst_obj.merged_syst_stat_down  = syst_merged_downs
                merged_syst_obj.merged_syst_stat_up    = syst_merged_ups
                sample_obj.systematics['merged_systs'] = merged_syst_obj 

            
        #4) merge combined systs across procs id we are stacking them. else (see below)
        
        if self.options.stack_bkgs: #add uncertanties across samples and plot a single band (inc on ratio plot)
            merged_downs = np.zeros(len(bins)-1)
            merged_ups   = np.zeros(len(bins)-1)
            for sample, sample_obj in self.samples.iteritems():
                if sample_obj.label == 'background':
                    merged_downs += (sample_obj.systematics['merged_systs'].merged_syst_stat_down)**2
                    merged_ups   += (sample_obj.systematics['merged_systs'].merged_syst_stat_up)**2

            final_down = np.sqrt(merged_downs)
            final_up = np.sqrt(merged_ups)

            #do the drawing           
            if self.options.ratio_plot:
                ratio_canv  = axes[1]
                axes        = axes[0] 

            total_mc = np.zeros(len(bins)-1) # n_bins = n_bin_edges -1 
            for sumw in self.mc_totals.values():
              total_mc += sumw

            up_yield   = total_mc + final_up
            #FIXME: fix this niche issue below with poiss err function
            up_yield[up_yield==np.inf] = 0
            down_yield = total_mc - final_down

            axes.fill_between(bins, list(down_yield)+[down_yield[-1]], list(up_yield)+[up_yield[-1]], alpha=0.3, step="post", color="lightcoral", lw=0, zorder=4, label='Simulation stat. $\oplus$ syst. unc.')
            if self.options.ratio_plot:
              sigma_tot_ratio_up   = final_up/total_mc
              sigma_tot_ratio_down = final_down/total_mc
                      
              #NOTE:below is for error in ratio (i.e. including error in data propped into ratio ratio that kept sep)
              ratio_up_excess      = np.ones(len(total_mc)) + sigma_tot_ratio_up
              ratio_down_excess    = np.ones(len(total_mc)) - sigma_tot_ratio_down
                      
              #1. if we have no entries, the upper limit is inf and lower is nan
              #2. hence we set both to nan, so they aren't plot in the ratio plot
              #3  BUT if we have [nan, nan, 1 ,2 ,,, ] and/or [1, 2, ... nan, nan] 
              #   i.e. multiple Nan's at each end, then we have to set to Nan closest
              #   to the filled numbers to 1, such that error on the closest filled value
              #   doesn't mysteriously disappear
              #EDIT: gave up and did this dumb fix:

              ratio_up_excess[ratio_up_excess==np.inf] = 1 
              ratio_down_excess = np.nan_to_num(ratio_down_excess)
              ratio_down_excess[ratio_down_excess==0] =1
            
              ratio_canv.fill_between(bins, list(ratio_up_excess)+[ratio_up_excess[-1]], list(ratio_down_excess)+[ratio_down_excess[-1]], alpha=0.3, step="post", color="lightcoral", lw=1 , zorder=2)

        #5) if not stacking keep uncertainties associated to their owns sample i.e. loop through sample objects. This also means we wont drawing a ratio plot...
        #NOTE: however that if we have one proc we can still defnie a ratio plot, so take care to account for this below
        else: 
            #do the drawing           
            if self.options.ratio_plot:
                ratio_canv  = axes[1]
                axes        = axes[0] 

            for sample, sample_obj in self.samples.iteritems():
                if sample_obj.label == 'background':


                        merged_downs = sample_obj.systematics['merged_systs'].merged_syst_stat_down
                        merged_ups   = sample_obj.systematics['merged_systs'].merged_syst_stat_up

                        up_yield   = self.mc_totals[sample] + merged_ups
                        #FIXME: fix this niche issue below with poiss err function
                        up_yield[up_yield==np.inf] = 0
                        down_yield = self.mc_totals[sample] - merged_downs

                        axes.fill_between(bins, list(down_yield)+[down_yield[-1]], list(up_yield)+[up_yield[-1]], alpha=0.3, step="post", color="lightcoral", lw=0, zorder=4, label='Simulation stat. $\oplus$ syst. unc.')

                        if (self.options.ratio_plot) and (len(self.bkg_df_map.keys()) == 1):
                            total_mc             = sum(self.mc_totals.values()) #only one bkg sample so no need to loop over dict
                            sigma_tot_ratio_up   = merged_ups/total_mc
                            sigma_tot_ratio_down = merged_downs/total_mc
                                    
                            ratio_up_excess      = np.ones(len(total_mc)) + sigma_tot_ratio_up
                            ratio_down_excess    = np.ones(len(total_mc)) - sigma_tot_ratio_down
                                    
                            #1. if we have no entries, the upper limit is inf and lower is nan
                            #2. hence we set both to nan, so they aren't plot in the ratio plot
                            #3  BUT if we have [nan, nan, 1 ,2 ,,, ] and/or [1, 2, ... nan, nan] 
                            #   i.e. multiple Nan's at each end, then we have to set to Nan closest
                            #   to the filled numbers to 1, such that error on the closest filled value
                            #   doesn't mysteriously disappear
                            #EDIT: gave up and did this dumb fix:

                            ratio_up_excess[ratio_up_excess==np.inf] = 1 
                            ratio_down_excess = np.nan_to_num(ratio_down_excess)
                            ratio_down_excess[ratio_down_excess==0] =1
            
                            ratio_canv.fill_between(bins, list(ratio_up_excess)+[ratio_up_excess[-1]], list(ratio_down_excess)+[ratio_down_excess[-1]], alpha=0.3, step="post", color="lightcoral", lw=1 , zorder=2)

                        else: raise Exception("Requested a ratio plot, but have given more than one background. Did you mean to set option 'stack_backgrounds' = True? ")
                          
                          

                



    #--------------------------------------------------------------
        
    def assemble_cut_string(self, var_to_plot):
        '''Form a string of cuts to query samples. Take care to remove
           the cuts on the variable being plot
        '''
        cut_dict = self.cut_map.copy()
        if var_to_plot in cut_dict.keys(): del cut_dict[var_to_plot]
        cut_list_non_null = [cut for cut in cut_dict.values() if cut != '']
        separator = ' and '
        cut_string = separator.join(cut_list_non_null)
        return cut_string

    #--------------------------------------------------------------

    def print_cuts(self, line_length):
        cut_list_non_null = [cut for cut in self.cut_map.values() if cut != '']
        for cut in self.cut_map.values():
            if cut!='': print '|' + cut.ljust(line_length) + '|'
    #FIXME: make this into a nice format for printing

    #--------------------------------------------------------------

    def print_plotting_info(self):
        #title is 60 spaces, divider is 62 spaces
        #so if you want to add '*' or other symbols at the start/end of a line
        #use the title length 
        
        title = 'Plotting Info'.center(76)
        divider = '+' + '-'*(len(title)) + '+'


        print 
        print '*'*len(divider)
        print '*' + title + '*'

        print '-'*len(divider)

        print '|'+'Variables'.center(len(title))+'|'

        print divider 

        for var in self.var_names:
    	    print '|' + var.ljust(len(title)) + '|'

        print divider

        print '|'+'Cuts'.center(len(title))+'|'
        print divider
        self.print_cuts(len(title))

        print divider

        self.options.print_options(len(title), divider)
       
        print divider

        if len(self.syst_tree_name_map)!=0:
            print '|'+'Systematic trees'.center(len(title))+'|'
            print divider
            for syst, tree_paths in self.syst_tree_name_map.iteritems():
                print '|'+(syst+':').ljust(len(title))+'|'
                print '|'+ ('  Up:  '+ str(tree_paths[0])).ljust(len(title)) + '|'
                print '|'+ ('  Down:'+ str(tree_paths[1])).ljust(len(title)) + '|'

        print divider

        print '|'+'Samples'.center(len(title))+'|'
        print divider
        self.print_samples(len(title), divider)

        print divider
        print '*' + 'End of Plotting Info'.center(76) + '*'
        print '*'*len(divider)
       

    def print_samples(self, title_length, divider):
        bkgs = [sample_obj for sample_obj in self.samples.values() if sample_obj.label == 'background']    
        sigs = [sample_obj for sample_obj in self.samples.values() if sample_obj.label == 'signal']    
        data = [sample_obj for sample_obj in self.samples.values() if sample_obj.label == 'data']    

        print '|' + 'Signal samples:'.ljust(title_length) + '|'
        for sig in sigs:
    	    print '|' + ('  name: '+sig.name+', set: '+sig.sample_set).ljust(title_length) + '|'
        print '|' + (''.ljust(title_length)) + '|'
        print '|' + 'Background samples:'.ljust(title_length) + '|'
        for bkg in bkgs:
    	    print '|' + ('  name: '+bkg.name+', set: '+bkg.sample_set).ljust(title_length) + '|'
        print '|' + (''.ljust(title_length)) + '|'
        print '|' + 'Data samples:'.ljust(title_length) + '|'
        for dat in data:
    	    print '|' + ('  name: '+dat.name+', set: '+dat.sample_set).ljust(title_length) + '|'
    

    #--------------------------------------------------------------

    def derive_scale_factors(self, var_key, bkg_df, data_df):
        '''
        Scale MC to data in bins of a given variable and return scale factors
        '''
        variable   = self.variables[var_key]
        bins       = np.linspace( eval(variable.bin_range)[0], eval(variable.bin_range)[1],
                                 variable.num_bins)

        preselected_mc         = bkg_df.query(self.assemble_cut_string(var_key))
        preselected_data       = data_df.query(self.assemble_cut_string(var_key))
        summed_mc, _           = np.histogram(preselected_mc[var_key].values, bins=bins, 
                                              weights=preselected_mc['weight'].values)
        summed_data, bin_edges = np.histogram(preselected_data[var_key].values, bins=bins)
        scale_factors          = summed_data/summed_mc
        return  scale_factors, bin_edges

    #--------------------------------------------------------------

    def apply_reweighting(self, var_key, scale_factors, bin_edges, df):
        scaled_list = []
        for ibin in range(len(scale_factors)):
            temp_total = df[df[var_key] > bin_edges[ibin]] 
            temp_total = temp_total[temp_total[var_key] < bin_edges[ibin+1]] 
            temp_total['weight'] *= scale_factors[ibin]
            scaled_list.append(temp_total) 
        scaled_df = pd.concat(scaled_list)
        del scaled_list           
        return scaled_df


    #--------------------------------------------------------------

    def get_var_names(self):
        return self.var_names
    #--------------------------------------------------------------

    def draw_tot_error_band(self):
        pass

    #--------------------------------------------------------------
    def draw_stat_error_band(self):
        pass
        #should take syst and stat

    #--------------------------------------------------------------
    def draw_syst_error_band(self):
        pass
        #should take syst and stat

    def draw_data_mc_ratio(self):
        pass
        #should take syst, stat, boolean from var object on whether to draw or not

    def poisson_interval(self, x, variance, level=0.68):
        neff = x**2/variance
        scale = x/neff
        
        # CMS statcomm recommendation
        l = scipy.stats.gamma.interval(
            level, neff, scale=scale,
        )[0]
        u = scipy.stats.gamma.interval(
            level, neff+1, scale=scale
        )[1]
        
        # protect against no effecitve entries
        l[neff==0] = 0.
        
        # protect against no variance
        l[variance==0.] = 0.
        u[variance==0.] = np.inf
        
        return l, u

    def check_dir(self, file_path):
        if not path.isdir(file_path):
            system('mkdir -p %s'%file_path)

class trackSystError():
    '''
    Keep track on syst error in each bin, for each sample
    (could make a dict and update but more transparent as a class).
    '''
    def __init__(self):
        self.up_tree   = ''
        self.down_tree = ''


class trackStatError():
    '''
    keep track on stat error in each bin, for each sample,
    (could make a dict and update but more transparent as a class).
    '''

    def __init__(self):
        pass

