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

    def __repr__(): pass

#--------------------------------------------------------------------------------------
    
        
class Sample(object):
    '''
    object for sample attributes
    '''

    def __init__(self, name, options={}):
        self.template={
        'file_path'  : '',
        'tree'       : '',
        'order'      : None,
        'label'      : '',
        'lumi'       : 1.0,
        'scale'      : 1.0,
        'year'       : None,
        'file_ext'   : '.root',
        'sample_set' : '',
        'colour'     : 'red'
        }
        self.__dict__ = self.template
        self.__dict__.update(options)
        if self.file_path.endswith('/'): self.file_path = self.file_path[:-1]
        self.name     = name
        self.label    = self.label.lower()
        
    def __repr__(): pass



#--------------------------------------------------------------------------------------

class Systematic(object):
    '''Object containing the systematic options'''
    def __init__(self, name, options={}):
        self.template = {
            'up_tree'      : None,
            'down_tree'    : None,
            'up_histo'     : [],
            'down_histo'   : [],
            'true_up'      : {},
            'true_down '   : {},
            'names '       : []
        }
        self.__dict__  = self.template
        self.__dict__.update(options)
        self.name = name

    def __repr__(): pass

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
            'stack_bkgs'         : False,
            'var_list '          : []

        }
        self.__dict__  = self.template
        self.__dict__.update(options)
        self.total_lumi          = 137 #FIXME: update based on the set of lumis that are provided in json config

    def __repr__(): pass

    def print_options(self):
        pass
      
#--------------------------------------------------------------------------------------

class Plotter(object):
    '''Class to read the plot card'''

    def __init__(self, card, output_dir=''):
        self.card                = card 
        self.variables           = {}
        self.samples             = {} #all samples in { key:df } format
        self.systematics         = {}
        self.options             = Options()
        self.ratio_range         = (0,2) #FIXME: could make changeable for each variable? (add att to var)
        self.cut_map             = {} 
        self.output_dir          = output_dir
        self.legend              = []
        self.var_names           = []
        self.sample_sets         = set() #unique set of name that span years 
        self.sample_years        = set() 
        self.sample_labels       = set() 
        self.signal_df_map       = {} 
        self.bkg_df_map          = {} 
        self.data_df_map         = {} 
        self.syst_df_map         = {} 
        self.colour_sample_map   = {} 
        self.rew_scale_factors   = []
        self.rew_bin_edges       = np.array([])

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
                    syst = Systematic(syst_name, syst_options)
                    self.systematics[syst.name] = syst
            if 'samples' in key:
                samples = {}
                for sample_name, sample_opts in config[key].iteritems():
                    sample = Sample(sample_name, sample_opts)
                    samples[sample_name]  = sample
                    self.sample_sets.add(sample.sample_set)
                    self.sample_years.add(sample.year)
                    self.colour_sample_map[sample.sample_set] = sample.colour
                self.samples = samples

    def trees_to_dfs(self):
        """
        Function to read in files for each sample and convert to a DataFrame.
        If concat across years, hold dataframes in dict that looks like this: 
        {'sample_set_1' : [df_for_year_1, df_for_year_2, ... ], 'sample_set_2' : [...], ...}
        The lists will then be concattenated across years i.e. flattened out
        Otherwise just store with one key per sample:
        {'sample_id_1' : df_1, 'sample_id_2' : ..., ...}
        #Dictionaries are split into signal, background, and Data
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
        Only reads systematics for processes labelled with bkg.
        If concat across years, hold dataframes in dict that looks like this: 
        {'sample_set_1' : {'syst_1_up': [syst_df_up_for_year_1, syst_df_up_for_year_2, ... ],
                           'syst_1_down': [syst_df_down_for_year_1, syst_df_down_for_year_2, ... ]
                           'syst_2': [...], ... }, 'sample_set_2' : {...}, ...}
        The lists  will then be concattenated across years i.e. flattened out
        Otherwise just store with one key per sample:
        {'sample_id_1' : {syst_1_up : df_up, syst_1_down : df_down, syst_2_up: ..., ...}, 'sample_id_2' : ..., ...}
        """

        syst_df_map = {}

        if self.options.concat_years: 
            for sample_obj in self.samples.values():
                if sample_obj.label=='background':
                  if sample_obj.sample_set not in syst_df_map.keys(): syst_df_map[sample_obj.sample_set] = {}
                  for syst in self.systematics.keys():
                      syst_df_map[sample_obj.sample_set][syst+'Up']   = []
                      syst_df_map[sample_obj.sample_set][syst+'Down'] = []
        else: #use the actual sample key rather than the key for the set as in above
            for sample_obj in self.samples.values():
                if sample_obj.label=='background':
                    for syst in self.systematics.keys():
                        syst_df_map[sample_obj.name][syst+'Up']   = pd.DataFrame()
                        syst_df_map[sample_obj.name][syst+'Down'] = pd.DataFrame()

        for sample_name, sample_obj in self.samples.iteritems():
            if sample_obj.label=='background':
                input_file = upr.open( '%s/%s%s' % (sample_obj.file_path, sample_obj.name, sample_obj.file_ext) )
                for syst, syst_obj in self.systematics.iteritems():
                    print 'reading systematic %s for sample %s' %(sample_name, syst)
                    up_tree            = input_file[syst_obj.up_tree]
                    up_df              = up_tree.pandas.df(self.var_names+['weight'])

                    up_df['weight']   *= (float(sample_obj.lumi) * float(sample_obj.scale))
                    if len(self.options.reweight_var)!=0: 
                        up_df = self.apply_reweighting(self.options.reweight_var,
                                                       self.rew_scale_factors, 
                                                       self.rew_bin_edges, 
                                                       up_df)

                    if self.options.concat_years:
                        syst_df_map[sample_obj.sample_set][syst+'Up'].append(up_df)
                    else:
                        syst_df_map[sample_obj.name][syst+'Up'] = up_df


                    down_tree          = input_file[syst_obj.down_tree]
                    down_df            = down_tree.pandas.df(self.var_names+['weight'])
                    down_df['weight'] *= (float(sample_obj.lumi) * float(sample_obj.scale))
                    if len(self.options.reweight_var)!=0: 
                        down_df = self.apply_reweighting(self.options.reweight_var,
                                                         self.rew_scale_factors,
                                                         self.rew_bin_edges,
                                                         down_df)

                    if self.options.concat_years:
                        syst_df_map[sample_obj.sample_set][syst+'Down'].append(down_df)
                    else:
                        syst_df_map[sample_obj.name][syst+'Down'] = down_df

        if self.options.concat_years: 
            for sample_set_name, syst_dict in syst_df_map.iteritems():
                for syst_name, syst_list  in syst_dict.iteritems():
                    syst_df_map[sample_set_name][syst_name] = pd.concat(syst_list)

        self.syst_df_map = syst_df_map

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

        if len( self.syst_df_map) != 0: 
            self.plot_systematics(cut_string, axes, variable, bins)

        #default axes dont have enough whitespace above plots, so add more. also set logy here
        # (not in hists where you will get weird behaviour if y~0)

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

            sumw, _           = np.histogram(var_to_plot, bins=bins, weights=var_weights)
            sumw2, _          = np.histogram(var_to_plot, bins=bins, weights=var_weights**2)
            #FIXME: store this in dict, inside systematics object
            stat_down, stat_up = self.poisson_interval(sumw, sumw2)
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
        axes[0].text(1, 1.005, r'%.0f~fb$^{-1}$ (13 TeV)'%self.options.total_lumi, ha='right', va='bottom', transform=axes[0].transAxes, size=14)
        #if variable.norm: axes[0].set_ylabel('1/N dN/d(%s) /%.2f' % (variable.xlabel,x_err, ha='right', y=1)
       
       
        axes[1].set_xlabel(variable.xlabel, size=14, ha='right', x=1)
        axes[1].set_ylim(self.ratio_range[0],self.ratio_range[1])
        axes[1].grid(True, linestyle='dotted')

        return axes
    

    #--------------------------------------------------------------

    def plot_systematics(self, cut_string, axes, variable, bins):
        if self.options.ratio_plot:
            ratio_panel = axes[1]
            axes = axes[0] 

        for sample_id, syst_dict in self.syst_df_map.iteritems():
            for syst_name, syst_df  in syst_dict.iteritems():
                syst_df_cut = syst_df.query(cut_string)

        #FIXME: use multiple functions to handle systemtic variations, 
        #FIXME  and append results to attributes of a systematics class!

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

    def print_cuts(self):
        print self.cut_map
        #FIXME: make this into a nice format for printing

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
        #should take syst and stat

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


   


