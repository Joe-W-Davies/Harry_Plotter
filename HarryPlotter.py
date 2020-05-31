import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import uproot as upr
import json


from os import system, path    
import math
                               
class Variable(object):              
    '''
    hold plotting information inside variable object
    '''

    def __init__(self, name, options={}):
        self.template = {
            'name'       : '',
            'range'      : '(0,100)',
            'num_bins'   : '30',
            'cuts'       : '',
            'lumi'       : 1.0,
            'xlabel'     : '',
            'log'        : False,
            'norm'       : False,
            'reweight'   : False, #FIXME: protect against this being defined for more than one var
            'ratio'      : False ,
            'colour'     : 'red'
        }
        self.__dict__ = self.template
        self.__dict__.update(options)
        self.name     = name
      

    def assemble_cut_string(self, varToPlot, cutMap):
        '''Form a string of cuts to query samples. Take care to remove
           the cuts on the variable being plot
        '''
        tempCuts = dict(cutMap)
        if varToPlot in tempCuts.keys(): del tempCuts[varToPlot]
        cutString = ''
        for iVar, (varName,cuts) in enumerate(tempCuts.iteritems()):
          for cut in cuts:  
            if iVar != len(tempCuts)-1:
              cutString += (varName + cut + ' and ')
            else: 
              cutString += (varName + cut)
        return cutString

    
        
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
        'sample_set' : ''
        }
        self.__dict__ = self.template
        self.__dict__.update(options)
        if self.file_path.endswith('/'): self.file_path = self.file_path[:-1]
        self.name     = name
        self.label    = self.label.lower()
        
    def check_dir(self, file_path):
        if not path.exists(file_path):
          raise ('Error: %s file path not found') %file_path

    def scale_to_lumi(self):
        '''
          Scale weights of given sample by lumi
        '''
        #apply self.lumi
        pass

    def scale_arbitrary(self):
        '''
          Scale sample weights by arbitrary ammount
        '''
        pass

    def reweight(self, reweight_var):
        '''
        scale MC to data in bins of a given variable
        '''

        bins              = np.linspace( 0, 450, 51) 
        preselected_mc    = mc_total.query(assemble_cut_string(reweight_var, cutMap))
        preselected_data  = data_total.query(assemble_cut_string(reweight_var, cutMap))
        var_mc            = preselected_mc[reweight_var].values
        weight            = preselected_mc['weight'].values
        var_data          = preselected_data[reweight_var].values
         
        summed_mc, _           = np.histogram(var_mc, bins=bins, weights=weight)
        summed_data, bin_edges = np.histogram(var_data, bins=bins)
        scale_factors          = summed_data/summed_mc
         
        scaled_list = []
        for iBin in range(len(scale_factors)):
            temp_total = mc_total[mc_total[reweight_var] > bin_edes[iBin]] 
            temp_total = temp_total[temp_total[reweight_var] < bin_edges[iBin+1]] 
            temp_total['weight'] *= scale_factors[iBin]
            scaled_list.append(temp_total) 
         
        return pd.concat(scaled_list)

class dataFrameTools(object):
    ''' 
    collection of functions to handle dataframe operations
    '''
    def __init__(self):
        pass

    #--------------------------------------------------------------
    def get_syst_trees(self):
        '''Read ROOT TTree with systematic variations'''
        pass

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

#--------------------------------------------------------------------------------------

class Options(object):
    '''Object containing the plotting options'''
    def __init__(self, options={}):
        self.template = {
            'ratio_range'  : (0,2),
            'ratio_plot'   : False,
            'total_lumi'   : False,
            'concat_years' : False,
            'var_list '    : []
        }
        self.__dict__  = self.template
        self.__dict__.update(options)

    def print_options(self):
        pass
      
#--------------------------------------------------------------------------------------

class Reader(object):
    '''Class to read the plot card'''

    def __init__(self, card, sample_dir='', output_dir=''):
        self.card          = card 
        self.variables     = {}
        self.samples       = {} #all samples in { key:df } format
        self.systematics   = {}
        self.options       = Options()
        self.ratio_range   = (0,2) #get from options
        self.total_lumi    = '137 fb$^({-1}$' #get from options
        self.cut_map       = {} #get all vars when read in
        #self.sample_dir    = sample_dir
        #self.output_dir    = output_dir

        self.legend        = ['']
        self.var_names     = []
        self.sample_sets   = set() #unique set of name that span years 
        self.sample_years  = set() 
        self.sample_labels = set() 
        self.signal_df_map = {} 
        self.bkg_df_map    = {} 
        self.data_df_map   = {} 
        self.syst_df_map   = {} 

    def read(self):
        config = json.loads(self.card)
        for key in config: 
            if 'variables' in key:
                  for var_name, var_options in config[key].iteritems():
                      var = Variable(var_name, var_options)
                      self.variables[var.name] = var
                      self.var_names.append(var.name)
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
            #for sample_set_name, sample_set_list in self.signal_df_map.iteritems():
            #    try:  assert( len(self.sample_years) == len(sample_set_list) )
            #    except: raise Exception('Number of years (%i) not equal to number of sub-samples (%i) in sample set: %s' %(len(self.sample_years), len(sample_set_list), sample_set_name) )
            #    self.sample_df_map[sample_set_name] = pd.concat(sample_set_list)

            #for sample_set_name, sample_set_list in self.bkg_df_map.iteritems():
            #    try:  assert( len(self.sample_years) == len(sample_set_list) )
            #    except: raise Exception('Number of years (%i) not equal to number of sub-samples (%i) in sample set: %s' %(len(self.sample_years), len(sample_set_list), sample_set_name) )
            #    self.sig_df_map[sample_set_name] = pd.concat(sample_set_list)

            #for sample_set_name, sample_set_list in self.data_df_map.iteritems():
            #    try:  assert( len(self.sample_years) == len(sample_set_list) )
            #    except: raise Exception('Number of years (%i) not equal to number of sub-samples (%i) in sample set: %s' %(len(self.sample_years), len(sample_set_list), sample_set_name) )
            #    self.data_df_map[sample_set_name] = pd.concat(sample_set_list)

            #refactor above code
            for df_dict in [sig_df_map, bkg_df_map, data_df_map]:
                for sample_set_name, sample_set_list in df_dict.iteritems():
                    try:  assert( len(self.sample_years) == len(sample_set_list) )
                    except: raise Exception('Number of years (%i) not equal to number of sub-samples (%i) in sample set: %s' %(len(self.sample_years), len(sample_set_list), sample_set_name) )
                    df_dict[sample_set_name] = pd.concat(sample_set_list)

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

        #if self.options.concat_years: 
        #    for sample_set in self.sample_sets:
        #        #FIXME: implement this condition in abetter way, by looking at dicts
        #        if sample_set.lower() != 'data':
        #          syst_df_map[sample_set] = {}
        #          for syst in self.systematics.keys():
        #              syst_df_map[sample_set][syst+'Up']   = []
        #              syst_df_map[sample_set][syst+'Down'] = []

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
                    up_df['weight']   *= float(sample_obj.lumi)
                    if self.options.concat_years:
                        syst_df_map[sample_obj.sample_set][syst+'Up'].append(up_df)
                    else:
                        syst_df_map[sample_obj.name][syst+'Up'] = up_df

                    down_tree          = input_file[syst_obj.down_tree]
                    down_df            = down_tree.pandas.df(self.var_names+['weight'])
                    down_df['weight'] *= float(sample_obj.lumi)
                    if self.options.concat_years:
                        syst_df_map[sample_obj.sample_set][syst+'Down'].append(down_df)
                    else:
                        syst_df_map[sample_obj.name][syst+'Down'] = down_df
        print 'before concatting'
        print syst_df_map

        if self.options.concat_years: 
            for sample_set_name, syst_dict in syst_df_map.iteritems():
                for syst_name, syst_list  in syst_dict.iteritems():
                        syst_df_map[sample_set_name][syst_name] = pd.concat(syst_list)

        self.syst_df_map = syst_df_map
        print 'after concatting'
        print self.syst_df_map

    def fill_cut_map(self):
        pass

    def print_cuts(self):
        pass

    def get_var_names(self):
        return self.var_names

class plotting():
    '''Class to perform various plot related functions'''

    def __init__(self):
        pass


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

    #--------------------------------------------------------------
    def draw(self):
        '''
        Main function to do the drawing
        '''
        #variable = self.variables[var_key]
        pass


class trackSystError():
    '''
    keep track on syst error in each bin, for each sample
    '''
    def __init__(self):
        self.up_tree   = ''
        self.down_tree = ''


class trackStatError():
    '''
    keep track on stat error in each bin, for each sample
    '''
    def __init__(self):
        pass

    def poisson_interval(x, variance, level=0.68):
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

   


