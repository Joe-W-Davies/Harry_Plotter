# Harry Plotter

A collection of code designed to generalise the process of taking a ROOT ntuple and producing a nice plot in pyplot. The aim is to allow the user to produce a generic configuration file to specify options for plotting.
This is still under development but will hopefully include:

* Handling of statistical variations         [DONE]

* Handling of systematic variations 
    * merging of systematics across years    [DONE]
    * combining multiple systematics for each plot [DONE]

* Merging of inputs samples across years of data taking     [DONE] 
    * (sample sets/processes specified in options)
    * (every sample within a given set is merged)

* MC Sample reweighting to data in bins of a given variable [DONE]
    * samples specified in options: can specify a a single MC and data sample, or sample set if concatting years

* Output figures generated with pyplot                      [DONE]

* Cutomisable aesthetics e.g. colour schemes, histo styles, marker sizes etc [DONE]
    * specify colour schemes in the plot card [DONE]
    * (can use own matplot stle file for more control over all other aesthetics)

* Ratio plots  [DONE]
    * for stacked backgrounds [DONE]
    * for a single background of choice [DONE]

* Stacking background distributions from different processes [DONE]
    * (utility has been added but need to double check it works properly)

* Output figures generated with ROOT


## open issues
* would be nice to synchronise the samples attribute of the Reader class such that it handles the sig and data in the same way it handles the bkgs
    * i.e. can have a sample object for each signal and data sample, and do the plotting by looping through these objects. Functionality is there for signal, but we just just rely on nested dictionaries at the moment


### other notes
* there is no reason to make a searate class for bkg MC. We can just use instances of the Sample class, and fill the {systematics} attribute with systematic objects. For data and signal, we just leave this systematic attribute empty
