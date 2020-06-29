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


