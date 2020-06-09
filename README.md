# Harry Plotter

A collection of code designed to generalise the process of taking a ROOT ntuple and producing a nice plot in pyplot. The aim is to allow the user to produce a generic configuration file to specify options for plotting.
This is still all under development but will hopefully include:

* Handling of statistical and systematic variations         [DONE]
* Merging of inputs samples across years of data taking     [DONE] 
    * sample sets specified in options. Every sample within a given set is merged
* MC Sample reweighting to data in bins of a given variable [DONE]
    * samples specified in options: can specify a a single MC and data sample, or sample set if concatting years
* Output figures generated with pyplot                      [DONE]
* Cutomisable aesthetics e.g. colour schemes, histo styles, marker sizes etc [DONE]
    * specify colour schemes in the plot cars
    * can use own matplot stle file for more control over all other aesthetics
* Ratio plots (including stacked histograms)
* Output figures generated with ROOT

