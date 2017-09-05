# TODO

* Put manuscript on overleaf (DONE)
* Update shiny app with normal plots (AL)
* Try plotly (AL)
* Read/comment manuscript (AL)
* Modify figure labels on summary (AL)
* Get clinical analysis updates from Irina (ER)

* Change tumor/normal to lower case in *_met (DONE)
* Send manuscript highlights email (DONE)
* Do inverse (DONE)
* Look through uniquemets (DONE)
* Do string comparison to KEGG removing parentheses content (DONE)

# Methods Notes

* For four studies (BLCA,LGG,PRAD,PAAD), normalized data in linear (i.e. not log-scale) units was not available. Since we have little information on the details of sample preparation and analysis for these samples, our approach to including these data was to apply the minimal transformation required to make the data comparable to other data in our study, and analyzable using non-parametric (i.e. rank-based) approaches.

* 1. In two studies (BLCA and PAAD), the data was normalized by the authors in the original publication, but reported in log-scaled units. In BLCA, data was normalized after log-transformation, and in PAAD, data was normalized prior to log-transformation. For both studies, we re-exponentiated all the data (i.e. returned it to linear scale from log scale). This made the data amenable to plotting and other visualization while retaining the rank-order of all measurements.

* 2. For two studies (LGG and PRAD), the data was only reported as raw ion counts. For a given metabolite in a given study, the median of all measured values of that metabolite was calculated. Next, all missing data for a given metabolite was imputed to the lowest measured value of that metabolite. Finally, all measurements for a given metabolite were  normalized by the previously calculated median measurement for that metabolite.

* In BLCA, sample P14N.EU is labeled as CIS, presumably carcinoma in situ. This label is ambiguous, as the suffix "N" generally indicates a normal sample. However, given that it is CIS, we have treated it as a tumor for the purposes of our analysis. We have also removed the P14 pair from any paired analysis.

* In PRAD, there are two measurements of 2-hydroxybutyrate, as quantified by either LC or GC. We retained only the LC measurement.

* Add comments about the transformation of prostrate grades
* ~~Example completed task~~
