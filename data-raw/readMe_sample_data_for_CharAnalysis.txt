Sample data files for CharAnalysis (May 2012):

The sample data files in this folder serve to (1) test the set up of CharAnalysis 
and (2) provide examples of theoretical records. If you cannot run these files, then
there is a problem with the way CharAnalysis is set up on the machine. The data 
are contained in both .xls and .csv format (the latter has two files associated 
with each record).

randPeaks.xls (.csv):
A random series with peaks that are proportional to background CHAR: i.e. peaks
get bigger as background gets bigger, half way through the record. The frequency
of peaks does not change throughout this record. Thus peak-detection methods 
should identify no difference in peak frequency between 1st and 2nd half of this
record. This is an extreme example, but try a global vs. local threshold on this 
record and see why a local threshold is advantageous. This record is similar to 
that displayed in Figure 2 ("Scenario 2") in:

	Higuera, P. E., D. G. Gavin, P. J. Bartlein, and D. J. Hallett. 2010. Peak 
	detection in sediment-charcoal records: impacts of alternative data
	analysis methods on fire-history interpretations. International Journal of
	Wildland Fire 19:996-1014.

redNoise.xls (.csv): 
A random series with no peaks, only red noise (random variation with 1st-order
autocorrelation of ca 0.3). Use signal-to-noise index from this series to compare 
to a real record (make sure randNoise is analyzed with the same parameters as the
real record). Any real record should have a median signal-to-noise index greater 
than that from the series of red noise. This record is similar to that displayed 
in Figure 2 ("RN") in:

	Kelly, R. F., P. E. Higuera, C. M. Barrett, and F. S. Hu. 2011. A signal-
	to-noise index to quantify the potential for peak detection in sediment-
	charcoal records. Quaternary Research 75:11-17.


COchar.xls (.csv): 
Charcoal data from Code Lake, Brooks Range, Alaska. Published in:

	Higuera, P. E., L. B. Brubaker, P. M. Anderson, F. S. Hu, and T. A. Brown.
	2009. Vegetation mediated the impacts of postglacial climate change on 
	fire regimes in the south-central Brooks Range, Alaska. Ecological 
	Monographs 79:201-219.


