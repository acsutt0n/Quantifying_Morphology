# Quantifying morphology

This documentation is intended to support the paper (Sutton et al., 2016).

## Pipeline

Hoc files are the currency of quantitative morphology. Other files are common, especially neurolucida and SWC, but hoc is (1) more widespread, (2) easier for scripting programs (and humans) to read and write, and (3) interacts natively with the NEURON simulation environment. 

I have included software to convert hocs and SWCs back and forth.

0. Our neurons are dye-filled with lucifer yellow and imaged on an SP5 confocal microscope. Then, they are traced manually by skilled tracers using the freeware Knossos. (In defense of Knossos: it is free, uses XML (making it very easy to parse with scripting languages), and fast (due to excellent cubing technology developed for electron microscopy).
1. .nml/.xml files are converted to hoc files (XmlToHoc_simple), loops are removed (neuron_nxRemoveLoops) and the hoc files are scaled to um-scale (this can be done in Knossos, too) (knossos_ScaleCoords).
2. Geometry objects are created from hoc files. Ted Brookings and I designed geometry objects to behave as nested graphs to facilitate easy examination of detailed properties. During exploration, geometry objects should be manipulated interactively to get a feel for the data and how the objects respond (and their attributes). 
3. Almost all of the analysis software interacts with geometry objects. From here almost any property (neuron_getProperties) can be queried. A very few (axon identification, subtree analysis) also require some of the Matlab files.
4. Plotting: some simple plots are included in neuron_Subtrees and neuron_getProperties, but most of the plotting code is in PlottingSuite-python-, one of my repos.

* Geometry objects and some sample code are outlined in the pdf _morphology.pdf_ 

## More

Future versions will include a walkthrough.
