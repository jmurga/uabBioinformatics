import allel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def variantsByWindows(variants, windowSize, density=True, title=None):
	'''
	Function to plot the total number of variants along the VCF.
	
	Parameters
	----------
	variants: variable containing total number of variants at the VCF file
	windowsSize: step to plot
	density: numberOfVariants/windowsSize
	title: string to plot a title

	Returns
	-------
	Plot SNPs by windows depending on density
	'''

	# setup windows 
	bins = np.arange(min(variants), max(variants), windowSize)
	
	# use window midpoints as x coordinate
	x = (bins[1:] + bins[:-1])/2
	
	# compute variant density in each window
	h, _ = np.histogram(variants, bins=bins)
	
	if(density == True):
		y = h / windowSize
	else:
		y = h
	
	# plot
	fig, ax = plt.subplots(figsize=(16, 5))
	ax.plot(x, y)
	ax.set_xlabel('Chromosome position (bp)')
	ax.set_ylabel('Variant density (bp$^{-1}$)')
	if title:
		ax.set_title(title)

def plotVariantHist(vcfInput,attribute, bins=30):
	'''
	Function to make an histogram of VCF attributes.
	
	Parameters
	----------
	vcfInput: variable containing raw VCF information
	attribute: str field name to plot 
	bins: number of bins to group the histogram
	
	Returns
	----------
	Histogram of attributes
	'''

	x = vcfInput['variants/' + attribute]
	fig, ax = plt.subplots(figsize=(7, 5))
	
	ax.hist(x, bins=bins)
	
	ax.set_xlabel(attribute)
	ax.set_ylabel('No. variants')
	ax.set_title('Variant %s distribution' % attribute)

def statisticsByWindows(positions,statisticValues,statisticName=None,title=None):
	'''
 	Function to make an line chart representing statistics values by windows
	
	Parameters
	----------
	positions: ndarray containing windows analyzed by scikit-allel
	statisticsValues: ndarray containing values calculated on genomic regions
	windowsSize: windows size analyzed 
	statiscticName: string to plot the statistic name analyzed
	title: string to plot a title
	Returns
	-------
	Line Chart
	'''
	
	# Define x axis through all positions
	x = positions[:,0]
	
	# Define y axis with statistics estimation
	y = statisticValues

	# plot
	fig, ax = plt.subplots(figsize=(16, 5))
	ax.plot(x, y)
	ax.set_xlabel('Chromosome position (bp)')

	if title:
		ax.set_title(title)
	if statisticName:
		ax.set_ylabel(statisticName)
		
def plotLd(gn, title):
	'''
	Function to estimate linkage disequilibrium parameter r for each pair of variants using the method of Rogers and Huff (2008) and plot a matrix of genotype linkage disequilibrium values between all pairs of variants.
	
	Parameters
	----------
	gn: genotype array
	title: string to plot a title
	
	Returns
	-------
	Scatter plot showing LD between all pairs of variants
	'''
	m = allel.rogers_huff_r(gn) ** 2
	ax = allel.plot_pairwise_ld(m)
	ax.set_title(title)

def ldPrune(gn, size, step, threshold=.1, n_iter=1):
	'''
	Function to remove correlated SNPs. Locate SNPs that are not correlated with each other, using the locate_unlinked() function from scikit-allel. Works by sliding a window along the data, computing pairwise LD between all SNPs within each window, then removing one SNP from each correlated pair.
	
	Parameters
	----------
	gn: genotype array
	size: string to plot a title
	size : int, optional
		The window size (number of bases).
	step : int, optional
		The distance between start positions of windows. If not given,
		defaults to the window size, i.e., non-overlapping windows.
	Returns
	-------
	Genotype variable without correlated SNPs
	'''
	for i in range(n_iter):
		loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
		n = np.count_nonzero(loc_unlinked)
		n_remove = gn.shape[0] - n
		print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
		gn = gn.compress(loc_unlinked, axis=0)
	return(gn)

def plotPcaCoords(coords, model, pc1, pc2, ax, sample_population):

	pop_colours = {
    	'CEU':'#2d74b2',
    	'CHB':'#33b033',
    	'YRI':'#f5f39f'
	}

	x = coords[:, pc1]
	y = coords[:, pc2]
	for pop in sample_population:
		flt = (sample_population == pop)
		ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], 
				label=pop, markersize=6, mec='k', mew=.5)
	ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
	ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))

def figPca(coords, model, title, sample_population=None):
	'''
	Parameter
	-------
    coords : ndarray, float, shape (n_samples, n_components)
        Transformed coordinates for the samples.
    model : GenotypePCA
        Model instance containing the variance ratio explained and the stored
        components (a.k.a., loadings). Can be used to project further data
        into the same principal components space via the transform() method.
	title : string to plot a title
	sample_populations : ndarray of populations analyzed

	Return
	-------
	PCA plot
	'''
	# plot coords for PCs 1 vs 2, 3 vs 4
	fig = plt.figure(figsize=(10, 5))
	ax = fig.add_subplot(1, 2, 1)
	plotPcaCoords(coords, model, 0, 1, ax, sample_population)
	ax = fig.add_subplot(1, 2, 2)
	plotPcaCoords(coords, model, 2, 3, ax, sample_population)

	fig.suptitle(title, y=1.02)
	fig.tight_layout()  

def ihsPlot(score,sites,ylabel,absolute=False):
	'''
	Function to make and scatter plot of iHS values
	
	Parameter
	-------
	score : ndarray of iHS value (unstandarized or standarized)
	sites: ndarray of analyzed sites

	Return
	-------
	Scatter plot
	'''
	if(absolute is True):
		score = np.abs(score)
		
	plt.figure(figsize=(16, 4))
	plt.plot(sites, score, linestyle=' ', marker='o', mfc='none')
	plt.grid(axis='y')
	plt.xlabel('Position (bp)')
	plt.ylabel(ylabel)
	plt.autoscale(axis='x', tight=True)
