# Segment and measure areas of dense organoids in RGB images.
# by Richard Butler, Gurdon Institute Imaging Facility

import math as maths

from ij import IJ, ImagePlus, ImageStack
from ij.plugin import Duplicator
from ij.plugin.filter import MaximumFinder, ThresholdToSelection, EDM, Convolver, BackgroundSubtracter
from ij.process import ImageProcessor, FloatProcessor, ByteProcessor, ColorProcessor, Blitter, AutoThresholder, FloodFiller
from ij.measure import Calibration, ResultsTable
from ij.gui import Roi, ShapeRoi, Overlay, TextRoi

from java.awt import Color, Rectangle, Font, BorderLayout
from java.awt.image import BufferedImage

from javax.swing import JFrame

from org.jfree.chart import JFreeChart, ChartFactory, ChartPanel, ChartFrame
from org.jfree.chart.plot import PlotOrientation
from org.jfree.data.statistics import HistogramDataset, HistogramType
from org.jfree.chart.renderer.xy import StandardXYBarPainter

# object area range to include as organoids (px^2)
minA = 100
maxA = 25000


def plotHistogram(values, label):
	nbins = 32
	dataset = HistogramDataset()
	dataset.setType(HistogramType.RELATIVE_FREQUENCY)
	dataset.addSeries(label, values, nbins)
	chart = ChartFactory.createHistogram("", label, "f", dataset, PlotOrientation.VERTICAL, False, True, False)
	plot = chart.getXYPlot()
	renderer = plot.getRenderer()
	renderer.setSeriesPaint(0, Color.BLACK)
	painter = StandardXYBarPainter()
	renderer.setBarPainter(painter)
	
	chartPanel = ChartPanel(chart)
	frame = JFrame(label)
	frame.setLayout(BorderLayout())
	frame.setSize(800, 800)
	frame.setLocationRelativeTo(None)
	frame.add(chartPanel, BorderLayout.CENTER)
	frame.pack()
	frame.setVisible(True)

def getMask2D(ip, sigma):
	ip.invert()

	BackgroundSubtracter().rollingBallBackground( ip, 100, False, False, False, False, True )

	n = 11	#kernel radius
	kernel = []
	for y in range(-n,n+1):
		for x in range(-n,n+1):
			t = maths.sqrt( x**2 + y**2 )
			k = -(1.0/maths.pi*sigma**4) * abs(1-(t/(2*sigma**2))) * maths.exp(-(t/(2*sigma**2)))
			kernel.append( k )
	ip.convolve(kernel, 2*n+1, 2*n+1)

	stats = ip.getStatistics()
	reallyInt = [ int(n) for n in stats.getHistogram()]
	thresh = AutoThresholder().getThreshold(AutoThresholder.Method.Triangle, reallyInt )	#Huang
	for i in range(ip.getWidth()*ip.getHeight()):
		if ip.get(i) >= thresh:
			ip.set(i, 255)
		else:
			ip.set(i, 0)
	ip = ip.convertToByte(False)
	
	fillHoles(ip)

	watershed(ip)

	ed = int(sigma)
	for d in range(ed):
		ip.dilate()

	return ip

def watershed(ip):
	tol = 0.3	#default in EDM is 0.5
	floatEdm = EDM().makeFloatEDM(ip, 0, False)
	maxIp = MaximumFinder().findMaxima(floatEdm, tol, ImageProcessor.NO_THRESHOLD, MaximumFinder.SEGMENTED, False, True)
	if (maxIp != None):
		ip.copyBits(maxIp, 0, 0, Blitter.AND)

def fillHoles(ip):
	width = ip.getWidth()
	height = ip.getHeight()
	ff = FloodFiller(ip)
	ip.setColor(127)
	foreground = 127
	background = 0
	for y in range(height):
	    if ip.getPixel(0,y)==background:
	    	ff.fill(0, y)
	    if ip.getPixel(width-1,y)==background:
	    	ff.fill(width-1, y)
	for x in range(width):
	    if ip.getPixel(x,0)==background:
	    	ff.fill(x, 0)
	    if ip.getPixel(x,height-1)==background:
	    	ff.fill(x, height-1)
	pixels = ip.getPixels()
	n = width*height
	for i in range(n):
		if ip.get(i)==127:
		    ip.set(i, 0)
		else:
		    ip.set(i, 255)

def getRoisFromMask(ip):
	tts = ThresholdToSelection()
	ip.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE)
	composite = tts.convert(ip)
	if composite is None:
		exit("No organoids found")
	
	rois = ShapeRoi(composite).getRois()

	if len(rois)==1:	#correct location for non-composite Roi
		offset = composite.getBounds()
		rr = rois[0].getBounds()
		rois[0].setLocation(rr.x+offset.x, rr.y+offset.y)
	
	return rois

imp = IJ.getImage()
userRoi = imp.getRoi()
if userRoi is None:
	userRoi = Roi(0, 0, imp.getWidth(), imp.getHeight())
	imp.setRoi(userRoi)

ip = imp.getProcessor().convertToByte(True)

mask = getMask2D(ip, 0.75)	#sigma in px

mask.setRoi(userRoi)
mask = mask.crop()

offsetUserRoi = userRoi.clone()
offsetUserRoi.setLocation(0,0)
mask.setColor(0)
mask.fillOutside(offsetUserRoi)

rois = getRoisFromMask(mask)

rect = userRoi.getBounds()

ol = Overlay()
font = Font(Font.SANS_SERIF, Font.PLAIN, 16)
fp = FloatProcessor(rect.width, rect.height)
rt = ResultsTable()
rt.showRowNumbers(False)
rt.setPrecision(3)
areas = []
for roi in rois:
	stats = roi.getStatistics()
	if stats.area >= minA and stats.area<=maxA:

		areas.append(stats.area)

		roi.setStrokeColor(Color.MAGENTA)
		ol.add(roi)
		
		perim = roi.getLength()
		circ = 4*maths.pi*(stats.area/(perim*perim))
		
		row = rt.getCounter()
		rt.setValue("Organoid", row, row)
		rt.setValue("X", row, stats.xCentroid)
		rt.setValue("Y", row, stats.yCentroid)
		rt.setValue("Area (px"+u'\u00b2'+")", row, stats.area)
		rt.setValue("Circularity", row, circ)
		
		label = TextRoi(stats.xCentroid-6, stats.yCentroid-12, str(row), font)
		label.setStrokeColor(Color.CYAN)
		ol.add(label)

plotHistogram(areas, "Organoid Area (px"+u'\u00b2'+")")

ol.translate(rect.x, rect.y)
imp.setOverlay(ol)
rt.show(imp.getTitle()+" "+str(rect)+" Organoids")
