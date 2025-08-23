import imageio.v3 as iio
import laserbeamsize as lbs
import numpy as np
from matplotlib import pyplot as plt

file = "2_pos.tiff"

beam = iio.imread(file)

x, y, d_major, d_minor, phi = lbs.beam_size(beam)
print("The center of the beam ellipse is at (%.0f, %.0f)" % (x, y))
print("The major axis (diameter) is %.0f pixels" % d_major)
print("The minor axis (diameter) is %.0f pixels" % d_minor)
print("The major axis is rotated %.0f° CCW from the horizontal" % (phi * 180/3.1416))

lbs.plot_image_analysis(beam)
plt.show()