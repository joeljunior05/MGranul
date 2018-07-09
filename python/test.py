import msgranul
import numpy as np
import matplotlib.pylab as plt

im = plt.imread("demo.png")

print ( dir(msgranul))

def to_grayscale(im, weights = np.c_[0.2989, 0.5870, 0.1140, 0]):
    """
    Transforms a colour image to a greyscale image by
    taking the mean of the RGB values, weighted
    by the matrix weights
    """
    tile = np.tile(weights, reps=(im.shape[0],im.shape[1],1))
    return np.sum(tile * im, axis=2)

im = to_grayscale(im)

kernels = msgranul.createkernels('c', 60, 60);
print (kernels)
locals = msgranul.correlate(im, kernels);