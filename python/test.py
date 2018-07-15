from module import msgranul
import numpy as np
from skimage.io import imread

im = imread("demo.png", as_grey=True)

print (msgranul)

kernels = msgranul.createkernels('c', 60, 60);
print (kernels)

locals = msgranul.correlate(im, kernels);
print (locals)