from module import msgranul
import numpy as np
from skimage.io import imread, imsave

im = imread("demo.jpg")

print (msgranul)

kernels = msgranul.createkernels('c', 15, 15)
print ( "\n kernels: " + str( len(kernels)))


locals = msgranul.correlate(im, kernels)
print ( "\n maxlocals before apply granulometry: " + str( len(locals)))


locals = msgranul.apply(im, locals, 0.2)[0]
print ( "\n maxlocals after apply granulometry: " + str( len(locals)))


out = msgranul.printmaxlocals(im, locals)
imsave("output_test_image.jpg", out)