# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 13:26:22 2024

@author: 96328
"""
import matplotlib.pyplot as plt

img0 = plt.imread('C:/Users/96328/Image 7.tif')
img1 = plt.imread('C:/Users/96328/gray_image.jpg')

fig,axs = plt.subplots(1, 2)

axs[0].imshow(img0,'grey')
axs[0].set_title('original')
axs[0].axis('off')

axs[1].imshow(img1,'grey')
axs[1].set_title('grey')
axs[1].axis('off')