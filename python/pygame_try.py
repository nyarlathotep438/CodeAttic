# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 21:54:18 2025

@author: 96328
"""
import pygame
import sys
 
# Initialize pygame
pygame.init()
 
# Define colors
BLACK = (0, 0, 0)
 
# Set the width and height of the table and its lines
TABLE_WIDTH = 400
TABLE_HEIGHT = 400
LINE_THICKNESS = 2
 
# Calculate the size of each cell
CELL_WIDTH = TABLE_WIDTH // 5
CELL_HEIGHT = TABLE_HEIGHT // 5
 
# Set up the display
screen = pygame.display.set_mode((TABLE_WIDTH, TABLE_HEIGHT))
pygame.display.set_caption("5x5 Table")
 
# Game loop
running = True
while running:
    # Handle events
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
 
    # Clear the screen
    screen.fill((255, 255, 255))  # Fill with white or any other background color
 
    # Draw the table lines
    for i in range(6):  # 6 horizontal lines for a 5x5 table
        pygame.draw.line(screen, BLACK, (0, i * CELL_HEIGHT), (TABLE_WIDTH, i * CELL_HEIGHT), LINE_THICKNESS)
    for i in range(6):  # 6 vertical lines for a 5x5 table
        pygame.draw.line(screen, BLACK, (i * CELL_WIDTH, 0), (i * CELL_WIDTH, TABLE_HEIGHT), LINE_THICKNESS)
 
    # Update the display
    pygame.display.flip()
 
# Quit pygame
pygame.quit()
sys.exit()
