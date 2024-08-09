import random
import matplotlib.pyplot as plt
import numpy as np

def draw(size=10):
    for i in range(size):
        for j in range(size):
            plt.scatter(i, j, color='g')
    plt.show()
    return

def u_update():

    bombs = []
    for i in range(size):
        bombs.append([random.randint(3, 9), random.randint(3, 9)])

    running = True
    while running:
        # Determine State
        u_input = input('Do something: ')
        if u_input in ['q', '-q', 'quit', 'quit()']:
            running = False
            continue





    return



if __name__ == '__main__':
    u_update()

