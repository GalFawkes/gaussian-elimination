import numpy as np
import sys


def orderer(a, b, u):
    pointer = 0
    for i in range(0, a.shape[0]):
        if a[i][0] > a[pointer][0]:
            pointer = i
    a[[0, pointer]] = a[[pointer, 0]]
    b[[0, pointer]] = b[[pointer, 0]]
    u[[0, pointer]] = u[[pointer, 0]]


def solver():
    print('this is a solver for Ax=b via Gaussian elimination.\nIt currently only accepts square matrices.')
    print('please enter rows as numbers separated by commas, no spaces.')
    nRows = int(input("How many rows?\n"))
    A = []
    U = []
    for i in range(0, nRows):
        row = (input('please enter row {} of A\n'.format(i+1)))
        row = row.split(',')
        while len(row) != nRows:
            print('invalid input, please re-enter row.')
            row = (input('please enter row {} of A\n'.format(i + 1)))
            row = row.split(',')
        A.append(row)
    A = np.array(A)
    U = np.identity(nRows)
    # print('A = {}'.format(A))
    b_string = input('please enter b\n')
    b = b_string.split(',')
    while len(b) != nRows:
        print('invalid b, please try again.')
        b_string = input('please enter b\n')
        b = b_string.split(',')
    b = np.array(b)
    # print('b: {}'.format(b))
    orderer(A, b, U)
    print('{}x = {}'.format(A, b))


if __name__ == '__main__':
    solver()
