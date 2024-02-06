from BandOfGraphene import BandOfGraphene
from DOSOfCNT import DOSOfCNT


def main():
    b = BandOfGraphene()
    b.plot()
    d = DOSOfCNT(6, 5)
    d.plot((-3, 3), direction='vertical')


if __name__ == '__main__':
    main()
