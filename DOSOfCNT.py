import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from tqdm import tqdm
from BandOfGraphene import BandOfGraphene

plt.rcParams["font.size"] = 30
plt.rcParams["font.family"] = 'Times New Roman'
rc('mathtext', **{'rm': 'serif', 'bf': 'serif:bold', 'fontset': 'cm'})
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.width'] = 2.0
plt.rcParams['ytick.major.width'] = 2.0
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['figure.subplot.bottom'] = 0.2
plt.rcParams['figure.subplot.left'] = 0.2


def delta(x, n_delta=10000):  # デルタ関数はガウス関数で近似. nを大きくするほどデルタ関数に近づく
    ret = np.exp(-1 * n_delta * x ** 2)
    return np.sqrt(n_delta / np.pi) * ret


class DOSOfCNT:
    a = 2.46
    a1 = np.array([np.sqrt(3) / 2, 1 / 2]) * a
    a2 = np.array([np.sqrt(3) / 2, - 1 / 2]) * a

    def __init__(self, n: int, m: int):
        self.n = n
        self.m = m
        # カイラルベクトル
        self.Ch = n * self.a1 + m * self.a2
        self.L = np.linalg.norm(self.Ch)
        # 直径
        self.d_tube = self.L / np.pi
        gcd = np.gcd(2 * m + n, 2 * n + m)
        self.t1 = (2 * m + n) / gcd
        self.t2 = -1 * (2 * n + m) / gcd
        # 並進ベクトル
        self.T = self.t1 * self.a1 + self.t2 * self.a2
        # 単位胞内の六角形の数
        self.N = int(2 * (n ** 2 + m ** 2 + n * m) / gcd)

        # 逆格子ベクトル
        self.b1 = np.array([1 / np.sqrt(3), 1]) * 2 * np.pi / self.a
        self.b2 = np.array([1 / np.sqrt(3), -1]) * 2 * np.pi / self.a
        self.K1 = 1 / self.N * (- self.t2 * self.b1 + self.t1 * self.b2)
        self.K2 = 1 / self.N * (m * self.b1 - n * self.b2)
        # 積分の範囲
        k_min = -1 * np.pi / np.linalg.norm(self.T)
        k_max = np.pi / np.linalg.norm(self.T)
        self.dk = 1 / 1000
        self.k_range = np.linspace(k_min, k_max, int(1 / self.dk))

        self.band_of_graphene = BandOfGraphene()
        self.band_of_graphene.s = 0  # 結合性バンド，半結合性バンドが対称になる

    def E_mu(self, k, mu, sign):  # グラフェンのエネルギー分散関係に，CNTの量子化条件を入れる
        x, y = k * self.K2 / np.linalg.norm(self.K2) + mu * self.K1
        return self.band_of_graphene.E_2g(x, y, sign)

    def DOS(self, E):
        ret = 0
        for k in tqdm(self.k_range):
            for mu in range(1, self.N + 1):
                ret += self.dk * (delta(E - self.E_mu(k, mu, '+')) + delta(E - self.E_mu(k, mu, '-')))  # 結合性バンドと反結合性バンドのDOSの和
        return ret / (2 * np.pi) / 2

    def plot(self, E_range: tuple, direction: str = 'vertical'):
        if len(E_range) != 2:
            print('E_range should be tuple of 2 elements.')
            return -1
        if E_range[0] > E_range[1]:
            E_range = (E_range[1], E_range[0])
        if direction not in ['vertical', 'horizontal']:
            print('wrong direction. direction should be "vertical" or "horizontal".')
            return -1

        plt.figure(figsize=(8, 10))
        E = np.linspace(*E_range, 1000)
        if direction == 'vertical':
            plt.plot(E, self.DOS(E), color='k')
            plt.xlabel('Energy [eV]')
            plt.ylabel('DOS')
            plt.xlim(*E_range)
            plt.ylim(0, )
            plt.xticks(range(E_range[0], E_range[1] + 1))
            plt.yticks([0])
        elif direction == 'horizontal':
            plt.plot(self.DOS(E), E, color='k')
            plt.xlabel('DOS')
            plt.ylabel('Energy [eV]')
            plt.xlim(0, )
            plt.ylim(*E_range)
            plt.xticks([])
            plt.yticks(range(E_range[0], E_range[1] + 1))
        else:
            print('wrong direction. direction should be "vertical" or "horizontal".')
            return 0
        plt.show()


def main():
    d = DOSOfCNT(6, 5)
    d.plot((-3, 3), direction='vertical')


if __name__ == '__main__':
    main()
